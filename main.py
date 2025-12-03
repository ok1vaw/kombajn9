
# ---------------------------------------------------------
# main_v9_0_from_8_0.py — MONITOR HLADINY / ALARMOVÝ SYSTÉM
# Sjednocená detekce A/U (verze 8.4) pro online i offline režim
#
# Vygenerováno: 01.12.2025 11:00:44 CET (UTC+1), tj. 10:00:44 UTC
# ---------------------------------------------------------

import os
import csv
import time
import math
from datetime import datetime, timedelta

# pokus o import optional knihoven
try:
    import requests
except Exception:
    requests = None

try:
    from serial import Serial, SerialException
    _HAS_SERIAL = True
except Exception:
    Serial = None
    SerialException = Exception
    _HAS_SERIAL = False

# ============================================================= bum
# KONFIGURACE
# =============================================================
BASE_DIR = r"D:\data\Vojta\jimka"       # adresář pro všechny výstupy
os.makedirs(BASE_DIR, exist_ok=True)

MERENI_FILENAME_BASE   = "mereni_jimka"
CYKLY_FILENAME_BASE    = "cykly_cerpadla"
ALARMY_TXT_BASE        = "alarmy"
PRUSECIKY_FILENAME_BASE = "pruseciky"   # případně později pro vizualizaci
SMARTWELL_URL          = "https://mereni.smartwell.cz/online/00200034?typ=2"

# Sampling / thresholds
SAMPLE_INTERVAL    = 30                  # s
PUMP_LIMIT_CM      = 90
FORCED_A_LIMIT = 89.0
SMS_ALARM_LIMIT    = 165                 # cm
SMS_ALARM_COUNT    = 6                   # aktivace: 6 posledních bodů nad limitem
SMS_CLEAR_COUNT    = 2                   # deaktivace: 2 poslední body pod limitem

# A/U detection parametry
REG_WIN               = 6                # délka regresního okna
ANGLE_THRESHOLD_DEG   = 15.0             # minimální |úhel| pro kandidáta A/U
BUFFER_LIMIT_ONLINE   = 2000             # max počet bodů pro online buffer

# SMS config
SMS_PHONE = "+420723537207"

# DEBUG
DEBUG = True

# OFFLINE přepínač (kvůli SMS apod.)
OFFLINE_FLAG = False

# =============================================================
# HELPER FUNCTIONS
# =============================================================
def now_cet_str():
    return datetime.now().strftime("%d.%m.%Y %H:%M:%S")


def dprint(msg):
    if DEBUG:
        print(msg)


def cm_to_liters(h):
    """
    Převod hladiny [cm] na objem [l].
    Funkce převzatá z původní verze 8.0.
    """
    try:
        h = float(h)
    except Exception:
        return 0.0

    if h <= 125:
        return 5.140906685 * h - 254.8316105 + 152.0134768
    elif h <= 140:
        return 0.1029655281 * h**2 - 12.78365812 * h + 376.9026136 + 152.0134768
    else:
        return -0.00001540016 * h**3 - 0.0801546650 * h**2 + 58.24673706 * h - 5935.938885 + 152.0134768


# =============================================================
# File utilities
# =============================================================
def get_output_path(base, by_date=None, offline=False, ext=".csv"):
    """
    base      – základ názvu souboru (např. "cykly_cerpadla")
    by_date   – datum pro vložení do názvu, default = dnešní
    offline   – pokud True, přidá suffix "_OFFLINE"
    """
    base_dir = os.path.abspath(os.path.expanduser(BASE_DIR))
    os.makedirs(base_dir, exist_ok=True)

    if offline:
        # offline výsledky do jednoho souboru
        return os.path.join(base_dir, f"{base}_OFFLINE{ext}")

    if by_date is None:
        date_part = datetime.now().strftime("%d%m%y")
    else:
        if isinstance(by_date, datetime):
            date_part = by_date.strftime("%d%m%y")
        else:
            date_part = str(by_date)

    return os.path.join(base_dir, f"{base}_{date_part}{ext}")


def safe_append_csv(path, row=None, header=None, retries=8, delay=0.5):
    """
    Bezpečný append do CSV s několika pokusy (kvůli otevřeným souborům v Excelu).
    """
    for attempt in range(retries):
        try:
            exists = os.path.exists(path)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "a", newline="", encoding="utf-8") as f:
                w = csv.writer(f, delimiter=";")
                if header and not exists:
                    w.writerow(header)
                if row is not None:
                    w.writerow(row)
            return True
        except PermissionError:
            if attempt == 0:
                dprint(f"[WARN] Soubor {os.path.basename(path)} otevřen – čekám...")
            time.sleep(delay)
        except Exception as e:
            dprint(f"[ERROR] Zápis do {path} selhal: {e}")
            return False

    dprint(f"[WARN] Nepodařilo se zapsat do {os.path.basename(path)} po {retries} pokusech.")
    return False


# =============================================================
# MODEM / SMS
# =============================================================
def sanitize_message(msg: str) -> str:
    try:
        safe = msg.encode('ascii', errors='ignore').decode()
        safe = safe.replace('–', '-').replace('—', '-').replace('’', "'")
        return safe.strip()
    except Exception:
        return ''.join(ch for ch in msg if ord(ch) < 128)


def find_modem(force_port=None, ports_range=range(1, 21)):
    if not _HAS_SERIAL:
        return None

    if force_port:
        try:
            ser = Serial(force_port, 115200, timeout=1)
            ser.write(b'AT\r')
            time.sleep(0.3)
            resp = ser.read(ser.in_waiting).decode(errors='ignore')
            ser.close()
            if 'OK' in resp:
                return force_port
        except Exception:
            pass

    for i in ports_range:
        port = f"COM{i}"
        try:
            ser = Serial(port, 115200, timeout=1)
            ser.write(b'AT\r')
            time.sleep(0.3)
            resp = ser.read(ser.in_waiting).decode(errors='ignore')
            ser.close()
            if 'OK' in resp:
                return port
        except Exception:
            continue

    return None


def send_sms_robust(phone_number: str, message: str, force_port=None, max_retries=3):
    if not _HAS_SERIAL:
        dprint("[SMS] pyserial not available; SMS disabled.")
        return False, "NO_SERIAL"

    port = find_modem(force_port)
    if not port:
        dprint("[SMS] Modem not found.")
        return False, "MODEM_NOT_FOUND"

    msg = sanitize_message(message)

    for attempt in range(1, max_retries + 1):
        try:
            dprint(f"[SMS] Pokus {attempt}/{max_retries} -> {phone_number}: '{msg}'")
            ser = Serial(port, 115200, timeout=5)
            ser.write(b'AT+CFUN=1\r')
            time.sleep(1.2)
            ser.write(b'AT+CMGF=1\r')
            time.sleep(0.6)
            ser.write(f'AT+CMGS="{phone_number}"\r'.encode())
            time.sleep(0.6)
            ser.write((msg + "\x1A").encode())
            time.sleep(6)
            resp = ser.read(ser.in_waiting).decode(errors='ignore')
            ser.close()

            if ('+CMGS' in resp) or ('OK' in resp and 'ERROR' not in resp):
                dprint("[SMS] Modem potvrdil odeslání.")
                return True, "SUCCESS"
            else:
                dprint("[SMS] Selhání, odpověď modemu: " + repr(resp.strip()))
                time.sleep(1.0)
        except Exception as e:
            dprint(f"[SMS] Výjimka při odesílání: {e}")
            try:
                ser.close()
            except Exception:
                pass
            time.sleep(1.0)

    return False, "SMS_SEND_FAILED"


# =============================================================
# A/U DETEKCE — verze 8.4 (sjednocená pro online/offline)
# =============================================================
# =============================================================
# A/U DETEKCE — PODLE TVÉHO ALGORITMU (UP, UP1, UP2…)
# =============================================================
ANALYZE_DEBUG = True  # nastav na False, až tě debug přestane zajímat


# =============================================================
# A/U DETEKCE – PŘESNĚ DLE ZADÁNÍ (UP, UP1, UP2…)
# =============================================================
ANALYZE_DEBUG = True  # případně přepni na False

def analyze_pump_cycles(data):
    """
    Sjednocená A/U detekce podle zadání.
    - UP → cluster → nejlepší šestice → finalize
    - posun o start_idx + 3
    - střídání typů (A→U→A…)
    - forced-A při překročení hladiny pumpy
    - výstup: events i cycles
    """

    events = []
    cycles = []
    last_type = None

    n = len(data)
    if n < REG_WIN:
        return events, cycles

    # ---------------------------
    # Pomocné funkce
    # ---------------------------
    def linreg3(pts):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        xm = sum(xs)/3
        ym = sum(ys)/3
        num = sum((xs[k]-xm)*(ys[k]-ym) for k in range(3))
        den = sum((xs[k]-xm)**2 for k in range(3))
        m = 0 if den == 0 else num/den
        b = ym - m*xm
        return m, b

    def compute_window(i):
        if i + REG_WIN > n:
            return None

        win = data[i:i+REG_WIN]
        t0 = win[0][0]

        pts = [((t - t0).total_seconds(), float(l)) for (t, l) in win]
        left  = pts[:3]
        right = pts[3:6]

        mL, bL = linreg3(left)
        mR, bR = linreg3(right)

        if abs(mL - mR) < 1e-12:
            return None

        x_int = (bR - bL) / (mL - mR)
        y_int = mL * x_int + bL

        phiL = math.atan(mL)
        phiR = math.atan(mR)
        dphi = phiR - phiL

        if dphi > math.pi:
            dphi -= 2*math.pi
        elif dphi < -math.pi:
            dphi += 2*math.pi

        ang_signed = math.degrees(dphi)
        ang = abs(ang_signed)

        if ang < ANGLE_THRESHOLD_DEG:
            return None

        typ = "U" if ang_signed > 0 else "A"

        return {
            "idx": i,
            "typ": typ,
            "up": ang,
            "sign": ang_signed,
            "x": x_int,
            "y": y_int,
            "times": [t for (t,l) in win],
            "vals": [l for (t,l) in win]
        }

    def finalize(best):
        typ = best["typ"]
        times = best["times"]
        vals  = best["vals"]

        t3 = times[2]
        t4 = times[3]

        t_int = times[0] + timedelta(seconds = best["x"])

        if t3 <= t_int <= t4:
            return typ, t_int, best["y"], vals

        if typ == "U":
            k = vals.index(min(vals))
        else:
            k = vals.index(max(vals))

        return typ, times[k], vals[k], vals

    # ---------------------------
    # HLAVNÍ SMYČKA — forced A INCLUDED
    # ---------------------------
    i = 0
    while i <= n - REG_WIN:

        # ---------- 1) FORCED-A ----------
        # pokud aktuální hladina překročí PUMP_LIMIT_CM → okamžité A
        if data[i][1] >= PUMP_LIMIT_CM:
            if last_type != "A":  # jen pokud není poslední také A
                events.append(("A", data[i][0], data[i][1], ["FORCED"]))
                last_type = "A"
            i += 1
            continue

        # ---------- 2) REGRESNÍ OKNO ----------
        first = compute_window(i)
        if not first:
            i += 1
            continue

        base = first
        j = i + 1

        # cluster
        while j <= n - REG_WIN:
            nxt = compute_window(j)
            if not nxt:
                break
            if nxt["typ"] != base["typ"]:
                break
            if nxt["up"] > base["up"]:
                base = nxt
                j += 1
                continue
            break

        typ, t_ev, lvl, win = finalize(base)

        # střídání typů
        if last_type is None or typ != last_type:
            events.append((typ, t_ev, lvl, win))
            last_type = typ
        else:
            # stejné dva po sobě se ignorují
            pass

        # posun o +3
        i = base["idx"] + 3

    # ---------------------------
    # KONSTRUKCE CYKLŮ (A→U / U→A)
    # ---------------------------
    for k in range(1, len(events)):
        e1 = events[k-1]
        e2 = events[k]

        typ_prev, t1, h1, _ = e1
        _, t2, h2, _ = e2

        stav = "CERPANI" if typ_prev == "A" else "PLNENI"

        dt = max((t2 - t1).total_seconds(), 0)
        v1 = cm_to_liters(h1)
        v2 = cm_to_liters(h2)
        dv = abs(v2 - v1)
        r  = dv/dt if dt > 0 else 0

        if stav == "PLNENI":
            nat = r; vyc = ""; vyk = 0; pct = 0
        else:
            vyc = r; nat = ""; vyk = r; pct = 50

        cycles.append({
            "typ": typ_prev,
            "zacatek": t1,
            "konec": t2,
            "hlad_zac": h1,
            "hlad_kon": h2,
            "trvani_s": dt,
            "objem_l": dv,
            "rychlost_natoku": nat,
            "rychlost_vycerpani": vyc,
            "vykon_cerpadla": vyk,
            "pomer_pct": pct,
            "stav": stav,
            "window": "|".join(str(int(x)) for x in e1[3])
        })

    return events, cycles



# =============================================================
# FETCH DAT ZE SERVERU
# =============================================================
def fetch_once(url):
    if requests is None:
        dprint("[FETCH] requests library not available.")
        return None, None

    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        txt = r.text.strip()
        parts = txt.split("|")
        if len(parts) >= 3:
            try:
                level_m = float(parts[0])
                level_cm = round(level_m * 100.0, 1)
            except Exception:
                return None, None
            ts = parts[2]
            return level_cm, ts
    except Exception as e:
        dprint(f"[FETCH] Chyba: {e}")

    return None, None


# =============================================================
# OFFLINE SIMULACE — čtení bod po bodu + sjednocená analýza
# =============================================================
def simulate_from_file(csv_path):
    """
    Offline režim:
    - načte celý CSV soubor (čas, hladina),
    - data jsou brána postupně v pořadí v souboru,
    - po načtení všech bodů jednou spustí analyze_pump_cycles(data_points),
    - všechny cykly zapíše do OFFLINE souboru cykly_cerpadla_OFFLINE.csv
    """

    global OFFLINE_FLAG
    OFFLINE_FLAG = True

    if not os.path.exists(csv_path):
        print(f"[SIM] Soubor nenalezen: {csv_path}")
        return

    data_points = []
    total_points = 0

    try:
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=";")
            header = next(reader, None)

            for row in reader:
                if not row or len(row) < 3:
                    continue

                cas_raw = row[0].strip()
                t = None
                for fmt in ("%d.%m.%Y %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%d.%m.%Y %H:%M"):
                    try:
                        t = datetime.strptime(cas_raw, fmt)
                        break
                    except Exception:
                        continue
                if t is None:
                    try:
                        t = datetime.fromisoformat(cas_raw)
                    except Exception:
                        continue

                try:
                    level = float(row[2].strip().replace(",", "."))
                except Exception:
                    continue

                data_points.append((t, level))
                total_points += 1

        if not data_points:
            print("[SIM] Žádná platná data.")
            return

        print(f"[SIM] Načteno {total_points} bodů. Spouštím A/U analýzu...")

        cycles = analyze_pump_cycles(data_points)

        print(f"[SIM] Detekováno {len(cycles)} cyklů. Zapisuji do OFFLINE CSV...")

        for c in cycles:
            write_cycle_to_csv(c, offline=True)

        print("[SIM] Hotovo.")

    except Exception as e:
        print("[SIM] Chyba při čtení/offline analýze:", e)

# =============================================================
# ONLINE LOOP — sjednocená A/U analýza 8.4
# =============================================================
def run_online_loop(poll_interval=SAMPLE_INTERVAL):
    """
    Online režim:
    - každých poll_interval sekund načte bod ze stanice,
    - zapisuje měření do MERENI CSV,
    - udržuje buffer posledních BUFFER_LIMIT_ONLINE bodů,
    - na každém kroku spustí analyze_pump_cycles(data_points),
    - nové cykly zapisuje do CYKLY CSV,
    - současně běží SMS alarm 6/2 podle hladiny.
    """
    global OFFLINE_FLAG
    OFFLINE_FLAG = False

    print("[START] Online monitoring – A/U detekce 8.4 (sjednocená)")

    current_date = datetime.now().strftime("%d%m%y")
    mer_path = get_output_path(MERENI_FILENAME_BASE, by_date=current_date, offline=False, ext=".csv")
    cyk_path = get_output_path(CYKLY_FILENAME_BASE,  by_date=current_date, offline=False, ext=".csv")

    safe_append_csv(mer_path, header=[
        "cas_ze_stanice", "diff_cas_sec", "hladina_cm",
        "diff_hladina", "objem_l", "diff_objem_l"
    ])

    safe_append_csv(cyk_path, header=[
        "datum", "stav", "ZACATEK", "hlad_ZAC[cm]", "objem_ZAC[l]",
        "KONEC", "hlad_KON[cm]", "objem_KON[l]", "trvani[s]", "objem_l",
        "v_NATOK[l/s]", "v_ODCERP[l/s]", "CERPADLO[l/s]", "ZATIZENI[%]", "body"
    ])

    data_points = []
    last_time = None
    last_level = None
    last_volume = None
    last_cycle_end = None

    alarm_buffer = []
    alarm_active = False

    try:
        while True:
            # 1) načtení bodu
            level_cm, ts_str = fetch_once(SMARTWELL_URL)
            if level_cm is None:
                dprint("[FETCH] Data nenačtena, čekám...")
                time.sleep(poll_interval)
                continue

            # 2) parsování času
            tstamp = None
            for fmt in ("%d.%m.%Y %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%d.%m.%Y %H:%M"):
                try:
                    tstamp = datetime.strptime(ts_str, fmt)
                    break
                except Exception:
                    continue
            if tstamp is None:
                try:
                    tstamp = datetime.fromisoformat(ts_str)
                except Exception:
                    tstamp = datetime.now()

            # 3) přepnutí dne → nové soubory
            now_date = datetime.now().strftime("%d%m%y")
            if now_date != current_date:
                print(f"[DATE CHANGE] {current_date} -> {now_date}")
                current_date = now_date

                # necháme si jen pár posledních bodů pro návaznost
                data_points = data_points[-10:]
                last_cycle_end = None

                mer_path = get_output_path(MERENI_FILENAME_BASE, by_date=current_date, offline=False, ext=".csv")
                cyk_path = get_output_path(CYKLY_FILENAME_BASE,  by_date=current_date, offline=False, ext=".csv")

                safe_append_csv(mer_path, header=[
                    "cas_ze_stanice", "diff_cas_sec", "hladina_cm",
                    "diff_hladina", "objem_l", "diff_objem_l"
                ])

                safe_append_csv(cyk_path, header=[
                    "datum", "stav", "ZACATEK", "hlad_ZAC[cm]", "objem_ZAC[l]",
                    "KONEC", "hlad_KON[cm]", "objem_KON[l]", "trvani[s]", "objem_l",
                    "v_NATOK[l/s]", "v_ODCERP[l/s]", "CERPADLO[l/s]", "ZATIZENI[%]", "body"
                ])

            # 4) potlačení duplicit časů
            if last_time and abs((tstamp - last_time).total_seconds()) < 0.5:
                time.sleep(poll_interval)
                continue

            # 5) výpočet objemu a rozdílů
            volume_l = cm_to_liters(level_cm)

            if last_time is None:
                dt_sec = poll_interval
                diff_h = 0.0
                diff_v = 0.0
            else:
                dt_sec = (tstamp - last_time).total_seconds()
                if dt_sec <= 0:
                    dt_sec = poll_interval
                diff_h = round(level_cm - last_level, 1)
                diff_v = round(volume_l - last_volume, 1)

            row_mer = [
                tstamp.strftime("%d.%m.%Y %H:%M:%S"),
                f"{round(dt_sec, 1)}",
                f"{round(level_cm, 1)}",
                f"{round(diff_h, 1)}",
                f"{round(volume_l, 1)}",
                f"{round(diff_v, 1)}",
            ]
            safe_append_csv(mer_path, row_mer)

            last_time = tstamp
            last_level = level_cm
            last_volume = volume_l

            # 6) SMS alarm 6/2
            alarm_buffer.append(level_cm)
            if len(alarm_buffer) > 300:
                alarm_buffer = alarm_buffer[-300:]

            # aktivace
            if (
                len(alarm_buffer) >= SMS_ALARM_COUNT
                and all(x > SMS_ALARM_LIMIT for x in alarm_buffer[-SMS_ALARM_COUNT:])
            ):
                if not alarm_active:
                    if not OFFLINE_FLAG:
                        ok, st = send_sms_robust(
                            SMS_PHONE,
                            f"ALARM! Hladina nad limitem ({int(level_cm)} cm) v {datetime.now().strftime('%H:%M')}"
                        )
                        print("[ALARM] Aktivace:", ok, st)
                    else:
                        print("[ALARM] (OFFLINE) aktivace")
                    alarm_active = True

            # deaktivace
            if (
                alarm_active
                and len(alarm_buffer) >= SMS_CLEAR_COUNT
                and all(x < SMS_ALARM_LIMIT for x in alarm_buffer[-SMS_CLEAR_COUNT:])
            ):
                if not OFFLINE_FLAG:
                    ok, st = send_sms_robust(
                        SMS_PHONE,
                        f"KONEC ALARMU, hladina klesla ({int(level_cm)} cm) v {datetime.now().strftime('%H:%M')}"
                    )
                    print("[ALARM] Deaktivace:", ok, st)
                else:
                    print("[ALARM] (OFFLINE) deaktivace")
                alarm_active = False

            # 7) přidání do bufferu + omezení velikosti
            data_points.append((tstamp, level_cm))
            if len(data_points) > BUFFER_LIMIT_ONLINE:
                data_points = data_points[-BUFFER_LIMIT_ONLINE:]

            # 8) A/U detekce — sjednocená (stejná jako offline)
                # 8) A/U detekce — právě jako offline, ale jen NOVÉ cykly
                if len(data_points) >= REG_WIN + 2:
                    evs, cyc = analyze_pump_cycles(data_points)
                    for c in cyc:
                        kon = c.get("konec")
                        if last_cycle_end is None or kon > last_cycle_end:
                            write_cycle_to_csv(c, offline=False)
                            last_cycle_end = kon

                for c in cycles:
                    end_t = c.get("konec") or c.get("cas")
                    if last_cycle_end is None or (end_t and end_t > last_cycle_end):
                        write_cycle_to_csv(c, offline=False)
                        last_cycle_end = end_t

            if DEBUG:
                dprint(f"[MEAS] {tstamp.strftime('%H:%M:%S')} | {level_cm:.1f} cm")

            time.sleep(poll_interval)

    except KeyboardInterrupt:
        print("\n[STOP] Ukončeno uživatelem.")
    except Exception as e:
        print("[ERROR] Smyčka online selhala:", e)
    finally:
        print("[EXIT] Ukládám a končím.")


# =============================================================
# TEST SMS MODE
# =============================================================
def run_sms_test(force_port=None):
    print("=" * 60)
    print("📡 TESTOVACÍ REŽIM SMS – ODESÍLÁM 2x ZKUŠEBNÍ SMS")
    print("=" * 60)
    ok1, s1 = send_sms_robust(SMS_PHONE, "TEST SMS 01 Aktivace alarmu", force_port)
    print("Výsledek aktivace:", ok1, s1)
    time.sleep(1)
    ok2, s2 = send_sms_robust(SMS_PHONE, "TEST SMS 02 Deaktivace alarmu", force_port)
    print("Výsledek deaktivace:", ok2, s2)
    print("Test SMS dokončen.")


# =============================================================
# MAIN MENU
# =============================================================
if __name__ == "__main__":
    print("Verze 9.0 (z 8.0) - generováno", datetime.now().strftime("%d.%m.%Y %H:%M:%S"), "CET")
    print("=" * 80)
    print("Zvol režim běhu:")
    print("1 = Offline simulace z CSV")
    print("2 = Online měření (s aktivními SMS)")
    print("3 = Testovací SMS režim")
    choice = input("👉 Volba (1/2/3): ").strip()

    if choice == "1":
        csv_file = os.path.join(BASE_DIR, "mereni_jimka.csv")
        if not os.path.exists(csv_file):
            csv_file = os.path.join(os.getcwd(), "mereni_jimka.csv")
        simulate_from_file(csv_file)

    elif choice == "2":
        run_online_loop()

    elif choice == "3":
        run_sms_test()

    else:
        print("Neplatná volba — končím.")
