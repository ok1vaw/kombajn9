
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

# =============================================================
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
SMS_ALARM_LIMIT    = 165                 # cm
SMS_ALARM_COUNT    = 6                   # aktivace: 6 posledních bodů nad limitem
SMS_CLEAR_COUNT    = 2                   # deaktivace: 2 poslední body pod limitem

# A/U detection parametry
REG_WIN               = 6                # délka regresního okna
ANGLE_THRESHOLD_DEG   = 25.0             # minimální |úhel| pro kandidáta A/U
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
def analyze_pump_cycles(data):
    """
    Verze 8.4 — sjednocená detekce A/U (stejná pro online i offline):

    - data: seznam [(datetime, hladina_cm), ...] v časovém pořadí
    - použije 6-bodová regresní okna (REG_WIN=6),
    - hledá shluky oken, kde |úhel| > ANGLE_THRESHOLD_DEG,
    - ze shluku vezme to okno, kde |úhel| je maximální,
    - z tohoto okna určí čas A/U:
        * pokud průsečík leží mezi 3. a 4. bodem → použije průsečík,
        * jinak vezme max/min z celé šestice dle znaménka úhlu (A = vrchol, U = dno),
    - z posloupnosti A/U bodů vytvoří cykly (A→U, U→A) stejně jako dříve.

    Návratová hodnota:
        list(cycles) — každý cyklus je dict kompatibilní s write_cycle_to_csv().
    """

    n = len(data)
    if n < REG_WIN + 2:
        return []

    # ---------------------------------------------------------
    # 1) Jedno 6-bodové okno → kandidát
    # ---------------------------------------------------------
    def compute_window_event(start_idx):
        if start_idx + REG_WIN > len(data):
            return None

        win = data[start_idx:start_idx + REG_WIN]
        t0 = win[0][0]

        pts = [((t - t0).total_seconds(), float(h)) for t, h in win]
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]

        left = pts[0:3]
        right = pts[3:6]

        def linreg(points):
            m = 0.0
            b = 0.0
            n_local = len(points)
            if n_local == 0:
                return 0.0, 0.0
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            xm = sum(xs) / n_local
            ym = sum(ys) / n_local
            num = sum((xs[i] - xm) * (ys[i] - ym) for i in range(n_local))
            den = sum((xs[i] - xm) ** 2 for i in range(n_local))
            if den != 0:
                m = num / den
            b = ym - m * xm
            return m, b

        mL, bL = linreg(left)
        mR, bR = linreg(right)

        # průsečík
        if abs(mL - mR) < 1e-12:
            return None

        x_int = (bR - bL) / (mL - mR)
        y_int = mL * x_int + bL

        # ořez y_int do rozsahu okna
        min_h = min(ys)
        max_h = max(ys)
        y_int = max(min(y_int, max_h), min_h)

        # orientovaný úhel
        phiL = math.atan(mL)
        phiR = math.atan(mR)
        delta = phiR - phiL
        if delta > math.pi:
            delta -= 2 * math.pi
        elif delta < -math.pi:
            delta += 2 * math.pi

        angle_signed_deg = math.degrees(delta)
        angle_abs_deg = abs(angle_signed_deg)

        if angle_abs_deg < ANGLE_THRESHOLD_DEG:
            return None

        # typ podle znaménka
        if angle_signed_deg > 0:
            typ = "U"
        elif angle_signed_deg < 0:
            typ = "A"
        else:
            return None

        event_time = t0 + timedelta(seconds=x_int)

        return {
            "typ": typ,
            "cas": event_time,
            "hladina_int": float(y_int),
            "angle_deg": angle_abs_deg,
            "angle_signed": angle_signed_deg,
            "start_idx": start_idx,
            "levels": ys,
            "times": [t for t, _ in win],
        }

    # ---------------------------------------------------------
    # 2) Finalizace eventu z nejlepšího okna ve shluku
    # ---------------------------------------------------------
    def finalize_event(ev):
        start_idx = ev["start_idx"]
        win = data[start_idx:start_idx + REG_WIN]
        times = [t for t, _ in win]
        levels = [float(h) for _, h in win]

        t3 = times[2]
        t4 = times[3]
        ev_time = ev["cas"]

        # průsečík mezi 3. a 4. bodem?
        if t3 <= ev_time <= t4:
            time_final = ev_time
            level_final = ev["hladina_int"]
        else:
            # fallback: lokální extrém v šestici
            if ev["typ"] == "U":
                idx_ext = min(range(REG_WIN), key=lambda i: levels[i])
            else:  # "A"
                idx_ext = max(range(REG_WIN), key=lambda i: levels[i])
            time_final = times[idx_ext]
            level_final = levels[idx_ext]

        return {
            "typ": ev["typ"],
            "cas": time_final,
            "hladina": round(level_final, 1),
            "window": "|".join(str(int(round(h))) for h in levels),
            "angle": round(ev["angle_deg"], 2),
            "angle_signed": round(ev["angle_signed"], 2),
        }

    # ---------------------------------------------------------
    # 3) Skenování dat: shluky kandidátů, z každého max |úhel|
    # ---------------------------------------------------------
    events = []
    i = 0
    while i <= n - REG_WIN:
        ev0 = compute_window_event(i)

        if ev0 is None:
            i += 1
            continue

        # máme první kandidát ve shluku
        best_ev = ev0
        j = i + 1

        while j <= n - REG_WIN:
            evj = compute_window_event(j)
            if evj is None:
                break  # konec shluku

            if evj["typ"] != best_ev["typ"]:
                # změna typu → ukončíme shluk, ať nemícháme A/U dohromady
                break

            if evj["angle_deg"] > best_ev["angle_deg"]:
                best_ev = evj

            j += 1

        # z nejlepší šestice vytvoříme finální event
        best_final = finalize_event(best_ev)

        # sloučení duplicit: nechceme A,A,A nebo U,U,U
        if events and events[-1]["typ"] == best_final["typ"]:
            # ponecháme ten s větším úhlem
            if best_final["angle"] > events[-1]["angle"]:
                events[-1] = best_final
        else:
            events.append(best_final)

        # posuneme se za shluk
        i = j

    # ---------------------------------------------------------
    # 4) Z A/U eventů sestavíme cykly
    # ---------------------------------------------------------
    if len(events) < 2:
        return []

    cycles = []

    # pomocná proxy na cm_to_liters, aby funkce byla použitelná i samostatně
    def cm_to_liters_local(h):
        try:
            return cm_to_liters(h)
        except Exception:
            return float(h)

    for k in range(1, len(events)):
        prev = events[k - 1]
        curr = events[k]

        t_start = prev["cas"]
        t_end   = curr["cas"]
        trvani  = max((t_end - t_start).total_seconds(), 0.0)

        h1 = float(prev["hladina"])
        h2 = float(curr["hladina"])
        v1 = cm_to_liters_local(h1)
        v2 = cm_to_liters_local(h2)
        diff_v = abs(v2 - v1)

        # stav podle předchozího eventu: A → ČERPÁNÍ, U → PLNĚNÍ
        stav = "CERPANI" if prev["typ"] == "A" else "PLNENI"
        rychlost = diff_v / trvani if trvani > 0 else 0.0

        if stav == "PLNENI":
            rych_nat = rychlost
            rych_vyc = ""
            vykon = 0.0
            pomer = 0.0
        else:  # CERPANI
            rych_vyc = rychlost
            rych_nat = prev.get("rychlost_natoku", 0.0)
            if rych_nat in ("", None):
                rych_nat = 0.0
            vykon = rych_vyc + rych_nat
            t_in = prev.get("trvani_s", trvani)
            t_out = trvani
            pomer = (t_out / (t_out + t_in) * 100.0) if (t_out + t_in) > 0 else 0.0

        cycle = {
            "typ": prev["typ"],
            "cas": prev["cas"],
            "zacatek": t_start,
            "konec": t_end,
            "hlad_zac": round(h1, 1),
            "hlad_kon": round(h2, 1),
            "trvani_s": trvani,
            "objem_l": round(diff_v, 1),
            "rychlost_natoku": round(rych_nat, 1) if stav == "PLNENI" else "",
            "rychlost_vycerpani": round(rych_vyc, 1) if stav == "CERPANI" else "",
            "vykon_cerpadla": round(vykon, 1),
            "pomer_pct": round(pomer, 1),
            "stav": stav,
            "window": prev["window"],
        }

        cycles.append(cycle)

    return cycles
# =============================================================
# ZÁPIS CYKLU DO CSV
# =============================================================
def write_cycle_to_csv(cycle, offline=False):
    """
    Zapíše jeden cyklus do CSV. Pokud offline=True, použije suffix _OFFLINE.
    """
    zacatek = cycle.get("zacatek") or cycle.get("cas") or datetime.now()
    out_path = get_output_path(CYKLY_FILENAME_BASE, by_date=zacatek, offline=offline, ext=".csv")

    header = [
        "datum", "stav", "ZACATEK", "hlad_ZAC[cm]", "objem_ZAC[l]",
        "KONEC", "hlad_KON[cm]", "objem_KON[l]", "trvani[s]", "objem_l",
        "v_NATOK[l/s]", "v_ODCERP[l/s]", "CERPADLO[l/s]", "ZATIZENI[%]", "body"
    ]

    datum = zacatek.strftime("%d.%m.%Y")
    stav = cycle.get("stav", "")

    zac_str = cycle.get("zacatek").strftime("%H:%M:%S") if cycle.get("zacatek") else ""
    hlad_zac = cycle.get("hlad_zac", "")
    if isinstance(hlad_zac, (int, float)):
        hlad_zac = round(hlad_zac, 1)
    else:
        hlad_zac = ""
    obj_zac = round(cm_to_liters(hlad_zac), 1) if hlad_zac != "" else ""

    konec_str = cycle.get("konec").strftime("%H:%M:%S") if cycle.get("konec") else ""
    hlad_kon = cycle.get("hlad_kon", "")
    if isinstance(hlad_kon, (int, float)):
        hlad_kon = round(hlad_kon, 1)
    else:
        hlad_kon = ""
    obj_kon = round(cm_to_liters(hlad_kon), 1) if hlad_kon != "" else ""

    trvani = round(cycle.get("trvani_s", 0), 1)
    obj_l  = round(cycle.get("objem_l", 0), 1)

    v_nat = cycle.get("rychlost_natoku", "")
    v_vyc = cycle.get("rychlost_vycerpani", "")
    vykon = cycle.get("vykon_cerpadla", "")
    pomer = cycle.get("pomer_pct", "")
    if isinstance(pomer, (int, float)):
        pomer = round(pomer, 1)
    else:
        pomer = ""

    body = cycle.get("window", "")

    # prázdné cykly ignorujeme
    if trvani == 0 and obj_l == 0:
        dprint("[WRITE] Prázdný cyklus (trvani=0 a obj_l=0) – nepíšu.")
        return False

    row = [
        datum, stav, zac_str, hlad_zac, obj_zac,
        konec_str, hlad_kon, obj_kon,
        trvani, obj_l, v_nat, v_vyc, vykon, pomer, body
    ]

    ok = safe_append_csv(out_path, row, header=header)
    if DEBUG:
        dprint(f"[WRITE_CYCLE] -> {os.path.basename(out_path)} | OK={ok} | {row}")
    return ok


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
    - čte vstupní CSV bod po bodu,
    - každý nový bod přidá do data_points,
    - na každém kroku spustí analyze_pump_cycles(data_points),
    - nové cykly (podle času konce) zapisuje do OFFLINE CSV.
    """
    global OFFLINE_FLAG
    OFFLINE_FLAG = True

    if not os.path.exists(csv_path):
        print(f"[SIM] Soubor nenalezen: {csv_path}")
        return

    print(f"[SIM] Offline simulace: {csv_path}")

    data_points = []
    last_time = None
    last_cycle_end = None
    total_points = 0
    total_cycles = 0

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

                # potlačení duplicitních časů
                if last_time and abs((t - last_time).total_seconds()) < 0.5:
                    continue

                data_points.append((t, level))
                total_points += 1
                last_time = t

                # zavoláme sjednocenou detekci
                cycles = analyze_pump_cycles(data_points)

                # zapíšeme jen nové cykly (podle koncového času)
                for c in cycles:
                    end_t = c.get("konec") or c.get("cas")
                    if last_cycle_end is None or (end_t and end_t > last_cycle_end):
                        write_cycle_to_csv(c, offline=True)
                        last_cycle_end = end_t
                        total_cycles += 1

        print(f"[SIM] Načteno {total_points} bodů, detekováno {total_cycles} cyklů (offline).")

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
            if len(data_points) >= REG_WIN + 2:
                cycles = analyze_pump_cycles(data_points)

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
