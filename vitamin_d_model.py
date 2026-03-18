"""
vitamin_d_model.py
==================
A computational implementation of the plasma 25(OH)D (vitamin D) model
described in the accompanying LaTeX document, based on the pharmacokinetic
framework of Diffey (2013), with demographic corrections from:
  - Hilger et al. (2014) — sex independence
  - Young (2020)         — Fitzpatrick skin-type scaling
  - Chalcraft et al. (2020) — age-related synthesis decline

Two simulation entry points are provided:

  simulate_subject()
      Legacy interface.  Accepts raw UVB (J/m²) + BSA% per day.
      Use when wearable sensor CSV data is NOT available.

  simulate_subject_from_logs()
      Sensor-data interface.  Accepts per-body-part SED logs produced by
      data_loader.load_subject_logs() via bsa.effective_doses_from_daily_logs().
      The dose-area product E(t)*A(t) is computed directly from the sensor
      readings instead of being approximated from a clothing survey.
      Use this path for all subjects with wearable devices.

All equation numbers in the comments refer to the source document.

Author  : Claude (Anthropic, claude-sonnet-4-6)
Purpose : Reference implementation for research / educational use
"""

import math


# =============================================================================
# FIXED BIOLOGICAL PARAMETERS
# (Table 1 in source document — from Diffey 2013)
# =============================================================================

f      = 0.15   # Fraction of vitamin D stored in tissue (dimensionless)
beta   = 25     # Plasma clearance half-life (days)
gamma  = 250    # Tissue-store clearance half-life (days)
alpha  = 0.6    # Half-time for UV-derived vitamin D uptake into plasma (days)
alpha_ = 1.5    # Half-time for oral vitamin D uptake into plasma (days)  [alpha']
A_uv   = 0.18   # UV scaling factor (nmol/L per SED per % body surface area)
S      = 0.023  # Oral scaling factor (nmol/L per ug)


# =============================================================================
# DEMOGRAPHIC SCALAR FUNCTIONS
# =============================================================================

def skin_factor(fitzpatrick_type: int) -> float:
    """
    Eq. 9  --  Fitzpatrick skin-type UV exposure correction factor.

    Darker skin contains more melanin, which competes with 7-dehydrocholesterol
    for UVB photons, reducing vitamin D synthesis per unit of UV dose.
    Young (2020) found that Types III-VI require ~1.35x more UV exposure to
    achieve equivalent synthesis compared to Types I-II.

    Parameters
    ----------
    fitzpatrick_type : int
        Skin type on the Fitzpatrick scale (1-6).

    Returns
    -------
    float
        f_skin -- dimensionless scalar in range (0, 1].
    """
    if fitzpatrick_type in (1, 2):
        return 1.00
    elif fitzpatrick_type in (3, 4, 5, 6):
        return 1.0 / 1.35   # approx 0.7407
    else:
        raise ValueError(f"Fitzpatrick type must be 1-6, got {fitzpatrick_type}")


def age_factor(age: float) -> float:
    """
    Eq. 10  --  Age-related decline in cutaneous vitamin D synthesis.

    Chalcraft et al. (2020) measured a ~1.3% per year decline in synthesis
    capacity, anchored at age 20 (f_age = 1.0).

    Parameters
    ----------
    age : float
        Subject's age in years.

    Returns
    -------
    float
        f_age -- dimensionless scalar (nominally in (0, 1] for ages 20-97).
    """
    return 1.0 - 0.013 * (age - 20)


# =============================================================================
# UV DOSE CONVERSION
# =============================================================================

def uvb_dose_to_sed(uvb_j_per_m2: float) -> float:
    """
    Eq. 1 (adapted)  --  Convert daily UVB dose to Standard Erythemal Doses.

    Conversion:  1 SED = 100 J/m²  (ISO 17166)
    Therefore:   SED = UVB_dose (J/m²) / 100

    UVA is excluded — it does not drive vitamin D photosynthesis.

    Parameters
    ----------
    uvb_j_per_m2 : float
        Total daily UVB irradiance dose in J/m².

    Returns
    -------
    float
        E(t) in SED.
    """
    return uvb_j_per_m2 / 100.0


# =============================================================================
# RESPONSE FUNCTIONS
# =============================================================================

def R_UV(t: float) -> float:
    """
    Eq. 2  --  UV-derived vitamin D impulse response function.

    Describes the rise-and-fall kinetics of plasma 25(OH)D following a single
    unit UV exposure (1 SED over 1% body area).

    Parameters
    ----------
    t : float
        Time since exposure (days, >= 0).

    Returns
    -------
    float
        Response in nmol/L per (SED × % body area).
    """
    return A_uv * (
        (1 - f) * 2 ** (-t / beta)
        + f     * 2 ** (-t / gamma)
        - 2     ** (-t / alpha)
    )


def R_UV_eff(t: float) -> float:
    """
    Sensor-mode impulse response for the effective dose-area formulation.

    Identical mathematics to R_UV but the input unit is
    (SED × fractional BSA) rather than (SED × % BSA).
    The scaling factor A_uv is multiplied by 100 to compensate for the
    fractional-vs-percent difference:

        R_UV_eff = 100 × R_UV   because   E*A(%) = 100 × E*A(fraction)

    Parameters
    ----------
    t : float
        Time since exposure (days, >= 0).

    Returns
    -------
    float
        Response in nmol/L per (SED × fractional BSA).
    """
    return 100.0 * R_UV(t)


def R_oral(t: float) -> float:
    """
    Eq. 3  --  Oral vitamin D impulse response function.

    Parameters
    ----------
    t : float
        Time since ingestion (days, >= 0).

    Returns
    -------
    float
        Response in nmol/L per ug.
    """
    return S * (
        (1 - f) * 2 ** (-t / beta)
        + f     * 2 ** (-t / gamma)
        - 2     ** (-t / alpha_)
    )


# =============================================================================
# CONVOLUTION ACCUMULATORS
# =============================================================================

def compute_C_oral(oral_doses: list) -> list:
    """
    Eq. 4  --  Cumulative oral contribution to plasma 25(OH)D.

    Parameters
    ----------
    oral_doses : list
        O(t) for t = 0 … N-1  (ug/day).

    Returns
    -------
    list
        C_oral(T) for T = 0 … N-1  (nmol/L).
    """
    N = len(oral_doses)
    C_oral = [0.0] * N
    for T in range(N):
        total = 0.0
        for t in range(T + 1):
            lag = T - t + 1
            total += oral_doses[t] * R_oral(lag)
        C_oral[T] = total
    return C_oral


def compute_C_sun(
    uv_doses:    list,
    body_areas:  list,
    age:         float,
    fitzpatrick: int,
) -> list:
    """
    Eq. 11  --  Cumulative UV contribution (survey / legacy path).

    Uses SED + % BSA inputs as in the original model.

        C_sun(T) = Σ_{t=0}^{T}  f_age × f_skin × E(t) × A(t) × R_UV(T-t+1)

    Parameters
    ----------
    uv_doses   : list  -- E(t) in SED per day.
    body_areas : list  -- A(t) in % BSA per day.
    age        : float        -- Subject age (years).
    fitzpatrick: int          -- Fitzpatrick skin type 1-6.

    Returns
    -------
    list
        C_sun(T) for T = 0 … N-1  (nmol/L).
    """
    f_age_val  = age_factor(age)
    f_skin_val = skin_factor(fitzpatrick)
    demo       = f_age_val * f_skin_val

    N = len(uv_doses)
    C_sun = [0.0] * N
    for T in range(N):
        total = 0.0
        for t in range(T + 1):
            lag    = T - t + 1
            total += demo * uv_doses[t] * body_areas[t] * R_UV(lag)
        C_sun[T] = total
    return C_sun


def compute_C_sun_effective(
    uv_eff_doses: list,
    age:          float,
    fitzpatrick:  int,
) -> list:
    """
    Sensor-mode cumulative UV contribution  (replaces compute_C_sun when
    per-body-part SED measurements are available from data_loader.py).

    The effective dose-area product  D_eff(t) = Σ_i SED_i(t) × BSA_fraction_i
    is read directly from the wearable sensor log via bsa.py and passed in as
    ``uv_eff_doses``.  This collapses the two-argument (E, A%) form into a
    single pre-computed scalar per day:

        C_sun(T) = Σ_{t=0}^{T}  f_age × f_skin × D_eff(t) × R_UV_eff(T-t+1)

    No BSA list is required — it is already embedded in D_eff.

    Parameters
    ----------
    uv_eff_doses : list
        D_eff(t) — effective dose-area product per day
        (SED × fractional BSA, from bsa.effective_doses_from_daily_logs()).
    age          : float  -- Subject age (years).
    fitzpatrick  : int    -- Fitzpatrick skin type 1-6.

    Returns
    -------
    list
        C_sun(T) for T = 0 … N-1  (nmol/L).
    """
    f_age_val  = age_factor(age)
    f_skin_val = skin_factor(fitzpatrick)
    demo       = f_age_val * f_skin_val

    N = len(uv_eff_doses)
    C_sun = [0.0] * N
    for T in range(N):
        total = 0.0
        for t in range(T + 1):
            lag    = T - t + 1
            total += demo * uv_eff_doses[t] * R_UV_eff(lag)
        C_sun[T] = total
    return C_sun


# =============================================================================
# SATURATION FACTOR
# =============================================================================

def saturation_factor(C_total_prev: float) -> float:
    """
    Eq. 5  --  Diminishing-returns saturation factor F(T).

        F(T) = exp( -0.01 × C_total(T-1) )

    Parameters
    ----------
    C_total_prev : float
        Plasma 25(OH)D on the previous day (nmol/L).

    Returns
    -------
    float
        F(T) in (0, 1].
    """
    return math.exp(-0.01 * C_total_prev)


# =============================================================================
# MASTER MODEL — total plasma 25(OH)D
# =============================================================================

def run_model(
    oral_doses:  list,
    uv_doses:    list,
    body_areas:  list,
    age:         float,
    fitzpatrick: int,
    C0:          float = 50.0,
) -> list:
    """
    Eq. 7  --  Total plasma 25(OH)D trajectory  (legacy / survey path).

    Uses raw SED + % BSA inputs.  For sensor-equipped subjects prefer
    run_model_from_effective_doses().

    Parameters
    ----------
    oral_doses   : list  -- O(t) in ug/day.
    uv_doses     : list  -- E(t) in SED/day.
    body_areas   : list  -- A(t) in % BSA/day.
    age          : float        -- Subject age (years).
    fitzpatrick  : int          -- Fitzpatrick skin type 1-6.
    C0           : float        -- Initial plasma 25(OH)D (nmol/L).

    Returns
    -------
    list
        C_total(T) for T = 0 … N-1  (nmol/L).
    """
    N = len(oral_doses)
    assert len(uv_doses)   == N, "All input lists must be the same length."
    assert len(body_areas) == N, "All input lists must be the same length."

    C_oral = compute_C_oral(oral_doses)
    C_sun  = compute_C_sun(uv_doses, body_areas, age, fitzpatrick)

    C_total = [0.0] * N
    for T in range(N):
        C_prev      = C0      if T == 0 else C_total[T - 1]
        C_oral_prev = 0.0     if T == 0 else C_oral[T - 1]
        C_sun_prev  = 0.0     if T == 0 else C_sun[T - 1]

        delta_oral = C_oral[T] - C_oral_prev
        delta_sun  = C_sun[T]  - C_sun_prev
        F          = saturation_factor(C_prev)

        C_total[T] = C_prev + delta_oral + F * delta_sun

    return C_total


def run_model_from_effective_doses(
    oral_doses:   list,
    uv_eff_doses: list,
    age:          float,
    fitzpatrick:  int,
    C0:           float = 50.0,
) -> list:
    """
    Eq. 7  --  Total plasma 25(OH)D trajectory  (sensor / CSV log path).

    Uses the pre-computed dose-area product from wearable sensor logs
    instead of the two-argument E(t) × A(t) form.  All other model
    mechanics (oral convolution, saturation, demographics) are unchanged.

    Parameters
    ----------
    oral_doses   : list  -- O(t) in ug/day.
    uv_eff_doses : list  -- D_eff(t) — effective dose-area product per day
                            (SED × fractional BSA), from
                            bsa.effective_doses_from_daily_logs() or
                            data_loader.extract_effective_doses().
    age          : float        -- Subject age (years).
    fitzpatrick  : int          -- Fitzpatrick skin type 1-6.
    C0           : float        -- Initial plasma 25(OH)D (nmol/L).

    Returns
    -------
    list
        C_total(T) for T = 0 … N-1  (nmol/L).
    """
    N = len(oral_doses)
    assert len(uv_eff_doses) == N, "oral_doses and uv_eff_doses must be the same length."

    C_oral = compute_C_oral(oral_doses)
    C_sun  = compute_C_sun_effective(uv_eff_doses, age, fitzpatrick)

    C_total = [0.0] * N
    for T in range(N):
        C_prev      = C0      if T == 0 else C_total[T - 1]
        C_oral_prev = 0.0     if T == 0 else C_oral[T - 1]
        C_sun_prev  = 0.0     if T == 0 else C_sun[T - 1]

        delta_oral = C_oral[T] - C_oral_prev
        delta_sun  = C_sun[T]  - C_sun_prev
        F          = saturation_factor(C_prev)

        C_total[T] = C_prev + delta_oral + F * delta_sun

    return C_total


# =============================================================================
# SIMULATION RUNNERS
# =============================================================================

def _print_header(name, age, fitzpatrick, C0, file=None):
    print("=" * 65, file=file)
    print(f"  {name}", file=file)
    print(f"    Age              : {age} years", file=file)
    print(f"    Fitzpatrick type : {fitzpatrick}", file=file)
    print(f"    f_age            : {age_factor(age):.3f}", file=file)
    print(f"    f_skin           : {skin_factor(fitzpatrick):.3f}", file=file)
    print(f"    Initial 25(OH)D  : {C0:.1f} nmol/L  (WARNING: assumed)", file=file)
    print("=" * 65, file=file)


def simulate_subject(
    name:          str,
    oral:          list,
    uvb_raw_j_m2:  list,
    body_area_pct: list,
    age:           float,
    fitzpatrick:   int,
    C0:            float,
    file=None,
) -> list:
    """
    Legacy interface: simulate using raw UVB (J/m²) + BSA% per day.

    Use this when no wearable sensor CSV data is available.
    For subjects with UV logger devices use simulate_subject_from_logs().

    Parameters
    ----------
    name           : Label for the printed output.
    oral           : Daily oral vitamin D intake (ug/day; 1 ug = 40 IU).
    uvb_raw_j_m2   : Daily UVB dose (J/m²).  UVA excluded.
    body_area_pct  : Body surface area exposed per day (%).
    age            : Subject age (years).
    fitzpatrick    : Fitzpatrick skin type 1-6.
    C0             : Initial plasma 25(OH)D (nmol/L).

    Returns
    -------
    list
        C_total(T) for each day  (nmol/L).
    """
    N = len(oral)
    assert len(uvb_raw_j_m2)  == N, "oral and uvb_raw_j_m2 must match in length."
    assert len(body_area_pct) == N, "oral and body_area_pct must match in length."

    uv_doses = [uvb_dose_to_sed(v) for v in uvb_raw_j_m2]
    result   = run_model(
        oral_doses  = oral,
        uv_doses    = uv_doses,
        body_areas  = body_area_pct,
        age         = age,
        fitzpatrick = fitzpatrick,
        C0          = C0,
    )

    _print_header(name, age, fitzpatrick, C0, file=file)
    print(f"  {'Day':>4}  {'Oral(ug)':>9}  {'UVB(J/m²)':>10}  {'BSA(%)':>7}  {'C_total(nmol/L)':>16}", file=file)
    print("  " + "-" * 57, file=file)
    for day in range(N):
        print(
            f"  {day:>4}  {oral[day]:>9.1f}  {uvb_raw_j_m2[day]:>10.1f}"
            f"  {body_area_pct[day]:>7.1f}  {result[day]:>16.2f}",
            file=file,
        )
    print(file=file)
    return result


def simulate_subject_from_logs(
    name:         str,
    oral:         list,
    daily_logs:   list,
    age:          float,
    fitzpatrick:  int,
    C0:           float,
    file=None,
) -> list:
    """
    Sensor-data interface: simulate using per-body-part SED CSV logs.

    ``daily_logs`` is the direct output of data_loader.load_subject_logs().
    The effective dose-area product is derived from the per-body-part SED
    measurements via bsa.effective_doses_from_daily_logs(), so no BSA survey
    is needed.

    Parameters
    ----------
    name        : Label for the printed output.
    oral        : Daily oral vitamin D intake (ug/day).
    daily_logs  : List of per-day summary dicts from data_loader.load_subject_logs().
    age         : Subject age (years).
    fitzpatrick : Fitzpatrick skin type 1-6.
    C0          : Initial plasma 25(OH)D (nmol/L).

    Returns
    -------
    list
        C_total(T) for each day  (nmol/L).
    """
    # Import here to avoid circular dependency at module level
    from bsa import effective_doses_from_daily_logs

    uv_eff_doses = effective_doses_from_daily_logs(daily_logs)
    N            = len(oral)

    assert len(uv_eff_doses) == N, (
        f"[{name}] oral ({len(oral)}) and daily_logs ({len(daily_logs)}) "
        "must cover the same number of days."
    )

    result = run_model_from_effective_doses(
        oral_doses   = oral,
        uv_eff_doses = uv_eff_doses,
        age          = age,
        fitzpatrick  = fitzpatrick,
        C0           = C0,
    )

    _print_header(name, age, fitzpatrick, C0, file=file)
    print(f"  {'Day':>4}  {'Oral(ug)':>9}  {'D_eff':>8}  {'%Outdoor':>9}  {'C_total(nmol/L)':>16}", file=file)
    print("  " + "-" * 57, file=file)
    for day in range(N):
        pct_out = daily_logs[day].get("fraction_outdoor", float("nan")) * 100
        print(
            f"  {day:>4}  {oral[day]:>9.1f}  {uv_eff_doses[day]:>8.4f}"
            f"  {pct_out:>8.1f}%  {result[day]:>16.2f}",
            file=file,
        )
    print(file=file)
    return result


# =============================================================================
# ORAL DOSE LOADER
# =============================================================================

def load_oral_doses_from_json(json_path, subject_key: str, dates: list) -> list:
    """
    Load daily oral vitamin D intakes from vitaminD-mapped.json.

    The JSON stores ``totalVitaminD_estimated`` in IU per day.
    This function converts to µg (÷ 40) to match the model's expected units.

    Parameters
    ----------
    json_path   : str or Path — path to vitaminD-mapped.json.
    subject_key : str         — top-level key in the JSON.
    dates       : list[str]   — calendar dates ("YYYY-MM-DD") derived from
                                the UV log file timestamps, one per day.

    Returns
    -------
    list[float]
        Oral dose in µg/day for each date, length len(dates).
    """
    import json

    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    iu_by_date = {rec["date"]: rec.get("totalVitaminD_estimated", 0.0) or 0.0
                  for rec in data.get(subject_key, [])}

    return [iu_by_date.get(d, 0.0) / 40.0 for d in dates]


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import sys
    import contextlib
    from pathlib import Path
    from datetime import datetime

    # -------------------------------------------------------------------------
    # Output file — written next to this script in formula/
    # -------------------------------------------------------------------------
    SCRIPT_DIR  = Path(__file__).parent
    OUTPUT_FILE = SCRIPT_DIR / "results.txt"

    # -------------------------------------------------------------------------
    # vitaminD-mapped.json — oral intake source (IU → µg conversion inside
    # load_oral_doses_from_json).  Expected to live next to this script.
    # -------------------------------------------------------------------------
    ORAL_JSON    = SCRIPT_DIR / "vitaminD-mapped.json"
    SCENARIO_UV_JSON = SCRIPT_DIR / "surface.json"

    # ---- Switch between real sensor data and scenario_exposure.json ----
    # Set USE_SCENARIO_UV = True  to use scenario_exposure.json for UV doses
    # Set USE_SCENARIO_UV = False to use the real resampled CSV sensor logs
    USE_SCENARIO_UV = True

    def run_all(out):
        """Run all simulations, writing every print() to `out`."""

        # =====================================================================
        # PATH A — Sensor data available
        # =====================================================================
        try:
            import json as _json
            from data_loader import load_subject_logs, print_log_summary

            DATA_ROOT = Path(__file__).parent / "data_resampled"

            sensor_subjects = [
                ("Oscar  -- Subject A -- ~21 y/o, Fitzpatrick II",
                 DATA_ROOT / "uv oscar",        "Oscar",   23, 2, 50.0),
                ("Nicole -- Subject B -- ~20 y/o, Fitzpatrick II",
                 DATA_ROOT / "Uv tests Nicole", "Nicole",  20, 2, 50.0),
                ("Subject C -- vi_logs",
                 DATA_ROOT / "vi_logs",         "vi_pdf",  22, 3, 50.0),
            ]

            # Load scenario UV data once if needed
            scenario_uv_data = {}
            if USE_SCENARIO_UV and SCENARIO_UV_JSON.exists():
                with open(SCENARIO_UV_JSON) as fh:
                    scenario_uv_data = _json.load(fh)
                print(f"[INFO] Using scenario UV data from {SCENARIO_UV_JSON.name}\n", file=out)

            for name, log_dir, subject_key, age, fitz, C0 in sensor_subjects:
                if not log_dir.is_dir():
                    print(f"[SKIP] Log directory not found: {log_dir}", file=out)
                    continue

                logs  = load_subject_logs(log_dir)
                dates = [d["date"] for d in logs]
                oral  = load_oral_doses_from_json(ORAL_JSON, subject_key, dates)

                print(
                    f"[INFO] {name}: loaded oral doses from JSON "
                    f"({len(oral)} days, key='{subject_key}')",
                    file=out,
                )

                # Patch sed_by_part from scenario_exposure.json if enabled
                if USE_SCENARIO_UV and subject_key in scenario_uv_data:
                    records_by_date = {r["date"]: r["sed_by_part"]
                                       for r in scenario_uv_data[subject_key]}
                    for log in logs:
                        if log["date"] in records_by_date:
                            log["sed_by_part"] = records_by_date[log["date"]]

                print_log_summary(logs, subject_name=name, file=out)
                simulate_subject_from_logs(
                    name        = name,
                    oral        = oral,
                    daily_logs  = logs,
                    age         = age,
                    fitzpatrick = fitz,
                    C0          = C0,
                    file        = out,
                )

        except ImportError:
            print("[INFO] data_loader not available — legacy mode only.\n", file=out)

        # =====================================================================
        # PATH B — Legacy mode (uncomment blocks below to use)
        # =====================================================================
        uvb_all = [300.0] * 6

        # oral_oscar  = load_oral_doses_from_json(ORAL_JSON, "Oscar",   dates_oscar)
        # oral_nicole = load_oral_doses_from_json(ORAL_JSON, "Nicole",  dates_nicole)
        # oral_vi     = load_oral_doses_from_json(ORAL_JSON, "vi_pdf",  dates_vi)

        # simulate_subject(
        #     name          = "Subject A -- 21 y/o, Fitzpatrick II  [legacy]",
        #     oral          = oral_oscar,
        #     uvb_raw_j_m2  = uvb_all,
        #     body_area_pct = [9.0, 9.0, 0.0, 9.0, 0.0, 9.0],
        #     age           = 21,
        #     fitzpatrick   = 2,
        #     C0            = 50.0,
        #     file          = out,
        # )
        # simulate_subject(
        #     name          = "Subject B -- 20 y/o, Fitzpatrick II  [legacy]",
        #     oral          = oral_nicole,
        #     uvb_raw_j_m2  = uvb_all,
        #     body_area_pct = [9.0, 0.0, 9.0, 9.0, 0.0, 9.0],
        #     age           = 20,
        #     fitzpatrick   = 2,
        #     C0            = 50.0,
        #     file          = out,
        # )
        # simulate_subject(
        #     name          = "Subject C -- 22 y/o, Fitzpatrick III  [legacy]",
        #     oral          = oral_vi,
        #     uvb_raw_j_m2  = uvb_all,
        #     body_area_pct = [0.0, 9.0, 9.0, 0.0, 15.0, 9.0],
        #     age           = 22,
        #     fitzpatrick   = 3,
        #     C0            = 50.0,
        #     file          = out,
        # )

    # -------------------------------------------------------------------------
    # Write to file AND echo to terminal simultaneously
    # -------------------------------------------------------------------------
    _real_stdout = sys.stdout   # capture before any redirection

    with open(OUTPUT_FILE, "w", encoding="utf-8") as outfile:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header = f"vitamin_d_model.py  —  results generated {timestamp}\n"
        _real_stdout.write(header + "\n")
        outfile.write(header + "\n")

        class Tee:
            """Write to both the output file and the real terminal stdout."""
            def write(self, msg):
                _real_stdout.write(msg)
                outfile.write(msg)
            def flush(self):
                _real_stdout.flush()
                outfile.flush()

        with contextlib.redirect_stdout(Tee()):
            run_all(out=sys.stdout)

    _real_stdout.write(f"\nResults saved to: {OUTPUT_FILE}\n")
    