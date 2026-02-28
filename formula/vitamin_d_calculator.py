"""
=============================================================================
VITAMIN D DEFICIENCY CALCULATOR
=============================================================================
Source of equations:
    Diffey, B.L. (2013). "Modelling vitamin D status due to oral intake and
    sun exposure in an adult British population."
    British Journal of Nutrition, 110(3), pp.569–577.
    DOI: https://doi.org/10.1017/S0007114512005466

Equations used from paper:
    Eq. 1  — Total plasma 25(OH)D on day T (master equation)
    Eq. 2  — Oral intake contribution to blood level
    Eq. 3  — Sun/UV contribution to blood level
    Eq. 4  — Diminishing returns factor F(T) for UV response
    Eq. 7  — UV response function R_UV(t)
    Eq. 9  — Oral intake response function R_oral(t)

Equations NOT used (and why):
    Eq. 5 & 6 — Skipped: those estimate UV from behaviour. We have real UV data.
    Eq. 10    — Skipped: that is for planning fixed supplements, not tracking.
    Eq. 11    — Skipped: that is for population stats. We only have 1 person.

=============================================================================
FIXED BIOLOGICAL PARAMETERS (from Diffey 2013, do not change these)
=============================================================================
    f   = 0.15       fraction of vitamin D that goes into tissue storage
    β   = 25  days   half-life of 25(OH)D in plasma (blood clearance)
    γ   = 250 days   half-life of 25(OH)D in tissue stores (slow release)
    α   = 0.6 days   half-time for UV-derived vitamin D to enter blood
    α'  = 1.5 days   half-time for oral vitamin D to enter blood (gut slower)
    A   = 0.18       UV scaling factor: nmol/L per SED per 1% body surface
    S   = 0.023      oral scaling factor: nmol/L per µg of vitamin D

=============================================================================
EXPECTED CSV FORMAT
=============================================================================
Your CSV must have these columns (column names are case-insensitive):

    date             — YYYY-MM-DD
    uvi              — UV Index for the day (used to compute SED)
    uva              — UVA irradiance (informational only, not used in model*)
    uvb              — UVB irradiance (informational only, not used in model*)
    oral_intake_ug   — vitamin D intake that day in micrograms (µg)
                       1000 IU pill = 25 µg | 400 IU pill = 10 µg
    skin_area_pct    — % of body surface area exposed to sun (0–100)
                       face+hands only ≈ 8%
                       + forearms      ≈ 15%
                       + legs          ≈ 30%

* NOTE on UVA vs UVB vs UVI:
    - UVA (315–400 nm) does NOT synthesise vitamin D in skin. Not used.
    - UVB (280–315 nm) IS responsible for vitamin D synthesis.
    - UVI (UV Index) is an erythemal-weighted measure that correlates well
      with UVB and is the standard metric used in the Diffey model.
      We convert UVI → SED (Standard Erythema Dose) to match the paper.

=============================================================================
"""

import pandas as pd
import numpy as np
import sys
import os


# =============================================================================
# FIXED BIOLOGICAL CONSTANTS (Diffey 2013, Table in Methods section)
# =============================================================================

f     = 0.15    # fraction of synthesised/absorbed vitamin D stored in tissue
beta  = 25.0    # β: plasma clearance half-life (days)
gamma = 250.0   # γ: tissue store clearance half-life (days)
alpha = 0.6     # α: half-time for UV vitamin D uptake into plasma (days)
alpha_prime = 1.5  # α': half-time for oral vitamin D uptake into plasma (days)
A     = 0.18    # UV scaling factor (nmol/L per SED per 1% BSA)
S     = 0.023   # oral scaling factor (nmol/L per µg)

# Deficiency thresholds (clinical standard, referenced in Diffey 2013)
THRESHOLD_DEFICIENT   = 25.0   # nmol/L — below this = DEFICIENT
THRESHOLD_INSUFFICIENT = 50.0  # nmol/L — below this = INSUFFICIENT


# =============================================================================
# UV INDEX → SED CONVERSION
# =============================================================================

def uvi_to_sed(uvi_value, daylight_hours=12.0):
    """
    Convert UV Index to SED (Standard Erythema Dose).

    SED is the UV unit used in the Diffey model (Eq. 3, 7).
    1 SED = 100 J/m² of erythemal-weighted UV radiation.

    Method:
        UVI × 0.025 = erythemal irradiance in W/m²
        Assuming a triangular distribution of UV over the day peaking at noon:
            daily energy = peak_irradiance × daylight_seconds / 2
        Then divide by 100 to convert J/m² → SED.

    Parameters:
        uvi_value      : float  — the UV Index reading for the day
        daylight_hours : float  — hours of daylight (default 12, adjust seasonally)

    Returns:
        float — daily SED value
    """
    if uvi_value <= 0:
        return 0.0
    peak_irradiance_W_per_m2 = uvi_value * 0.025   # W/m² erythemal
    daylight_seconds = daylight_hours * 3600
    daily_J_per_m2 = peak_irradiance_W_per_m2 * daylight_seconds / 2  # triangle
    sed = daily_J_per_m2 / 100.0
    return sed


# =============================================================================
# RESPONSE FUNCTIONS (Eq. 7 and Eq. 9 from Diffey 2013)
# =============================================================================

def R_UV(t):
    """
    UV Response Function — Eq. 7 from Diffey (2013).

    Represents the increase in plasma 25(OH)D (nmol/L) at time t days
    after a single UV exposure of 1 SED to 1% of body surface area.

    This is a two-compartment model:
        - Fast compartment: (1-f) of vitamin D stays in plasma, clears with half-life β
        - Slow compartment: f of vitamin D goes to tissue stores, clears with half-life γ
        - Uptake delay: −2^(−t/α) subtracts until vitamin D has entered blood

    R_UV(t) = A × [ (1−f)×2^(−t/β) + f×2^(−t/γ) − 2^(−t/α) ]

    Parameters:
        t : int/float — number of days since the UV exposure

    Returns:
        float — nmol/L increase per SED per 1% BSA
    """
    if t <= 0:
        return 0.0
    fast_compartment  = (1 - f) * (2 ** (-t / beta))    # plasma clearance
    slow_compartment  = f       * (2 ** (-t / gamma))   # tissue store release
    uptake_delay      =            2 ** (-t / alpha)     # absorption lag
    return A * (fast_compartment + slow_compartment - uptake_delay)


def R_oral(t):
    """
    Oral Intake Response Function — Eq. 9 from Diffey (2013).

    Represents the increase in plasma 25(OH)D (nmol/L) at time t days
    after a single oral dose of 1 µg of vitamin D (cholecalciferol).

    Same structure as R_UV but:
        - Scaled by S instead of A
        - Uptake delay uses α' (1.5 days) instead of α (0.6 days)
          because gut absorption is slower than skin synthesis

    R_oral(t) = S × [ (1−f)×2^(−t/β) + f×2^(−t/γ) − 2^(−t/α') ]

    Parameters:
        t : int/float — number of days since the oral dose

    Returns:
        float — nmol/L increase per µg of vitamin D
    """
    if t <= 0:
        return 0.0
    fast_compartment  = (1 - f) * (2 ** (-t / beta))
    slow_compartment  = f       * (2 ** (-t / gamma))
    uptake_delay      =            2 ** (-t / alpha_prime)
    return S * (fast_compartment + slow_compartment - uptake_delay)


# =============================================================================
# MAIN MODEL — runs day by day through the CSV
# =============================================================================

def run_vitamin_d_model(df):
    """
    Simulates daily plasma 25(OH)D for one person using Diffey (2013).

    Implements:
        Eq. 2: C_oral(T)  = sum over all past days of O(t) × R_oral(T−t+1)
        Eq. 3: C_sun(T)   = sum over all past days of E(t) × A(t) × R_UV(T−t+1)
        Eq. 4: F(T)       = exp(−0.01 × C_total(T−1))
        Eq. 1: C_total(T) = C_total(T−1) + ΔC_oral(T) + F(T) × ΔC_sun(T)

    Variables per day t:
        O(t)  = oral_intake_ug  — daily vitamin D intake in µg
        E(t)  = SED             — UV dose in Standard Erythema Doses
        A(t)  = skin_area_pct   — % body surface area exposed (0–100)

    Parameters:
        df : pandas DataFrame — must have columns: sed, oral_intake_ug, skin_area_pct

    Returns:
        list of float — C_total for each day
    """
    n_days = len(df)

    # Pre-compute response function values for all possible time lags
    # (avoids recomputing inside the inner loop)
    max_lag = n_days + 1
    r_uv_cache   = np.array([R_UV(t)   for t in range(max_lag + 1)])
    r_oral_cache = np.array([R_oral(t) for t in range(max_lag + 1)])

    C_total = np.zeros(n_days)  # total plasma 25(OH)D (nmol/L) each day
    C_sun   = np.zeros(n_days)  # sun contribution each day
    C_oral  = np.zeros(n_days)  # oral contribution each day

    for T in range(n_days):

        # --- Eq. 2: oral contribution on day T ---
        # Sum all past oral doses weighted by how much effect they still have today
        oral_sum = 0.0
        for t in range(T + 1):
            lag = T - t + 1          # how many days ago was dose t
            O_t = df['oral_intake_ug'].iloc[t]   # O(t): dose on day t in µg
            oral_sum += O_t * r_oral_cache[lag]
        C_oral[T] = oral_sum

        # --- Eq. 3: sun contribution on day T ---
        # Sum all past UV exposures weighted by skin area and UV response
        sun_sum = 0.0
        for t in range(T + 1):
            lag = T - t + 1
            E_t = df['sed'].iloc[t]              # E(t): UV dose in SED
            A_t = df['skin_area_pct'].iloc[t]    # A(t): % body surface exposed
            sun_sum += E_t * A_t * r_uv_cache[lag]
        C_sun[T] = sun_sum

        # --- Eq. 4: diminishing returns factor F(T) ---
        # Higher existing blood level → less benefit from additional UV
        C_prev = C_total[T - 1] if T > 0 else 0.0
        F_T = np.exp(-0.01 * C_prev)

        # --- Eq. 1: total plasma 25(OH)D on day T ---
        delta_oral = C_oral[T] - (C_oral[T-1] if T > 0 else 0.0)
        delta_sun  = C_sun[T]  - (C_sun[T-1]  if T > 0 else 0.0)
        C_total[T] = C_prev + delta_oral + F_T * delta_sun

        # Clamp to 0 (blood levels can't be negative)
        C_total[T] = max(0.0, C_total[T])

    return C_total


# =============================================================================
# DEFICIENCY ASSESSMENT — single person output
# =============================================================================

def assess_deficiency(C_total_values, dates):
    """
    For each day, print whether the person is deficient, insufficient,
    or sufficient — and by how much.

    Thresholds (Diffey 2013, Pearce & Cheetham cited within):
        < 25 nmol/L  → DEFICIENT
        < 50 nmol/L  → INSUFFICIENT
        ≥ 50 nmol/L  → SUFFICIENT
    """
    print("\n" + "="*65)
    print("  VITAMIN D STATUS REPORT")
    print("="*65)
    print(f"  {'Date':<14} {'25(OH)D (nmol/L)':>18}  {'Status':<15} {'Gap'}")
    print("-"*65)

    for i, (date, level) in enumerate(zip(dates, C_total_values)):
        if level < THRESHOLD_DEFICIENT:
            status = "⚠️  DEFICIENT"
            gap = f"{THRESHOLD_DEFICIENT - level:.1f} nmol/L below deficiency threshold"
        elif level < THRESHOLD_INSUFFICIENT:
            status = "⚡ INSUFFICIENT"
            gap = f"{THRESHOLD_INSUFFICIENT - level:.1f} nmol/L below sufficient level"
        else:
            status = "✅ SUFFICIENT"
            gap = f"{level - THRESHOLD_INSUFFICIENT:.1f} nmol/L above sufficient level"

        # Print every day, or just print summary every 30 days to avoid huge output
        # Change to: if True  — to print every single day
        if i % 30 == 0 or i == len(C_total_values) - 1:
            print(f"  {str(date):<14} {level:>18.2f}  {status:<20} {gap}")

    # --- Final summary ---
    final_level = C_total_values[-1]
    n_deficient    = sum(1 for v in C_total_values if v < THRESHOLD_DEFICIENT)
    n_insufficient = sum(1 for v in C_total_values if THRESHOLD_DEFICIENT <= v < THRESHOLD_INSUFFICIENT)
    n_sufficient   = sum(1 for v in C_total_values if v >= THRESHOLD_INSUFFICIENT)

    print("="*65)
    print(f"\n  FINAL DAY LEVEL : {final_level:.2f} nmol/L")
    print()

    if final_level < THRESHOLD_DEFICIENT:
        shortfall = THRESHOLD_DEFICIENT - final_level
        print(f"  STATUS : ⚠️  DEFICIENT")
        print(f"  You are {shortfall:.1f} nmol/L BELOW the deficiency threshold of {THRESHOLD_DEFICIENT} nmol/L.")
        print(f"  You would need to raise levels by {shortfall:.1f} nmol/L to exit deficiency.")
    elif final_level < THRESHOLD_INSUFFICIENT:
        shortfall = THRESHOLD_INSUFFICIENT - final_level
        print(f"  STATUS : ⚡ INSUFFICIENT")
        print(f"  You are {shortfall:.1f} nmol/L BELOW the sufficiency threshold of {THRESHOLD_INSUFFICIENT} nmol/L.")
        print(f"  You are above deficiency but not yet at optimal levels.")
    else:
        surplus = final_level - THRESHOLD_INSUFFICIENT
        print(f"  STATUS : ✅ SUFFICIENT")
        print(f"  You are {surplus:.1f} nmol/L ABOVE the sufficiency threshold of {THRESHOLD_INSUFFICIENT} nmol/L.")

    print()
    print(f"  Days deficient    : {n_deficient}  ({100*n_deficient/len(C_total_values):.1f}% of tracked period)")
    print(f"  Days insufficient : {n_insufficient}  ({100*n_insufficient/len(C_total_values):.1f}% of tracked period)")
    print(f"  Days sufficient   : {n_sufficient}  ({100*n_sufficient/len(C_total_values):.1f}% of tracked period)")
    print("="*65)


# =============================================================================
# LOAD CSV AND RUN
# =============================================================================

def load_and_validate_csv(filepath):
    """
    Loads the CSV and normalises column names.
    Expected columns: date, uvi, uva, uvb, oral_intake_ug, skin_area_pct
    """
    df = pd.read_csv(filepath)
    df.columns = df.columns.str.strip().str.lower()

    required = ['date', 'uvi', 'oral_intake_ug', 'skin_area_pct']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {missing}\n"
            f"Found columns: {list(df.columns)}\n\n"
            f"Required: date, uvi, oral_intake_ug, skin_area_pct\n"
            f"Optional: uva, uvb (noted but not used in Diffey model)\n\n"
            f"Note: oral_intake_ug is daily vitamin D in micrograms\n"
            f"      1000 IU pill = 25 µg  |  400 IU pill = 10 µg\n"
            f"      skin_area_pct: % body exposed (face+hands=8, +arms=15, +legs=30)"
        )

    df['date'] = pd.to_datetime(df['date'])
    df = df.sort_values('date').reset_index(drop=True)

    # Note UVA/UVB in output but explain why they aren't used
    if 'uva' in df.columns or 'uvb' in df.columns:
        print("\n  NOTE on UV columns:")
        print("  ├── UVA : detected in CSV — NOT used in model.")
        print("  │         UVA (315–400 nm) does not synthesise vitamin D in skin.")
        print("  ├── UVB : detected in CSV — NOT directly used in model.")
        print("  │         UVB (280–315 nm) drives vitamin D but the Diffey model")
        print("  │         uses SED (erythemal dose) as its UV metric, which we")
        print("  │         derive from UVI since it correlates best with UVB activity.")
        print("  └── UVI : USED — converted to SED (Standard Erythema Doses).\n")

    # Convert UVI → SED
    df['sed'] = df['uvi'].apply(uvi_to_sed)

    return df


def main(csv_path):
    print(f"\n  Loading: {csv_path}")
    df = load_and_validate_csv(csv_path)
    print(f"  Loaded {len(df)} days of data ({df['date'].iloc[0].date()} → {df['date'].iloc[-1].date()})")

    print("  Running Diffey (2013) vitamin D model...")
    C_total = run_vitamin_d_model(df)

    assess_deficiency(C_total, df['date'].dt.date)

    # Also save results to CSV
    df['C_total_nmol_L'] = C_total
    out_path = csv_path.replace('.csv', '_results.csv')
    df[['date', 'uvi', 'sed', 'oral_intake_ug', 'skin_area_pct', 'C_total_nmol_L']].to_csv(out_path, index=False)
    print(f"\n  Full daily results saved to: {out_path}\n")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("\n  Usage:  python vitamin_d_calculator.py your_data.csv")
        print("\n  Example CSV format:")
        print("  date,uvi,uva,uvb,oral_intake_ug,skin_area_pct")
        print("  2024-01-01,2.1,3500,45,25,8")
        print("  2024-01-02,1.8,3200,38,25,8")
        print("  ...")
        print("\n  oral_intake_ug: 1000 IU pill = 25 µg | 400 IU pill = 10 µg")
        print("  skin_area_pct : face+hands=8 | +forearms=15 | +legs=30\n")
        sys.exit(1)

    main(sys.argv[1])
