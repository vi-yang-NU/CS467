"""
vitamin_d_model.py
==================
A computational implementation of the plasma 25(OH)D (vitamin D) model
described in the accompanying LaTeX document, based on the pharmacokinetic
framework of Diffey (2013), with demographic corrections from:
  - Hilger et al. (2014) — sex independence
  - Young (2020)         — Fitzpatrick skin-type scaling
  - Chalcraft et al. (2020) — age-related synthesis decline

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
    achieve equivalent synthesis compared to Types I-II, but that Types III-VI
    are not statistically distinguishable from each other.

    The factor is inverted here because we want to *scale down* the effective
    UV dose for darker skin types, with Type I-II as the baseline (= 1.0).

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
    capacity, anchored at age 20 (f_age = 1.0).  By the seventies this
    amounts to roughly a 50% reduction relative to a 20-year-old under
    identical UV conditions.

    Note: this function can return negative values for very old ages
    (>~97 years).  Callers should clamp to [0, 1] if needed.

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
# UV DOSE (SED) CALCULATOR
# =============================================================================
#
# NOTE ON INPUT DATA -- UVB vs UVI vs UVA
# ----------------------------------------
# The model (Eq. 1) was originally written in terms of UV Index (UVI), but
# UVI is itself *derived* from UVB: it is the erythemally-weighted UVB
# irradiance divided by a reference value of 0.025 W/m2.  In other words:
#
#     UVI = E_eff (W/m2) / 0.025
#
# If your data source provides total daily UVB dose in J/m2, you already have
# E_eff integrated over the day, so no daylight-hours parameter is needed.
# The conversion to SED is then simply:
#
#     SED = UVB_dose (J/m2) / 100
#
# because by definition 1 SED = 100 J/m2 of erythemally effective UV.
#
# UVA (315-400 nm) plays NO meaningful role in vitamin D synthesis.
# The photolysis of 7-dehydrocholesterol to pre-vitamin D3 requires photons
# in the UVB range (290-315 nm).  UVA can therefore be discarded entirely.

def uvb_dose_to_sed(uvb_j_per_m2: float) -> float:
    """
    Eq. 1 (adapted)  --  Convert daily UVB dose to Standard Erythemal Doses.

    Replaces the UVI-based formulation with a direct conversion from measured
    daily UVB irradiance (J/m2).  The two representations are equivalent
    because UVI is defined as erythemally-weighted UVB irradiance / 0.025 W/m2,
    so integrating over a day and dividing by 100 J/m2 per SED yields the same
    result as Eq. 1 in the source document.

    Conversion derivation:
        1 SED  =  100 J/m2  of erythemally effective UV  (ISO 17166 definition)
        therefore:  SED  =  UVB_dose (J/m2) / 100

    UVA is excluded because it does not drive vitamin D photosynthesis; only
    UVB photons (290-315 nm) have sufficient energy to cleave the B-ring of
    7-dehydrocholesterol and initiate the vitamin D3 cascade.

    Parameters
    ----------
    uvb_j_per_m2 : float
        Total daily UVB irradiance dose in J/m2.

    Returns
    -------
    float
        E(t) in SED.
    """
    return uvb_j_per_m2 / 100.0


# =============================================================================
# RESPONSE FUNCTIONS  (impulse responses / Green's functions)
# =============================================================================

def R_UV(t: float) -> float:
    """
    Eq. 2  --  UV-derived vitamin D impulse response function.

    Describes the rise-and-fall kinetics of plasma 25(OH)D following a single
    unit UV exposure (1 SED over 1% body area).  The two-exponential decay
    captures redistribution into a slow tissue compartment (fraction f,
    half-life gamma) alongside the faster plasma compartment (fraction 1-f,
    half-life beta).  The negative term (half-time alpha) models the initial
    uptake lag before pre-vitamin D3 is converted to vitamin D3 and appears
    in plasma.

    Parameters
    ----------
    t : float
        Time since exposure (days).  Must be >= 0.

    Returns
    -------
    float
        Response in nmol/L per (SED * % body area).
    """
    return A_uv * (
        (1 - f) * 2 ** (-t / beta)       # fast plasma compartment
        + f     * 2 ** (-t / gamma)      # slow tissue compartment
        - 2     ** (-t / alpha)          # uptake lag (negative)
    )


def R_oral(t: float) -> float:
    """
    Eq. 3  --  Oral vitamin D impulse response function.

    Analogous to R_UV but for an ingested dose.  The uptake lag uses alpha'
    (1.5 days) rather than alpha (0.6 days), reflecting the slower
    gastrointestinal absorption pathway compared with cutaneous synthesis.
    The scaling factor S replaces A_uv to convert from ug to nmol/L.

    Parameters
    ----------
    t : float
        Time since ingestion (days).  Must be >= 0.

    Returns
    -------
    float
        Response in nmol/L per ug.
    """
    return S * (
        (1 - f) * 2 ** (-t / beta)       # fast plasma compartment
        + f     * 2 ** (-t / gamma)      # slow tissue compartment
        - 2     ** (-t / alpha_)         # uptake lag (slower for oral route)
    )


# =============================================================================
# CONVOLUTION ACCUMULATORS
# (running cumulative plasma contributions up to day T)
# =============================================================================

def compute_C_oral(oral_doses: list) -> list:
    """
    Eq. 4  --  Cumulative oral contribution to plasma 25(OH)D.

    Computes the discrete convolution of daily oral intake O(t) with the
    oral impulse response R_oral.  Each day's dose is spread forward in
    time according to R_oral, and contributions from all past doses are
    summed.

        C_oral(T) = sum_{t=0}^{T}  O(t) * R_oral(T - t + 1)

    The "+1" offset in the lag argument means the response on the same day
    as ingestion uses R_oral(1), i.e. one day of processing has elapsed.

    Parameters
    ----------
    oral_doses : list
        O(t) for t = 0, 1, ..., N-1  (ug per day).

    Returns
    -------
    list
        C_oral(T) for T = 0, 1, ..., N-1  (nmol/L).
    """
    N = len(oral_doses)
    C_oral = [0.0] * N
    for T in range(N):
        total = 0.0
        for t in range(T + 1):
            lag = T - t + 1          # days since dose t was taken (>= 1)
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
    Eq. 11  --  Cumulative UV contribution to plasma 25(OH)D
                (replaces Eq. 6 from an earlier model revision).

    Extends the basic sun convolution (Eq. 6) by incorporating the two
    demographic scalars f_age and f_skin directly into the synthesis term,
    so that cutaneous production is modulated by age-related decline and
    melanin-dependent UV attenuation.

        C_sun(T) = sum_{t=0}^{T}  f_age * f_skin * E(t) * A(t) * R_UV(T - t + 1)

    Both f_age and f_skin are time-invariant for a given subject, so they
    are computed once and factored through the sum.

    Parameters
    ----------
    uv_doses   : list  -- E(t) in SED for each day.
    body_areas : list  -- A(t) in % body surface area exposed each day.
    age        : float        -- Subject age in years (for f_age).
    fitzpatrick: int          -- Fitzpatrick skin type 1-6 (for f_skin).

    Returns
    -------
    list
        C_sun(T) for T = 0, 1, ..., N-1  (nmol/L).
    """
    f_age_val  = age_factor(age)
    f_skin_val = skin_factor(fitzpatrick)
    demo_scalar = f_age_val * f_skin_val      # combined demographic multiplier

    N = len(uv_doses)
    C_sun = [0.0] * N
    for T in range(N):
        total = 0.0
        for t in range(T + 1):
            lag = T - t + 1
            total += demo_scalar * uv_doses[t] * body_areas[t] * R_UV(lag)
        C_sun[T] = total
    return C_sun


# =============================================================================
# SATURATION FACTOR
# =============================================================================

def saturation_factor(C_total_prev: float) -> float:
    """
    Eq. 5  --  Diminishing-returns (saturation) factor F(T).

    At high circulating 25(OH)D concentrations the body downregulates
    further synthesis and absorption.  This is modelled as an exponential
    decay of the UV contribution's incremental effect, parameterised so
    that F ~= 1 at low/normal levels and decreases smoothly as levels rise.

        F(T) = exp( -0.01 * C_total(T-1) )

    Parameters
    ----------
    C_total_prev : float
        Plasma 25(OH)D on the *previous* day, C_total(T-1)  (nmol/L).

    Returns
    -------
    float
        F(T) -- dimensionless, in (0, 1].
    """
    return math.exp(-0.01 * C_total_prev)


# =============================================================================
# MASTER MODEL -- total plasma 25(OH)D
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
    Eq. 7  --  Master equation: total plasma 25(OH)D on day T.

    Integrates oral and UV contributions day-by-day, applying the saturation
    factor F(T) to the incremental UV-derived change at each step:

        C_total(T) = C_total(T-1)
                     + [C_oral(T)  - C_oral(T-1)]
                     + F(T) * [C_sun(T) - C_sun(T-1)]

    The incremental formulation (delta rather than absolute values) means the
    saturation factor only suppresses *new* UV-derived input, not the
    carry-over of previously accumulated levels.

    Parameters
    ----------
    oral_doses   : list  -- O(t) in ug for each simulation day.
    uv_doses     : list  -- E(t) in SED for each simulation day.
    body_areas   : list  -- A(t) in % BSA for each simulation day.
    age          : float        -- Subject age in years.
    fitzpatrick  : int          -- Fitzpatrick skin type 1-6.
    C0           : float        -- Initial plasma 25(OH)D level (nmol/L).

    Returns
    -------
    list
        C_total(T) for T = 0, 1, ..., N-1  (nmol/L).
    """
    N = len(oral_doses)
    assert len(uv_doses)   == N, "All input lists must be the same length."
    assert len(body_areas) == N, "All input lists must be the same length."

    # Pre-compute full oral and UV contribution curves (Eqs. 4 and 11)
    C_oral = compute_C_oral(oral_doses)
    C_sun  = compute_C_sun(uv_doses, body_areas, age, fitzpatrick)

    C_total = [0.0] * N
    for T in range(N):
        if T == 0:
            C_prev      = C0
            C_oral_prev = 0.0
            C_sun_prev  = 0.0
        else:
            C_prev      = C_total[T - 1]
            C_oral_prev = C_oral[T - 1]
            C_sun_prev  = C_sun[T - 1]

        # Incremental contributions on day T
        delta_oral = C_oral[T] - C_oral_prev
        delta_sun  = C_sun[T]  - C_sun_prev

        # Saturation factor based on previous day's level (Eq. 5)
        F = saturation_factor(C_prev)

        # Master update (Eq. 7)
        C_total[T] = C_prev + delta_oral + F * delta_sun

    return C_total


# =============================================================================
# SIMULATION RUNNER
# =============================================================================

def simulate_subject(
    name:          str,
    oral:          list,
    uvb_raw_j_m2:  list,
    body_area_pct: list,
    age:           float,
    fitzpatrick:   int,
    C0:            float,
) -> list:
    """
    Convenience wrapper that accepts raw inputs, converts units, runs the
    model, prints a formatted results table, and returns the daily plasma
    25(OH)D trajectory.

    This function handles the unit-conversion step between raw measured data
    and the model's internal representation:
        UVB (J/m2)  ->  SED    via uvb_dose_to_sed()

    Body surface area is a per-day list because clothing and activity level
    vary from day to day (e.g. indoors vs outdoors, light vs heavy clothing).

    Parameters
    ----------
    name           : str          -- Label used in the printed header.
    oral           : list  -- Daily oral vitamin D intake (ug/day).
                                     1 ug = 40 IU.
    uvb_raw_j_m2   : list  -- Daily UVB dose (J/m2).  UVA excluded.
    body_area_pct  : list  -- Body surface area exposed per day (%).
                                     Must be the same length as oral.
                                     Reference values per day:
                                       ~15%  dressed outdoors (face/hands/arms)
                                       ~35%  t-shirt + shorts
                                       ~80%  swimwear
    age            : float        -- Subject age in years.
    fitzpatrick    : int          -- Fitzpatrick skin type (1-6).
    C0             : float        -- Initial plasma 25(OH)D (nmol/L).
                                     WARNING: assumed if no blood test available.

    Returns
    -------
    list
        C_total(T) for each day T  (nmol/L).
    """
    assert len(oral) == len(uvb_raw_j_m2), (
        f"[{name}] oral and uvb_raw_j_m2 must be the same length "
        f"(got {len(oral)} vs {len(uvb_raw_j_m2)})"
    )
    assert len(oral) == len(body_area_pct), (
        f"[{name}] oral and body_area_pct must be the same length "
        f"(got {len(oral)} vs {len(body_area_pct)})"
    )

    N = len(oral)

    # Convert UVB J/m2 -> SED for each day
    uv_doses = [uvb_dose_to_sed(uvb) for uvb in uvb_raw_j_m2]

    # Run the full pharmacokinetic model (Eqs. 4, 5, 7, 11)
    result = run_model(
        oral_doses  = oral,
        uv_doses    = uv_doses,
        body_areas  = body_area_pct,
        age         = age,
        fitzpatrick = fitzpatrick,
        C0          = C0,
    )

    # ------------------------------------------------------------------
    # Print results table
    # ------------------------------------------------------------------
    print("=" * 62)
    print(f"  {name}")
    print(f"    Age              : {age} years")
    print(f"    Fitzpatrick type : {fitzpatrick}")
    print(f"    f_age            : {age_factor(age):.3f}")
    print(f"    f_skin           : {skin_factor(fitzpatrick):.3f}")
    print(f"    Initial 25(OH)D  : {C0:.1f} nmol/L  (WARNING: assumed)")
    print("=" * 62)
    print(f"  {'Day':>4}  {'Oral (ug)':>10}  {'UVB (J/m2)':>12}  {'BSA (%)':>8}  {'C_total (nmol/L)':>18}")
    print("  " + "-" * 57)
    for day in range(N):
        print(
            f"  {day:>4}  {oral[day]:>10.1f}  {uvb_raw_j_m2[day]:>12.1f}"
            f"  {body_area_pct[day]:>8.1f}  {result[day]:>18.2f}"
        )
    print()

    return result


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # =========================================================================
    # SHARED INPUT DATA
    # Both oral intake and UVB are the same measured data for all subjects --
    # they were exposed to the same environment and diet over these 6 days.
    # Replace each value with the real measured observations.
    #
    # Oral intake units : ug/day  (1 ug = 40 IU)
    # UVB units         : J/m2   (UVA excluded -- irrelevant to vitamin D)
    # =========================================================================

    oral_all = [
        15.0,   # day 0  -- replace with real dietary + supplement total (ug)
        15.0,   # day 1
        15.0,   # day 2
        15.0,   # day 3
        15.0,   # day 4
        15.0,   # day 5
    ]

    uvb_all = [
        300.0,  # day 0  -- replace with real measured UVB dose (J/m2)
        300.0,  # day 1
        300.0,  # day 2
        300.0,  # day 3
        300.0,  # day 4
        300.0,  # day 5
    ]

    # =========================================================================
    # SUBJECT PROFILES
    # Subjects share the same oral and UVB environment but differ in age,
    # Fitzpatrick skin type, and daily clothing / activity (body area).
    #
    # body_area_pct is a per-day list -- update each entry to reflect what
    # the subject was actually wearing or doing on that day.
    #
    # C0 notes:
    #   0.0  -- treat baseline as unknown (conservative lower bound)
    #  50.0  -- approximate population mean for young adults
    #  75.0  -- mid-range sufficiency threshold
    # =========================================================================

    simulate_subject(
        name          = "Subject A -- 21 y/o, Fitzpatrick II",
        oral          = oral_all,
        uvb_raw_j_m2  = uvb_all,
        body_area_pct = [
            9.0,   # day 0  -- replace with real daily BSA (%)
            9.0,   # day 1
            0.0,   # day 2
            9.0,   # day 3
            0.0,   # day 4
            9.0,   # day 5
        ],
        age           = 21,
        fitzpatrick   = 2,
        C0            = 50.0,
    )

    simulate_subject(
        name          = "Subject B -- 20 y/o, Fitzpatrick II",
        oral          = oral_all,
        uvb_raw_j_m2  = uvb_all,
        body_area_pct = [
            9.0,   # day 0  -- replace with real daily BSA (%)
            0.0,   # day 1
            9.0,   # day 2
            9.0,   # day 3
            0.0,   # day 4
            9.0,   # day 5
        ],
        age           = 20,
        fitzpatrick   = 2,
        C0            = 50.0,
    )

    simulate_subject(
        name          = "Subject C -- 22 y/o, Fitzpatrick III",
        oral          = oral_all,
        uvb_raw_j_m2  = uvb_all,
        body_area_pct = [
            0.0,   # day 0  -- replace with real daily BSA (%)
            9.0,   # day 1
            9.0,   # day 2
            0.0,   # day 3
            15.0,   # day 4
            9.0,   # day 5
        ],
        age           = 22,
        fitzpatrick   = 3,
        C0            = 50.0,
    )