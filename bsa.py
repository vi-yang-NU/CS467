"""
bsa.py
======
Body Surface Area (BSA) utilities for use with vitamin_d_model.py.

Two operational modes are supported:

  1. SURVEY MODE  (original)
     Binary per-day survey of which body regions were exposed.
     Used when no wearable sensor data is available.
     Returns body_area_pct  — a % BSA value per day.

  2. SENSOR MODE  (new — for CSV log data)
     Per-body-part SED readings from the wearable UV sensor CSV logs.
     Returns uv_eff_doses  — the dose-area product Σ_i(SED_i × BSA_fraction_i)
     per day, ready to pass directly to run_model_from_effective_doses() in
     vitamin_d_model.py.  This mode is more accurate than the survey approach
     because it uses *measured* exposure at each body site rather than a
     binary clothed/unclothed approximation.

BSA values are derived from the Lund-Browder chart (standard clinical
reference for body surface area proportions).

Author Credit : Claude (Anthropic, claude-sonnet-4-6)
"""


# =============================================================================
# BSA CONTRIBUTIONS PER BODY REGION  (% of total body surface area)
# Source: Lund-Browder chart, adapted for sun-exposure relevance.
# The keys align with the CSV sensor column names used by data_loader.py.
# =============================================================================

BSA_REGIONS = {
    # ---- head / face / neck ----
    "sed_scalp":     1.5,   # Top of head (scalp only; hair usually blocks some UV)
    "sed_forehead":  1.5,   # Upper face
    "sed_ears":      0.5,   # Auricles (both)
    "sed_jaw":       0.5,   # Lower face / chin
    "sed_neck":      1.0,   # Anterior + posterior neck

    # ---- torso ----
    "sed_chest":     9.0,   # Upper anterior trunk
    "sed_abdomen":   9.0,   # Lower anterior trunk
    "sed_back":     18.0,   # Entire posterior trunk

    # ---- arms ----
    "sed_shoulders": 3.0,   # Both shoulder caps
    "sed_upperArms": 6.0,   # Both upper arms (shoulder to elbow)
    "sed_forearms":  6.0,   # Both forearms (elbow to wrist)
    "sed_hands":     5.0,   # Both hands (dorsal + palmar)

    # ---- legs ----
    "sed_thighs":    9.5,   # Both anterior thighs
    "sed_knees":     1.0,   # Both knee caps
    "sed_calves":    7.0,   # Both calves
    "sed_ankles":    1.5,   # Both ankle regions
    "sed_feet":      7.0,   # Both feet
}

# Survey-mode convenience lookup (aggregated regions for binary surveys).
# Keys match the human-readable survey parameters in compute_bsa_from_survey().
_SURVEY_BSA = {
    "head":         BSA_REGIONS["sed_forehead"] + BSA_REGIONS["sed_jaw"],   # 2.0 — face/neck; scalp excluded
    "hands":        BSA_REGIONS["sed_hands"],       # 5.0
    "forearms":     BSA_REGIONS["sed_forearms"],    # 6.0
    "full_arms":    BSA_REGIONS["sed_shoulders"] + BSA_REGIONS["sed_upperArms"] + BSA_REGIONS["sed_forearms"],  # 15.0
    "lower_legs":   BSA_REGIONS["sed_knees"] + BSA_REGIONS["sed_calves"],   # 8.0
    "full_legs":    BSA_REGIONS["sed_thighs"] + BSA_REGIONS["sed_knees"] + BSA_REGIONS["sed_calves"] + BSA_REGIONS["sed_ankles"],  # 19.0
    "feet":         BSA_REGIONS["sed_feet"],        # 7.0
    "torso_front":  BSA_REGIONS["sed_chest"] + BSA_REGIONS["sed_abdomen"],  # 18.0
    "torso_back":   BSA_REGIONS["sed_back"],        # 18.0
}

# Quick-reference clothing presets (for documentation / sanity checks).
PRESETS = {
    "fully_covered":    _SURVEY_BSA["head"],                                        # ~2.0%
    "dressed_outdoors": _SURVEY_BSA["head"] + _SURVEY_BSA["hands"] + _SURVEY_BSA["forearms"],  # ~13.0%
    "t_shirt_shorts":   _SURVEY_BSA["head"] + _SURVEY_BSA["hands"] + _SURVEY_BSA["full_arms"] + _SURVEY_BSA["lower_legs"],  # ~35.0%
    "swimwear":         sum(_SURVEY_BSA.values()) - _SURVEY_BSA["full_arms"],       # ~80%
}


# =============================================================================
# MODE 1 — SURVEY MODE
# =============================================================================

def compute_bsa_from_survey(
    head:        bool = False,
    hands:       bool = False,
    forearms:    bool = False,
    full_arms:   bool = False,
    lower_legs:  bool = False,
    full_legs:   bool = False,
    feet:        bool = False,
    torso_front: bool = False,
    torso_back:  bool = False,
) -> float:
    """
    Convert a binary exposure survey into a total % BSA value.

    Each parameter indicates whether that body region was directly exposed
    to sunlight on a given day.  The returned value is the ``body_area_pct``
    expected by the original simulate_subject() interface.

    Arm / leg mutual-exclusion guards
    ----------------------------------
    - Set ``forearms=True``  for short sleeves  (elbow to wrist only).
    - Set ``full_arms=True`` for sleeveless      (whole arm exposed).
    - Setting both raises ValueError.
    Same logic applies to ``lower_legs`` vs ``full_legs``.

    Parameters
    ----------
    head        : Face + neck exposed (collar down, no scarf).
    hands       : Both hands exposed (no gloves).
    forearms    : Both forearms exposed (short sleeves).
    full_arms   : Both full arms exposed (sleeveless / tank top).
    lower_legs  : Both lower legs exposed (shorts / skirt below knee).
    full_legs   : Both full legs exposed (very short shorts / swimwear).
    feet        : Both feet exposed (sandals or barefoot).
    torso_front : Front torso exposed (shirt off / open).
    torso_back  : Back torso exposed (backless top / prone sunbathing).

    Returns
    -------
    float
        Total exposed BSA as a percentage (0–100).
    """
    if forearms and full_arms:
        raise ValueError(
            "forearms and full_arms are mutually exclusive — "
            "use full_arms alone if the entire arm is exposed."
        )
    if lower_legs and full_legs:
        raise ValueError(
            "lower_legs and full_legs are mutually exclusive — "
            "use full_legs alone if the entire leg is exposed."
        )

    total = 0.0
    if head:        total += _SURVEY_BSA["head"]
    if hands:       total += _SURVEY_BSA["hands"]
    if forearms:    total += _SURVEY_BSA["forearms"]
    if full_arms:   total += _SURVEY_BSA["full_arms"]
    if lower_legs:  total += _SURVEY_BSA["lower_legs"]
    if full_legs:   total += _SURVEY_BSA["full_legs"]
    if feet:        total += _SURVEY_BSA["feet"]
    if torso_front: total += _SURVEY_BSA["torso_front"]
    if torso_back:  total += _SURVEY_BSA["torso_back"]
    return total


def compute_bsa_for_days(daily_surveys: list) -> list:
    """
    Apply compute_bsa_from_survey() across a list of per-day survey dicts.

    Returns a ``body_area_pct`` list for use with the original
    ``simulate_subject()`` interface.

    Parameters
    ----------
    daily_surveys : list of dict
        One dict per day, with keyword args for compute_bsa_from_survey().
        Omitted keys default to False (region is covered).
        Example::

            [
                {"head": True, "hands": True, "forearms": True},  # day 0
                {},                                                # day 1 — indoors
            ]

    Returns
    -------
    list of float
        BSA percentage for each day.
    """
    return [compute_bsa_from_survey(**s) for s in daily_surveys]


# =============================================================================
# MODE 2 — SENSOR MODE  (per-body-part SED from CSV logs)
# =============================================================================

def effective_dose_from_sed_readings(sed_by_part: dict) -> float:
    """
    Compute the dose-area product from a per-body-part SED measurement dict.

    This replaces the ``E(t) * A(t)`` term in the pharmacokinetic model with
    a more accurate estimate that weights each body site's actual measured
    SED dose by that site's fractional BSA contribution:

        effective_dose = Σ_i  SED_i  ×  (BSA_i / 100)

    The result has units of (SED × fractional BSA) — identical to the
    ``E(t) × A(t)`` product used in Eq. 11 of the source document, so it
    can be passed directly to ``compute_C_sun_effective()`` in
    vitamin_d_model.py.

    Parameters
    ----------
    sed_by_part : dict
        {column_name: daily_total_SED} for each body-part sensor.
        Keys must be a subset of BSA_REGIONS (e.g. "sed_forehead": 0.14).
        Unknown keys are silently ignored.

    Returns
    -------
    float
        Effective dose-area product (SED × fractional BSA).
    """
    return sum(
        sed_by_part.get(col, 0.0) * bsa_pct / 100.0
        for col, bsa_pct in BSA_REGIONS.items()
    )


def effective_doses_from_daily_logs(daily_logs: list) -> list:
    """
    Batch-convert data_loader log summaries to effective dose-area products.

    Convenience wrapper around effective_dose_from_sed_readings() for use
    with the output of data_loader.load_subject_logs().

    Parameters
    ----------
    daily_logs : list of dict
        Output of data_loader.load_subject_logs(); each dict must contain
        a ``sed_by_part`` key.

    Returns
    -------
    list of float
        Effective dose for each day, ready as ``uv_eff_doses`` input to
        vitamin_d_model.run_model_from_effective_doses().
    """
    return [effective_dose_from_sed_readings(d["sed_by_part"]) for d in daily_logs]


# =============================================================================
# ENTRY POINT — example usage
# =============================================================================

if __name__ == "__main__":

    # ------------------------------------------------------------------
    # MODE 1: Survey-based BSA (original approach, no sensor data)
    # ------------------------------------------------------------------
    print("=== Mode 1: Survey-based BSA ===\n")

    surveys_A = [
        {"head": True, "hands": True},                   # day 0
        {"head": True, "hands": True},                   # day 1
        {},                                              # day 2 — indoors
        {"head": True, "hands": True},                   # day 3
        {},                                              # day 4 — indoors
        {"head": True, "hands": True},                   # day 5
    ]
    surveys_B = [
        {"head": True, "hands": True},
        {},
        {"head": True, "hands": True},
        {"head": True, "hands": True},
        {},
        {"head": True, "hands": True},
    ]
    surveys_C = [
        {},
        {"head": True, "hands": True},
        {"head": True, "hands": True},
        {},
        {"head": True, "hands": True, "forearms": True},
        {"head": True, "hands": True},
    ]

    for label, surveys in [("A", surveys_A), ("B", surveys_B), ("C", surveys_C)]:
        bsa_list = compute_bsa_for_days(surveys)
        print(f"# Subject {label}")
        print("body_area_pct = [")
        for i, (survey, bsa) in enumerate(zip(surveys, bsa_list)):
            exposed = [k for k, v in survey.items() if v]
            note    = ", ".join(exposed) if exposed else "fully indoors"
            print(f"    {bsa:.1f},   # day {i}  ({note})")
        print("]\n")

    # ------------------------------------------------------------------
    # MODE 2: Sensor-based effective dose (from CSV log data)
    # ------------------------------------------------------------------
    print("=== Mode 2: Sensor-based effective dose ===\n")

    # Example: one day's worth of per-body-part SED totals (SED units)
    # as would be provided by data_loader.read_uv_log()["sed_by_part"]
    example_sed_day = {
        "sed_scalp":     0.00,
        "sed_forehead":  0.12,
        "sed_ears":      0.08,
        "sed_jaw":       0.10,
        "sed_neck":      0.09,
        "sed_chest":     0.00,   # shirt on
        "sed_abdomen":   0.00,   # shirt on
        "sed_shoulders": 0.05,
        "sed_upperArms": 0.04,
        "sed_forearms":  0.11,
        "sed_hands":     0.13,
        "sed_thighs":    0.00,   # trousers
        "sed_knees":     0.00,
        "sed_calves":    0.00,
        "sed_ankles":    0.00,
        "sed_feet":      0.00,   # shoes
        "sed_back":      0.00,   # shirt on
    }

    eff = effective_dose_from_sed_readings(example_sed_day)
    print(f"Example day effective dose : {eff:.4f}  (SED × fractional BSA)")
    print("(pass this value in uv_eff_doses[] to run_model_from_effective_doses)\n")
