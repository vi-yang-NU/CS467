"""
bsa.py
======
Body Surface Area (BSA) survey classifier for use with vitamin_d_model.py.

Takes a per-day survey of which body parts were exposed to sunlight and
converts it into a list of daily BSA percentages — the body_area_pct input
expected by simulate_subject().

BSA values are derived from the Lund-Browder chart (standard clinical
reference for body surface area proportions).  Only the portions realistically
exposed to direct sunlight are included — e.g. the scalp is excluded from
"head" because it is almost always covered by hair.

Author  : Claude (Anthropic, claude-sonnet-4-6)
"""


# =============================================================================
# BSA CONTRIBUTIONS PER BODY REGION  (% of total body surface area)
# Source: Lund-Browder chart, adapted for sun-exposure relevance.
# =============================================================================

BSA_REGIONS = {
    # Key            % BSA   Notes
    "head":           4.0,   # Face + neck only; scalp excluded (usually covered by hair)
    "hands":          5.0,   # Both hands, dorsal + palmar surfaces
    "forearms":       6.0,   # Both forearms (elbow to wrist) — short-sleeve exposure
    "full_arms":      9.0,   # Both full arms (shoulder to wrist) — sleeveless exposure
                             #   Note: use forearms OR full_arms, not both
    "lower_legs":    14.0,   # Both lower legs (knee to ankle) — shorts exposure
    "full_legs":     18.0,   # Both full legs (hip to ankle) — swimwear / skirt exposure
                             #   Note: use lower_legs OR full_legs, not both
    "feet":           7.0,   # Both feet (sandals / barefoot)
    "torso_front":   18.0,   # Anterior trunk — t-shirt off / open shirt
    "torso_back":    18.0,   # Posterior trunk — backless clothing / sunbathing
}

# Quick-reference clothing presets (combinations of the above regions)
# These are provided as documentation; compute_bsa_from_survey() builds
# the value from scratch based on the actual survey answers.
PRESETS = {
    "fully_covered":    4.0,   # only face + neck exposed (e.g. coat, gloves, hat off)
    "dressed_outdoors": 15.0,  # head + hands + forearms (typical everyday clothing)
    "t_shirt_shorts":   38.0,  # head + hands + full_arms + lower_legs
    "swimwear":         80.0,  # head + hands + full_arms + full_legs + feet + torso_front + torso_back
}


# =============================================================================
# SURVEY FUNCTION
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

    Each parameter corresponds to a body region being exposed to direct
    sunlight.  Pass True for any region that was uncovered on a given day.

    Arm exposure notes
    ------------------
    - Set forearms=True for short sleeves (elbow to wrist only).
    - Set full_arms=True for sleeveless / tank top (whole arm exposed).
    - Setting both forearms=True and full_arms=True will double-count;
      the function raises an error if both are set.

    Leg exposure notes
    ------------------
    - Set lower_legs=True for shorts/skirts above the knee.
    - Set full_legs=True for very short shorts, swimwear, or no trousers.
    - Same double-count guard applies.

    Parameters
    ----------
    head        : bool — Face + neck exposed (no scarf, collar down).
    hands       : bool — Both hands exposed (no gloves).
    forearms    : bool — Both forearms exposed (short sleeves).
    full_arms   : bool — Both full arms exposed (sleeveless).
    lower_legs  : bool — Both lower legs exposed (shorts / skirt).
    full_legs   : bool — Both full legs exposed (very short / swimwear).
    feet        : bool — Both feet exposed (sandals or barefoot).
    torso_front : bool — Front of torso exposed (no shirt / open shirt).
    torso_back  : bool — Back of torso exposed (backless top / sunbathing).

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
    if head:        total += BSA_REGIONS["head"]
    if hands:       total += BSA_REGIONS["hands"]
    if forearms:    total += BSA_REGIONS["forearms"]
    if full_arms:   total += BSA_REGIONS["full_arms"]
    if lower_legs:  total += BSA_REGIONS["lower_legs"]
    if full_legs:   total += BSA_REGIONS["full_legs"]
    if feet:        total += BSA_REGIONS["feet"]
    if torso_front: total += BSA_REGIONS["torso_front"]
    if torso_back:  total += BSA_REGIONS["torso_back"]

    return total


def compute_bsa_for_days(daily_surveys: list) -> list:
    """
    Apply compute_bsa_from_survey() across a list of per-day survey dicts
    and return the body_area_pct list expected by simulate_subject().

    Each dict in daily_surveys should contain keyword arguments for
    compute_bsa_from_survey() — any omitted keys default to False
    (i.e. that region is covered).

    Parameters
    ----------
    daily_surveys : list
        One dict per day, e.g.:
            [
                {"head": True, "hands": True, "forearms": True},   # day 0
                {"head": True, "hands": True, "full_arms": True},  # day 1
            ]

    Returns
    -------
    list
        BSA percentage for each day, ready to pass as body_area_pct.
    """
    return [compute_bsa_from_survey(**survey) for survey in daily_surveys]


# =============================================================================
# ENTRY POINT — example usage
# =============================================================================

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # One survey dict per day, per subject.
    # Only set True for regions actually exposed to direct sunlight that day.
    # Days with no outdoor exposure should have all keys omitted (or False),
    # which will yield a BSA of 0.0 for that day.
    # -------------------------------------------------------------------------

    # Subject A — 21 y/o, Fitzpatrick II
    # Generally indoors; exposed on days when briefly outside.
    surveys_A = [
        {"head": True, "hands": True},   # day 0 — brief time outdoors
        {"head": True, "hands": True},   # day 1 — brief time outdoors
        {},                              # day 2 — fully indoors
        {"head": True, "hands": True},   # day 3 — brief time outdoors
        {},                              # day 4 — fully indoors
        {"head": True, "hands": True},   # day 5 — brief time outdoors
    ]

    # Subject B — 20 y/o, Fitzpatrick II
    # Similar pattern; one day with slightly more exposure (short sleeves).
    surveys_B = [
        {"head": True, "hands": True},            # day 0
        {},                                        # day 1 — fully indoors
        {"head": True, "hands": True},            # day 2
        {"head": True, "hands": True},            # day 3
        {},                                        # day 4 — fully indoors
        {"head": True, "hands": True},            # day 5
    ]

    # Subject C — 22 y/o, Fitzpatrick III
    # Mostly indoors; one day with short sleeves.
    surveys_C = [
        {},                                               # day 0 — fully indoors
        {"head": True, "hands": True},                   # day 1
        {"head": True, "hands": True},                   # day 2
        {},                                               # day 3 — fully indoors
        {"head": True, "hands": True, "forearms": True}, # day 4 — short sleeves
        {"head": True, "hands": True},                   # day 5
    ]

    # -------------------------------------------------------------------------
    # Compute and print BSA for each subject in simulate_subject() format
    # -------------------------------------------------------------------------
    for label, surveys in [("A", surveys_A), ("B", surveys_B), ("C", surveys_C)]:
        bsa_list = compute_bsa_for_days(surveys)
        print(f"# Subject {label}")
        print("body_area_pct = [")
        for i, (survey, bsa) in enumerate(zip(surveys, bsa_list)):
            exposed = [k for k, v in survey.items() if v]
            note = ", ".join(exposed) if exposed else "fully indoors"
            print(f"    {bsa:.1f},   # day {i}  ({note})")
        print("]")
        print()