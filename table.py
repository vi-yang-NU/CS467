"""
table.py
========
Prints a LaTeX tabular showing, for each study day and each participant:

  SED     — total effective dose-area product from the real sensor logs
             (D_eff = Σ SED_i × BSA_i/100, units: SED × fractional BSA)
  BSA%    — exposed body-surface-area percentage from surface.json
             (sum of BSA_REGIONS[col] for body parts where sed_by_part[col] > 0)
  µg      — total daily oral vitamin D intake (food + supplement) in µg
  nmol/L  — running plasma 25(OH)D estimate C_total from the pharmacokinetic model

Run:
    python table.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from data_loader import load_subject_logs
from bsa import BSA_REGIONS, effective_doses_from_daily_logs
from vitamin_d_model import compute_C_oral, compute_C_sun_effective, saturation_factor

# ---------------------------------------------------------------------------
# Config  (mirrors draw.py / vitamin_d_model.py)
# ---------------------------------------------------------------------------
SCRIPT_DIR      = Path(__file__).parent
DATA_ROOT       = SCRIPT_DIR / "data_resampled"
ORAL_JSON       = SCRIPT_DIR / "vitaminD-mapped.json"
SURFACE_JSON    = SCRIPT_DIR / "surface.json"
USE_SCENARIO_UV = True   # True = use surface.json for model; False = use real sensor

SUPP_KEYWORDS = {"supplement", "vitamin d", "gummies", "pill", "capsule", "tablet", "d3", "d2"}

# (display_name, log_dir, json_key, age, fitzpatrick, C0)
SUBJECTS = [
    ("P1 (Nicole)",  DATA_ROOT / "Uv tests Nicole", "Nicole", 20, 2, 50.0),
    ("P2 (Oscar)",   DATA_ROOT / "uv oscar",        "Oscar",  23, 2, 50.0),
    ("P3 (Vincent)", DATA_ROOT / "vi_logs",         "vi_pdf", 22, 3, 50.0),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_oral_by_date(json_path: Path, subject_key: str, dates: list):
    """Return (food_ug, supp_ug, total_ug) parallel lists."""
    with open(json_path, encoding="utf-8") as fh:
        data = json.load(fh)
    records_by_date = {rec["date"]: rec for rec in data.get(subject_key, [])}

    food_ug, supp_ug = [], []
    for d in dates:
        rec = records_by_date.get(d, {})
        f = s = 0.0
        for item in rec.get("foods", []):
            name = item.get("name", "").lower()
            vd   = item.get("vitaminD_estimated", 0) or 0
            if any(k in name for k in SUPP_KEYWORDS):
                s += vd
            else:
                f += vd
        food_ug.append(f / 40.0)
        supp_ug.append(s / 40.0)

    total_ug = [f + s for f, s in zip(food_ug, supp_ug)]
    return food_ug, supp_ug, total_ug


def _bsa_pct_from_surface(surface_entry: dict) -> float:
    """
    Given one surface.json entry's sed_by_part dict (0/1 values),
    return the total exposed BSA %.
    """
    sed_by_part = surface_entry.get("sed_by_part", {})
    return sum(
        bsa_pct
        for col, bsa_pct in BSA_REGIONS.items()
        if sed_by_part.get(col, 0.0) > 0
    )


def _compute_subject_rows(log_dir, json_key, age, fitzpatrick, C0):
    """
    Returns a list of dicts, one per day:
      date, sed_eff, bsa_pct, oral_ug, c_total
    """
    # --- load real sensor logs (for SED column) ---
    logs  = load_subject_logs(log_dir)
    dates = [d["date"] for d in logs]
    n     = len(logs)

    # Real effective dose (from actual sensor data, before any scenario override)
    real_eff_doses = [d["effective_dose"] for d in logs]

    # --- load surface.json for BSA% column and (optionally) model UV ---
    surface_by_date = {}
    if SURFACE_JSON.exists():
        with open(SURFACE_JSON, encoding="utf-8") as fh:
            surface_data = json.load(fh)
        surface_by_date = {r["date"]: r for r in surface_data.get(json_key, [])}

    # BSA% per day (from surface.json)
    bsa_pcts = []
    for d in dates:
        entry = surface_by_date.get(d)
        bsa_pcts.append(_bsa_pct_from_surface(entry) if entry else 0.0)

    # --- oral doses ---
    _, _, oral_ug = _load_oral_by_date(ORAL_JSON, json_key, dates)

    # --- model UV: scenario or real ---
    if USE_SCENARIO_UV:
        # Override sed_by_part from surface.json (same as vitamin_d_model.py)
        for log in logs:
            if log["date"] in surface_by_date:
                log["sed_by_part"] = surface_by_date[log["date"]]["sed_by_part"]
        model_eff_doses = effective_doses_from_daily_logs(logs)
    else:
        model_eff_doses = real_eff_doses

    # --- run model to get C_total per day ---
    C_oral = compute_C_oral(oral_ug)
    C_sun  = compute_C_sun_effective(model_eff_doses, age, fitzpatrick)

    c_total_series = []
    C_prev = C0
    for T in range(n):
        C_oral_prev = 0.0 if T == 0 else C_oral[T - 1]
        C_sun_prev  = 0.0 if T == 0 else C_sun[T - 1]
        delta_oral  = C_oral[T] - C_oral_prev
        delta_sun   = C_sun[T]  - C_sun_prev
        F           = saturation_factor(C_prev)
        C_prev      = C_prev + delta_oral + F * delta_sun
        c_total_series.append(C_prev)

    rows = []
    for i in range(n):
        rows.append({
            "date":     dates[i],
            "sed_eff":  real_eff_doses[i],   # real sensor D_eff
            "bsa_pct":  bsa_pcts[i],          # surface.json exposed BSA %
            "oral_ug":  oral_ug[i],           # daily oral µg
            "c_total":  c_total_series[i],    # model C_total nmol/L
        })
    return rows


# ---------------------------------------------------------------------------
# Build table data
# ---------------------------------------------------------------------------

all_subject_rows = []
for disp_name, log_dir, json_key, age, fitz, C0 in SUBJECTS:
    if not log_dir.is_dir():
        print(f"[SKIP] {disp_name} — directory not found: {log_dir}", file=sys.stderr)
        all_subject_rows.append([])
        continue
    rows = _compute_subject_rows(log_dir, json_key, age, fitz, C0)
    all_subject_rows.append(rows)
    print(f"Loaded {disp_name}: {len(rows)} days", file=sys.stderr)

# ---------------------------------------------------------------------------
# Print LaTeX table
# ---------------------------------------------------------------------------

n_days = max(len(r) for r in all_subject_rows)
names  = [s[0] for s in SUBJECTS]

print()
print(r"\begin{tabular}{c|cccc|cccc|cccc}")
print(r"\hline")
print(r"& \multicolumn{4}{c|}{\textbf{" + names[0] + r"}} "
      r"& \multicolumn{4}{c|}{\textbf{" + names[1] + r"}} "
      r"& \multicolumn{4}{c}{\textbf{" + names[2] + r"}} \\")
print(
    r"\textbf{Day}"
    r" & \textbf{SED} & \textbf{BSA\%} & \textbf{$\mu$g} & \textbf{nmol/L}"
    r" & \textbf{SED} & \textbf{BSA\%} & \textbf{$\mu$g} & \textbf{nmol/L}"
    r" & \textbf{SED} & \textbf{BSA\%} & \textbf{$\mu$g} & \textbf{nmol/L} \\"
)
print(r"\hline")

for day_idx in range(n_days):
    cells = [str(day_idx + 1)]
    for subj_rows in all_subject_rows:
        if day_idx < len(subj_rows):
            r = subj_rows[day_idx]
            cells += [
                f"{r['sed_eff']:.3f}",
                f"{r['bsa_pct']:.1f}",
                f"{r['oral_ug']:.1f}",
                f"{r['c_total']:.1f}",
            ]
        else:
            cells += ["—", "—", "—", "—"]
    print(" & ".join(cells) + r" \\")

print(r"\hline")
print(r"\end{tabular}")

# ---------------------------------------------------------------------------
# Also print a plain-text preview
# ---------------------------------------------------------------------------
print()
print(f"{'Day':>4}  {'':>46}  {'':>46}  {'':>46}")
header = (
    f"{'Day':>4}  "
    + "  ".join(
        f"{'SED':>6} {'BSA%':>5} {'µg':>6} {'nmol/L':>7}"
        for _ in SUBJECTS
    )
)
print(header)
print("-" * len(header))

for day_idx in range(n_days):
    line = f"{day_idx + 1:>4}  "
    parts = []
    for subj_rows in all_subject_rows:
        if day_idx < len(subj_rows):
            r = subj_rows[day_idx]
            parts.append(
                f"{r['sed_eff']:>6.3f} {r['bsa_pct']:>5.1f} {r['oral_ug']:>6.1f} {r['c_total']:>7.1f}"
            )
        else:
            parts.append(f"{'—':>6} {'—':>5} {'—':>6} {'—':>7}")
    print(line + "  ".join(parts))
