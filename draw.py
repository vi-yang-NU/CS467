"""
draw.py
=======
Stacked bar chart: total vitamin D contribution by source for each participant.

Three segments per bar:
  - Bottom (green) : food intake contribution     (nmol/L)
  - Middle (blue)  : supplement intake contribution (nmol/L)
  - Top    (orange): UV synthesis contribution     (nmol/L)
  - Label above bar: combined total

Run:
    python draw.py
Output:
    vitamin_d_contributions.png  (saved next to this script)
"""

import json
import sys
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))

from vitamin_d_model import (
    compute_C_oral,
    compute_C_sun_effective,
    saturation_factor,
    load_oral_doses_from_json,
    run_model_from_effective_doses,
)
from data_loader import load_subject_logs
from bsa import effective_doses_from_daily_logs

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SCRIPT_DIR       = Path(__file__).parent
DATA_ROOT        = SCRIPT_DIR / "data_resampled"
ORAL_JSON        = SCRIPT_DIR / "vitaminD-mapped.json"
SURFACE_JSON     = SCRIPT_DIR / "surface.json"
USE_SCENARIO_UV  = True   # True = use surface.json; False = use real sensor data

SUBJECTS = [
    ("Oscar (Subject A)",  DATA_ROOT / "uv oscar",        "Oscar",  23, 2, 50.0),
    ("Nicole (Subject B)", DATA_ROOT / "Uv tests Nicole", "Nicole", 20, 2, 50.0),
    ("Vincent (Subject C)",          DATA_ROOT / "vi_logs",         "vi_pdf",    22, 3, 50.0),
]

# Keywords that identify a supplement (vs food)
SUPP_KEYWORDS = {"supplement", "vitamin d", "gummies", "pill", "capsule", "tablet", "d3", "d2"}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_split_oral(json_path, subject_key: str, dates: list):
    """
    Returns (food_ug, supp_ug) — parallel lists in µg/day.
    Looks up each date in the JSON and splits items into food vs supplement.
    """
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    records_by_date = {rec["date"]: rec for rec in data.get(subject_key, [])}

    food_ug = []
    supp_ug = []
    for date_str in dates:
        rec     = records_by_date.get(date_str, {})
        food_iu = 0.0
        supp_iu = 0.0
        for item in rec.get("foods", []):
            name = item.get("name", "").lower()
            vd   = item.get("vitaminD_estimated", 0) or 0
            if any(k in name for k in SUPP_KEYWORDS):
                supp_iu += vd
            else:
                food_iu += vd
        food_ug.append(food_iu / 40.0)
        supp_ug.append(supp_iu / 40.0)

    return food_ug, supp_ug


# ---------------------------------------------------------------------------
# Compute totals for one subject — mirrors run_model_from_effective_doses()
# ---------------------------------------------------------------------------

def compute_totals(log_dir, json_key, age, fitzpatrick, C0):
    logs  = load_subject_logs(log_dir)
    dates = [d["date"] for d in logs]
    n     = len(logs)

    # Patch sed_by_part from surface.json if enabled (mirrors vitamin_d_model.py)
    if USE_SCENARIO_UV and SURFACE_JSON.exists():
        with open(SURFACE_JSON) as fh:
            surface_data = json.load(fh)
        records_by_date = {r["date"]: r["sed_by_part"]
                           for r in surface_data.get(json_key, [])}
        for log in logs:
            if log["date"] in records_by_date:
                log["sed_by_part"] = records_by_date[log["date"]]

    uv_eff_doses = effective_doses_from_daily_logs(logs)

    # Split oral into food vs supplement, then convolve each separately
    # (safe because compute_C_oral is linear — food + supp = total oral)
    food_doses, supp_doses = _load_split_oral(ORAL_JSON, json_key, dates)
    oral_doses = [f + s for f, s in zip(food_doses, supp_doses)]

    C_food = compute_C_oral(food_doses)
    C_supp = compute_C_oral(supp_doses)
    C_oral = compute_C_oral(oral_doses)          # combined — drives C_prev like the model
    C_sun  = compute_C_sun_effective(uv_eff_doses, age, fitzpatrick)

    total_food = 0.0
    total_supp = 0.0
    total_uv   = 0.0
    C_total_prev = C0

    for T in range(n):
        # Exactly the same as run_model_from_effective_doses
        C_oral_prev = 0.0 if T == 0 else C_oral[T - 1]
        C_sun_prev  = 0.0 if T == 0 else C_sun[T - 1]
        C_food_prev = 0.0 if T == 0 else C_food[T - 1]
        C_supp_prev = 0.0 if T == 0 else C_supp[T - 1]

        delta_oral = C_oral[T] - C_oral_prev
        delta_sun  = C_sun[T]  - C_sun_prev
        F          = saturation_factor(C_total_prev)

        total_food += C_food[T] - C_food_prev
        total_supp += C_supp[T] - C_supp_prev
        total_uv   += F * delta_sun

        C_total_prev = C_total_prev + delta_oral + F * delta_sun

    return total_food, total_supp, total_uv


# ---------------------------------------------------------------------------
# Longitudinal C_total series for one subject
# ---------------------------------------------------------------------------

def compute_c_total_series(log_dir, json_key: str, age: int, fitzpatrick: int, C0: float):
    """
    Returns (dates, c_total) — per-day lists.
    dates    : list of str  (YYYY-MM-DD, one per log day)
    c_total  : list of float (nmol/L, cumulative plasma 25(OH)D)
    """
    logs  = load_subject_logs(log_dir)
    dates = [d["date"] for d in logs]

    if USE_SCENARIO_UV and SURFACE_JSON.exists():
        with open(SURFACE_JSON) as fh:
            surface_data = json.load(fh)
        records_by_date = {r["date"]: r["sed_by_part"]
                           for r in surface_data.get(json_key, [])}
        for log in logs:
            if log["date"] in records_by_date:
                log["sed_by_part"] = records_by_date[log["date"]]

    uv_eff_doses = effective_doses_from_daily_logs(logs)
    food_doses, supp_doses = _load_split_oral(ORAL_JSON, json_key, dates)
    oral_doses = [f + s for f, s in zip(food_doses, supp_doses)]

    c_total = run_model_from_effective_doses(oral_doses, uv_eff_doses, age, fitzpatrick, C0)
    return dates, c_total


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def main():
    names      = []
    food_vals  = []
    supp_vals  = []
    uv_vals    = []

    for name, log_dir, json_key, age, fitz, C0 in SUBJECTS:
        if not log_dir.is_dir():
            print(f"[SKIP] {name} — directory not found")
            continue
        food, supp, uv = compute_totals(log_dir, json_key, age, fitz, C0)
        names.append(name)
        food_vals.append(food)
        supp_vals.append(supp)
        uv_vals.append(uv)

    x = list(range(len(names)))

    fig, ax = plt.subplots(figsize=(7, 5))

    b_food = ax.bar(x, food_vals, label="Food",        color="#4CAF50", zorder=3)
    b_supp = ax.bar(x, supp_vals, label="Supplement",  color="#4C9BE8",
                    bottom=food_vals, zorder=3)
    b_uv   = ax.bar(x, uv_vals,   label="UV synthesis", color="#F4A243",
                    bottom=[f + s for f, s in zip(food_vals, supp_vals)], zorder=3)

    # Total label above each bar
    max_total = max(f + s + u for f, s, u in zip(food_vals, supp_vals, uv_vals))
    for i, (f, s, u) in enumerate(zip(food_vals, supp_vals, uv_vals)):
        total = f + s + u
        ax.text(i, total + max_total * 0.02, f"{total:.2f} nmol/L",
                ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Value labels inside each segment (skip if too small)
    threshold = max_total * 0.04
    for i, (f, s, u) in enumerate(zip(food_vals, supp_vals, uv_vals)):
        if f > threshold:
            ax.text(i, f / 2, f"{f:.2f}",
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        if s > threshold:
            ax.text(i, f + s / 2, f"{s:.2f}",
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        if u > threshold:
            ax.text(i, f + s + u / 2, f"{u:.2f}",
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=11)
    ax.set_ylabel("Total vitamin D contribution (nmol/L)", fontsize=10)
    ax.set_title("Total Vitamin D Contribution by Source", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.5, zorder=1)
    ax.set_axisbelow(True)
    ax.set_ylim(0, max_total * 1.18)

    plt.tight_layout()

    out_path = SCRIPT_DIR / "vitamin_d_contributions.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()

    # -----------------------------------------------------------------------
    # Longitudinal C_total line chart
    # -----------------------------------------------------------------------
    COLORS = ["#4CAF50", "#4C9BE8", "#F4A243"]

    fig2, ax2 = plt.subplots(figsize=(8, 5))

    for (name, log_dir, json_key, age, fitz, C0), color in zip(SUBJECTS, COLORS):
        if not log_dir.is_dir():
            continue
        dates, c_total = compute_c_total_series(log_dir, json_key, age, fitz, C0)
        day_labels = [f"Day {i+1}" for i in range(len(dates))]
        ax2.plot(day_labels, c_total, marker="o", label=name, color=color, linewidth=2)
        # annotate final value
        ax2.annotate(
            f"{c_total[-1]:.1f}",
            xy=(day_labels[-1], c_total[-1]),
            xytext=(4, 4), textcoords="offset points",
            fontsize=9, color=color,
        )

    ax2.set_xlabel("Study Day", fontsize=10)
    ax2.set_ylabel("Plasma 25(OH)D — C_total (nmol/L)", fontsize=10)
    ax2.set_title("Longitudinal Vitamin D Estimate", fontsize=13, fontweight="bold")
    ax2.legend(fontsize=10)
    ax2.grid(axis="y", linestyle=":", alpha=0.5)
    ax2.set_axisbelow(True)

    plt.tight_layout()

    out_path2 = SCRIPT_DIR / "vitamin_d_longitudinal.png"
    plt.savefig(out_path2, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path2}")
    plt.show()


if __name__ == "__main__":
    main()
