"""
data_loader.py
==============
Ingests the per-body-part UV sensor CSV logs produced by the wearable device
and aggregates them into the per-day summary dictionaries consumed by bsa.py
and vitamin_d_model.py.

Expected CSV columns (from device firmware):
    timestamp, uva, uvb, uvc, uvi,
    sed_raw,
    sed_scalp, sed_forehead, sed_ears, sed_jaw, sed_neck,
    sed_chest, sed_abdomen, sed_shoulders, sed_upperArms,
    sed_forearms, sed_hands, sed_thighs, sed_knees, sed_calves,
    sed_ankles, sed_feet, sed_back,
    i_o, notification, response, offline

One CSV file covers one calendar day.  Each row is a timestamped sensor
reading.  The sed_* columns contain incremental SED values (dose accumulated
since the previous reading, in SED units); summing all rows for a given body
part yields the total daily SED dose at that location.

Author  : Claude (Anthropic, claude-sonnet-4-6)
"""

import csv
import os
from pathlib import Path
from datetime import datetime


# =============================================================================
# COLUMN DEFINITIONS
# =============================================================================

# All per-body-part SED column names as they appear in the CSV header.
SED_COLUMNS = [
    "sed_scalp",
    "sed_forehead",
    "sed_ears",
    "sed_jaw",
    "sed_neck",
    "sed_chest",
    "sed_abdomen",
    "sed_shoulders",
    "sed_upperArms",
    "sed_forearms",
    "sed_hands",
    "sed_thighs",
    "sed_knees",
    "sed_calves",
    "sed_ankles",
    "sed_feet",
    "sed_back",
]

# Columns that carry scalar UV environment measurements (not body-part specific).
UV_ENV_COLUMNS = ["uva", "uvb", "uvc", "uvi", "sed_raw"]


# =============================================================================
# BSA WEIGHTS FOR EACH SENSOR LOCATION
# (% of total body surface area — Lund-Browder chart)
# These map each sensor site to its fractional contribution to the
# dose-weighted effective input E(t)*A(t) used in the pharmacokinetic model.
# =============================================================================

BSA_WEIGHT = {
    # Sensor location   % BSA
    "sed_scalp":         1.5,   # Top of head only (hair-covered area)
    "sed_forehead":      1.5,   # Upper face
    "sed_ears":          0.5,   # Auricles
    "sed_jaw":           0.5,   # Lower face / chin
    "sed_neck":          1.0,   # Anterior + posterior neck
    "sed_chest":         9.0,   # Upper anterior trunk
    "sed_abdomen":       9.0,   # Lower anterior trunk
    "sed_shoulders":     3.0,   # Both shoulder caps
    "sed_upperArms":     6.0,   # Both upper arms
    "sed_forearms":      6.0,   # Both forearms
    "sed_hands":         5.0,   # Both hands (dorsal + palmar)
    "sed_thighs":        9.5,   # Both anterior thighs
    "sed_knees":         1.0,   # Both knee caps
    "sed_calves":        7.0,   # Both calves
    "sed_ankles":        1.5,   # Both ankle regions
    "sed_feet":          7.0,   # Both feet
    "sed_back":         18.0,   # Entire posterior trunk
}

# Sanity-check: BSA weights should add to ~100 if the whole body is covered.
_TOTAL_BSA = sum(BSA_WEIGHT.values())   # ≈ 86.5% (remaining ~13.5% = genitals/gluteal)


# =============================================================================
# SINGLE-FILE READER
# =============================================================================

def _safe_float(value: str) -> float:
    """Convert a CSV cell to float; return 0.0 for empty / non-numeric cells."""
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0


def read_uv_log(filepath: str | Path) -> dict:
    """
    Read a single daily UV-log CSV and return a summary dict for that day.

    The summary dict contains:

    ``uvb_j_m2``
        Total daily UVB dose in J/m² (sum of the ``uvb`` column).
        Used as fallback input for the original uvb_dose_to_sed() path.

    ``sed_raw``
        Total raw SED accumulated over the day (sum of ``sed_raw`` column).

    ``sed_by_part``
        {column_name: daily_total_SED} for every sed_* column.

    ``effective_dose``
        Σ_i (SED_i × BSA_weight_i / 100)  —  the dose-area product
        (units: SED × fractional body area) ready to substitute for the
        ``E(t) * A(t)`` term in compute_C_sun().  This is the primary
        output used by the updated model.

    ``fraction_outdoor``
        Fraction of readings where i_o == 1 (outdoor), 0.0–1.0.

    ``n_rows``
        Number of data rows parsed (useful for QA).

    Parameters
    ----------
    filepath : str or Path
        Path to the daily CSV log file.

    Returns
    -------
    dict
        Summary as described above.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"UV log not found: {filepath}")

    sed_totals   = {col: 0.0 for col in SED_COLUMNS}
    uvb_total    = 0.0
    sed_raw_total = 0.0
    outdoor_count = 0
    n_rows       = 0

    import io
    with open(filepath, newline="", encoding="utf-8") as fh:
        content = fh.read()   # read entire file at once — avoids network timeout on row-by-row reads
    reader = csv.DictReader(io.StringIO(content))
    for row in reader:
            n_rows += 1

            # Environmental UV
            uvb_total     += _safe_float(row.get("uvb",     "0"))
            sed_raw_total += _safe_float(row.get("sed_raw", "0"))

            # Per-body-part SED (incremental)
            for col in SED_COLUMNS:
                sed_totals[col] += _safe_float(row.get(col, "0"))

            # Indoor / outdoor flag  (1 = outdoor, 0 = indoor)
            io_val = row.get("i_o", "0").strip()
            if io_val in ("1", "1.0", "True", "true", "outdoor"):
                outdoor_count += 1

    # Dose-area weighted effective input: Σ_i (SED_i * BSA_fraction_i)
    # BSA_fraction_i = BSA_weight[i] / 100  (convert % to decimal)
    effective_dose = sum(
        sed_totals[col] * BSA_WEIGHT[col] / 100.0
        for col in SED_COLUMNS
    )

    return {
        "filepath":         str(filepath),
        "uvb_j_m2":         uvb_total,
        "sed_raw":          sed_raw_total,
        "sed_by_part":      dict(sed_totals),
        "effective_dose":   effective_dose,
        "fraction_outdoor": outdoor_count / n_rows if n_rows > 0 else 0.0,
        "n_rows":           n_rows,
    }


# =============================================================================
# MULTI-FILE / SUBJECT LOADER
# =============================================================================

def load_subject_logs(
    log_dir: str | Path,
    date_range: list[str] | None = None,
) -> list[dict]:
    """
    Load all daily UV-log CSVs for one subject from a directory.

    Files are expected to follow the naming convention:
        YYYY-MM-DD_UV_LOG.csv

    They are loaded in chronological order.  If ``date_range`` is given,
    only files whose date falls within the list are included.

    Parameters
    ----------
    log_dir    : str or Path
        Directory containing the subject's CSV files
        (e.g. ``data/uv oscar`` or ``data/Uv tests Nicole``).
    date_range : list of str, optional
        Ordered list of date strings  ``["YYYY-MM-DD", ...]``.
        If None, all CSV files in the directory are loaded in sorted order.

    Returns
    -------
    list of dict
        One summary dict per day (from read_uv_log()), sorted chronologically.
        Pass the ``effective_dose`` values from each dict as ``uv_eff_doses``
        to ``simulate_subject_from_logs()``.
    """
    log_dir = Path(log_dir)
    if not log_dir.is_dir():
        raise NotADirectoryError(f"Log directory not found: {log_dir}")

    # Collect candidate files
    csv_files = sorted(log_dir.glob("*_UV_LOG.csv"))

    if date_range is not None:
        date_set   = set(date_range)
        csv_files  = [
            f for f in csv_files
            if f.name[:10] in date_set   # file names start with YYYY-MM-DD
        ]

    if not csv_files:
        raise FileNotFoundError(
            f"No matching UV_LOG CSV files found in {log_dir}"
        )

    return [read_uv_log(f) for f in csv_files]


def extract_effective_doses(daily_logs: list[dict]) -> list[float]:
    """
    Pull the ``effective_dose`` scalar from each day's log summary.

    This is the primary adapter between the data-loader output and the
    pharmacokinetic model.  Pass the result directly as the ``uv_eff_doses``
    argument to ``run_model_from_effective_doses()`` in vitamin_d_model.py.

    Parameters
    ----------
    daily_logs : list of dict
        Output of load_subject_logs().

    Returns
    -------
    list of float
        Effective dose-area product for each day  (SED × fractional BSA).
    """
    return [d["effective_dose"] for d in daily_logs]


def extract_uvb_doses(daily_logs: list[dict]) -> list[float]:
    """
    Pull the raw UVB column totals (J/m²) for use with the legacy
    uvb_dose_to_sed() code path in vitamin_d_model.py.

    Parameters
    ----------
    daily_logs : list of dict
        Output of load_subject_logs().

    Returns
    -------
    list of float
        Total UVB dose for each day in J/m².
    """
    return [d["uvb_j_m2"] for d in daily_logs]


# =============================================================================
# DIAGNOSTIC PRINTING
# =============================================================================

def print_log_summary(daily_logs: list[dict], subject_name: str = "", file=None) -> None:
    """Print a human-readable summary table of loaded log data."""
    label = f"  Subject: {subject_name}" if subject_name else "  Log Summary"
    print("=" * 70, file=file)
    print(label, file=file)
    print(f"  {'Day':>4}  {'File':>30}  {'Rows':>5}  {'Eff.Dose':>10}  {'%Outdoor':>9}", file=file)
    print("  " + "-" * 63, file=file)
    for i, d in enumerate(daily_logs):
        fname    = Path(d["filepath"]).name[:30]
        rows     = d["n_rows"]
        eff      = d["effective_dose"]
        pct_out  = d["fraction_outdoor"] * 100
        print(f"  {i:>4}  {fname:>30}  {rows:>5}  {eff:>10.4f}  {pct_out:>8.1f}%", file=file)
    print(file=file)


# =============================================================================
# ENTRY POINT — example usage
# =============================================================================

if __name__ == "__main__":
    import sys

    DATA_ROOT = Path(__file__).parent / "data"

    # Map subject labels to their log subdirectories
    SUBJECTS = {
        "Oscar":  DATA_ROOT / "uv oscar",
        "Nicole": DATA_ROOT / "Uv tests Nicole",
        "vi_logs (shared)": DATA_ROOT / "vi_logs",
    }

    DATE_RANGE = [
        "2026-03-04",
        "2026-03-05",
        "2026-03-06",
        "2026-03-07",
        "2026-03-08",
        "2026-03-09",
    ]

    for name, directory in SUBJECTS.items():
        if not directory.is_dir():
            print(f"  [SKIP] Directory not found: {directory}")
            continue
        try:
            logs = load_subject_logs(directory, date_range=DATE_RANGE)
            print_log_summary(logs, subject_name=name)

            # Ready-to-use values for the pharmacokinetic model
            eff_doses = extract_effective_doses(logs)
            uvb_doses = extract_uvb_doses(logs)

            print(f"  effective_doses = {[round(v, 4) for v in eff_doses]}")
            print(f"  uvb_j_m2        = {[round(v, 1) for v in uvb_doses]}")
            print()

        except (FileNotFoundError, NotADirectoryError) as exc:
            print(f"  [ERROR] {exc}")
