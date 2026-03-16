"""
resample_5s.py
==============
Resamples all UV sensor CSV logs to a uniform 5-second window.

Problem
-------
The wearable device logs at inconsistent intervals:
  - Nicole  : ~1-second cadence
  - Oscar   : mixed 1-second / 5-second, with out-of-order timestamps
  - vi_logs : ~5-second cadence

This script normalises all three to 5-second windows so downstream
analysis (data_loader.py, vitamin_d_model.py) sees consistent data.

Resampling rules
----------------
  1. Sort rows by timestamp (fixes Oscar's out-of-order data).
  2. Bin rows into 5-second windows:
         bin_start = (timestamp // 5) * 5
  3. For each bin:
       - Numeric columns   → mean of all rows in the bin.
       - Non-numeric cols  → first non-empty value in the bin.
       - Output timestamp  → bin_start (integer).
  4. Windows with only one reading are passed through unchanged
     (already at 5-second resolution).

Output
------
  data_resampled/
    uv oscar/              ← resampled copies of data/uv oscar/
    Uv tests Nicole/       ← resampled copies of data/Uv tests Nicole/
    vi_logs/               ← resampled copies of data/vi_logs/

Each output file keeps the same filename as the original (YYYY-MM-DD_UV_LOG.csv).

Usage
-----
  python resample_5s.py

A summary is printed for every file processed.
"""

import csv
import io
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

SCRIPT_DIR   = Path(__file__).parent
DATA_ROOT    = SCRIPT_DIR / "data"
OUTPUT_ROOT  = SCRIPT_DIR / "data_resampled"

SUBJECTS = [
    "uv oscar",
    "Uv tests Nicole",
    "vi_logs",
]

WINDOW_SEC = 5   # target window size in seconds

# Instantaneous sensor readings — average across rows in the window.
AVG_COLS = {
    "uva", "uvb", "uvc", "uvi",
}

# Incremental dose values — sum across rows in the window to preserve total dose.
SUM_COLS = {
    "sed_raw",
    "sed_scalp", "sed_forehead", "sed_ears", "sed_jaw", "sed_neck",
    "sed_chest", "sed_abdomen", "sed_shoulders", "sed_upperArms",
    "sed_forearms", "sed_hands", "sed_thighs", "sed_knees", "sed_calves",
    "sed_ankles", "sed_feet", "sed_back",
}

# Integer flags — round the average (majority vote).
FLAG_COLS = {"i_o", "offline"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _try_float(value: str):
    """Return (float, True) if parseable, else ('', False)."""
    try:
        return float(value), True
    except (ValueError, TypeError):
        return None, False


def resample_file(src: Path, dst: Path) -> dict:
    """
    Resample a single CSV to 5-second windows.

    Returns a stats dict with keys:
        rows_in, rows_out, files_with_negatives (bool), subject
    """
    # --- read ---
    with open(src, newline="", encoding="utf-8") as fh:
        content = fh.read()

    reader   = csv.DictReader(io.StringIO(content))
    fieldnames = reader.fieldnames or []
    rows     = list(reader)

    if not rows:
        dst.parent.mkdir(parents=True, exist_ok=True)
        with open(dst, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
        return {"rows_in": 0, "rows_out": 0, "had_negatives": False}

    # --- sort by timestamp (fixes out-of-order data) ---
    def _ts(row):
        v, ok = _try_float(row.get("timestamp", "0"))
        return v if ok else 0.0

    rows.sort(key=_ts)

    # Check whether the original file had negative gaps (for reporting only)
    ts_vals = [_ts(r) for r in rows]
    had_negatives = any(
        ts_vals[i + 1] < ts_vals[i] for i in range(len(ts_vals) - 1)
    )

    # --- bin into 5-second windows ---
    # bins[bin_start] = list of rows
    bins: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        v, ok = _try_float(row.get("timestamp", "0"))
        bin_start = int(v // WINDOW_SEC) * WINDOW_SEC if ok else 0
        bins[bin_start].append(row)

    # --- aggregate each bin ---
    out_rows = []
    for bin_start in sorted(bins.keys()):
        group = bins[bin_start]

        agg: dict[str, str] = {}

        for col in fieldnames:
            if col == "timestamp":
                agg[col] = str(bin_start)
                continue

            if col in AVG_COLS or col in SUM_COLS or col in FLAG_COLS:
                nums = []
                for r in group:
                    v, ok = _try_float(r.get(col, ""))
                    if ok and v is not None:
                        nums.append(v)
                if nums:
                    if col in SUM_COLS:
                        result = sum(nums)
                    elif col in FLAG_COLS:
                        result = int(round(sum(nums) / len(nums)))
                    else:  # AVG_COLS
                        result = sum(nums) / len(nums)
                    agg[col] = str(result)
                else:
                    agg[col] = ""
            else:
                # categorical: first non-empty value
                val = ""
                for r in group:
                    candidate = (r.get(col) or "").strip()
                    if candidate:
                        val = candidate
                        break
                agg[col] = val

        out_rows.append(agg)

    # --- write ---
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    return {
        "rows_in":       len(rows),
        "rows_out":      len(out_rows),
        "had_negatives": had_negatives,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"Resampling UV logs → {OUTPUT_ROOT}\n")
    print(f"{'Subject':<22}  {'File':<30}  {'Rows in':>8}  {'Rows out':>9}  {'Sorted?':>7}")
    print("-" * 85)

    total_in = total_out = 0

    for subject in SUBJECTS:
        src_dir = DATA_ROOT / subject
        dst_dir = OUTPUT_ROOT / subject

        if not src_dir.is_dir():
            print(f"  [SKIP] {src_dir} — directory not found")
            continue

        csv_files = sorted(src_dir.glob("*_UV_LOG.csv"))
        if not csv_files:
            print(f"  [SKIP] {src_dir} — no *_UV_LOG.csv files")
            continue

        for i, src in enumerate(csv_files, start=1):
            dst  = dst_dir / f"day{i}_UV_LOG.csv"
            info = resample_file(src, dst)

            sorted_flag = "yes" if info["had_negatives"] else "—"
            print(
                f"  {subject:<20}  {src.name:<30}  "
                f"{info['rows_in']:>8,}  {info['rows_out']:>9,}  {sorted_flag:>7}"
            )
            total_in  += info["rows_in"]
            total_out += info["rows_out"]

    print("-" * 85)
    print(f"  {'TOTAL':<50}  {total_in:>8,}  {total_out:>9,}")
    print(f"\nDone.  Resampled files written to: {OUTPUT_ROOT}/")


if __name__ == "__main__":
    main()
