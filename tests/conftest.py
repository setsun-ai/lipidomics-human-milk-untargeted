"""Synthetic input data in the same format as the restricted course dataset.

The real feature table and metadata cannot be published, so the tests build a small
fake study: 5 paired donors (M1 + M12) plus one donor with M12 samples only,
2 biological x 2 extraction replicates, 10 pooled QCs and 1 extraction blank.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from openpyxl import Workbook

PAIRED = ["AB", "CD", "EFG", "HI", "JK"]  # made-up donor codes
UNPAIRED = "XY"
N_FEATURES = 60


def _sample_columns():
    cols, meta, suffix = [], [], 1000
    for donor in PAIRED + [UNPAIRED]:
        months = [12] if donor == UNPAIRED else [1, 12]
        for m in months:
            for b in (1, 2):
                for e in ("A", "B"):
                    suffix += 1
                    cols.append(f"[Vol] {donor} {m}_{b}{e} {suffix}")
                    meta.append({
                        "Sample_ID": f"{donor}_M{m}_B{b}_E{e}", "File_name": f"{donor}_{m}_{b}{e}_{suffix}.mzML",
                        "Sample_type": "Sample", "Donor_ID": donor, "Timepoint": f"M{m}", "Lactation_month": m,
                        "BioRep_ID": b, "ExtractionRep_ID": e, "Analytical_batch": "Batch_1",
                        "Measurement_day": "day6", "Notes": None,
                    })
    return cols, meta


def _write_sheet(wb, name, frame):
    ws = wb.create_sheet(name)
    ws.append([f"{name} (synthetic)"])  # descriptive title row, as in the real workbook
    ws.append(list(frame.columns))
    for row in frame.itertuples(index=False):
        ws.append([None if (isinstance(v, float) and np.isnan(v)) else v for v in row])


@pytest.fixture(scope="session")
def synthetic_input(tmp_path_factory) -> Path:
    rng = np.random.default_rng(0)
    d = tmp_path_factory.mktemp("input")
    sample_cols, sample_meta = _sample_columns()
    qc_cols = [f"[Vol] QC_day {day}_{k}" for day in (6, 7, 8, 9, 16) for k in (1, 2)]
    blank_col = "[Vol] Blank extraction _day 9"

    base = rng.lognormal(mean=13, sigma=1.0, size=N_FEATURES)
    table = {"Mass (avg)": rng.uniform(300, 900, N_FEATURES).round(4),
             "RT (avg)": rng.uniform(1, 35, N_FEATURES).round(2),
             "Height (avg)": (base * 0.5).round(0)}
    for col in sample_cols:
        shift = 1.5 if " 12_" in col else 1.0
        vals = base * shift * rng.lognormal(0, 0.25, N_FEATURES)
        vals[rng.random(N_FEATURES) < 0.05] = np.nan  # a few missing values
        table[col] = vals
    for col in qc_cols:
        table[col] = base * rng.lognormal(0, 0.08, N_FEATURES)
    table[blank_col] = base * 0.01
    pd.DataFrame(table).to_csv(d / "nontargeted table 5 000 threshold.csv", index=False)

    qc_meta = [{"Sample_ID": f"QC_{i:02d}", "File_name": c.replace("[Vol] ", "").replace(" ", "_") + ".mzML",
                "Sample_type": "QC", "Donor_ID": "POOL", "Timepoint": "QC"} for i, c in enumerate(qc_cols, 1)]
    blank_meta = [{"Sample_ID": "BL_01", "File_name": "Blank_extraction_day_9.mzML", "Sample_type": "Blank",
                   "Donor_ID": "NA", "Timepoint": "Blank"}]
    # the real workbook lists only donors with both timepoints
    paired_meta = [m for m in sample_meta if m["Donor_ID"] != UNPAIRED]
    sample_df = pd.DataFrame(paired_meta + qc_meta + blank_meta)
    seq = sample_df[["File_name", "Sample_ID", "Sample_type"]].copy()
    seq.insert(0, "Injection_order", range(1, len(seq) + 1))

    wb = Workbook()
    wb.remove(wb.active)
    _write_sheet(wb, "Sample_metadata", sample_df)
    _write_sheet(wb, "Sequence", seq)
    # demographic sheet present in the real workbook; the analysis must not need it
    _write_sheet(wb, "Donor_metadata", pd.DataFrame({"Donor_ID": PAIRED, "Maternal_age_at_delivery": [30] * 5}))
    wb.save(d / "metadata_M1_M12 2 05 2026.xlsx")
    return d
