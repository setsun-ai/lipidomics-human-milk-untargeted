#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Non-targeted lipidomics analysis for the project dataset.

Targeted analysis is intentionally omitted. This script uses only:
- nontargeted table 5 000 threshold.csv
- metadata_M1_M12 2 05 2026.xlsx

Outputs are written to figures/ and processed_data/.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import warnings
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
from zipfile import ZipFile

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=RuntimeWarning)


XLSX_NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
}


@dataclass
class ProjectPaths:
    input_dir: Path
    output_dir: Path
    figures_dir: Path
    processed_dir: Path


def column_letters_to_index(col: str) -> int:
    idx = 0
    for char in col:
        idx = idx * 26 + ord(char.upper()) - 64
    return idx - 1


def cell_reference_to_position(ref: str) -> Tuple[int, int]:
    match = re.match(r"([A-Z]+)(\d+)", ref)
    if not match:
        raise ValueError(f"Unsupported cell reference: {ref}")
    col, row = match.groups()
    return int(row) - 1, column_letters_to_index(col)


def read_xlsx_sheets(xlsx_path: Path) -> Dict[str, pd.DataFrame]:
    """Read a simple .xlsx workbook without using openpyxl.

    The metadata workbook contains simple rectangular sheets. This parser supports
    shared strings and numeric/string cells, which is enough for this dataset.
    """
    sheets: Dict[str, pd.DataFrame] = {}
    with ZipFile(xlsx_path) as archive:
        shared_strings: List[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for si in root.findall("main:si", XLSX_NS):
                text = "".join(t.text or "" for t in si.iter(f"{{{XLSX_NS['main']}}}t"))
                shared_strings.append(text)

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall("rel:Relationship", XLSX_NS)
        }

        for sheet in workbook.findall("main:sheets/main:sheet", XLSX_NS):
            sheet_name = sheet.attrib["name"]
            rid = sheet.attrib[f"{{{XLSX_NS['r']}}}id"]
            target = rid_to_target[rid].lstrip("/")  # targets may be absolute ("/xl/worksheets/...")
            xml_path = target if target.startswith("xl/") else f"xl/{target}"
            root = ET.fromstring(archive.read(xml_path))

            cells = []
            max_row = 0
            max_col = 0
            for cell in root.findall(".//main:c", XLSX_NS):
                ref = cell.attrib.get("r")
                if not ref:
                    continue
                row_idx, col_idx = cell_reference_to_position(ref)
                max_row = max(max_row, row_idx)
                max_col = max(max_col, col_idx)
                cell_type = cell.attrib.get("t")
                value_node = cell.find("main:v", XLSX_NS)
                value = None

                if cell_type == "s" and value_node is not None:
                    value = shared_strings[int(value_node.text)]
                elif cell_type == "inlineStr":
                    value = "".join(t.text or "" for t in cell.iter(f"{{{XLSX_NS['main']}}}t"))
                elif value_node is not None:
                    raw_value = value_node.text
                    if cell_type == "b":
                        value = bool(int(raw_value))
                    else:
                        try:
                            as_float = float(raw_value)
                            value = int(as_float) if as_float.is_integer() else as_float
                        except (TypeError, ValueError):
                            value = raw_value
                cells.append((row_idx, col_idx, value))

            grid = [[None] * (max_col + 1) for _ in range(max_row + 1)]
            for row_idx, col_idx, value in cells:
                grid[row_idx][col_idx] = value

            header_row_index = None
            for idx, row in enumerate(grid):
                # First real header after the descriptive title row.
                if idx > 0 and sum(value is not None for value in row) >= 2:
                    header_row_index = idx
                    break
            if header_row_index is None:
                sheets[sheet_name] = pd.DataFrame(grid)
                continue

            header = [str(value) if value is not None else "" for value in grid[header_row_index]]
            rows = [row[: len(header)] for row in grid[header_row_index + 1 :] if any(v is not None for v in row)]
            sheets[sheet_name] = pd.DataFrame(rows, columns=header)
    return sheets


def parse_volume_label(column: str) -> Dict[str, object]:
    """Infer sample metadata from a [Vol] column header."""
    label = column.split("] ", 1)[1].strip()
    compact = " ".join(label.split())

    if compact.lower().startswith("blank"):
        return {
            "column": column,
            "label": compact,
            "file_name_inferred": "Blank_extraction_day_9.mzML",
            "sample_type": "Blank",
            "donor": "NA",
            "timepoint": "Blank",
            "sample_id_inferred": "BL_01",
            "bio_rep": None,
            "extraction_rep": None,
            "measurement_day": "day9",
            "numeric_suffix": None,
        }

    if compact.startswith("QC_day"):
        match = re.match(r"QC_day\s*([0-9]+)_([0-9]+)", compact)
        day, suffix = match.groups() if match else (None, None)
        return {
            "column": column,
            "label": compact,
            "file_name_inferred": f"QC_day_{day}_{suffix}.mzML" if day else compact.replace(" ", "_") + ".mzML",
            "sample_type": "QC",
            "donor": "POOL",
            "timepoint": "QC",
            "sample_id_inferred": f"QC_day_{day}_{suffix}" if day else compact,
            "bio_rep": None,
            "extraction_rep": None,
            "measurement_day": f"day{day}" if day else None,
            "numeric_suffix": suffix,
        }

    match = re.match(r"([A-Za-z][A-Za-z0-9]*)\s+([0-9]+)_([0-9]+)([A-Za-z])\s+([0-9]+)", compact)
    if match:
        donor, month, bio_rep, extraction_rep, suffix = match.groups()
        return {
            "column": column,
            "label": compact,
            "file_name_inferred": f"{donor}_{month}_{bio_rep}{extraction_rep}_{suffix}.mzML",
            "sample_type": "Sample",
            "donor": donor,
            "timepoint": f"M{month}",
            "sample_id_inferred": f"{donor}_M{month}_B{bio_rep}_E{extraction_rep}",
            "bio_rep": int(bio_rep),
            "extraction_rep": extraction_rep,
            "measurement_day": None,
            "numeric_suffix": suffix,
        }

    return {
        "column": column,
        "label": compact,
        "file_name_inferred": compact.replace(" ", "_") + ".mzML",
        "sample_type": "Unknown",
        "donor": None,
        "timepoint": None,
        "sample_id_inferred": None,
        "bio_rep": None,
        "extraction_rep": None,
        "measurement_day": None,
        "numeric_suffix": None,
    }


def pseudonymise_donors(
    feature_table: pd.DataFrame, sample_metadata: pd.DataFrame, sequence: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Replace donor codes with neutral aliases (D1, D2, ...) before any processing.

    The source files use short, initials-like donor codes. Aliases are assigned at run time in
    sorted order, so the original codes never appear in outputs, figures or this source file.
    """
    codes = {
        parse_volume_label(col)["donor"]
        for col in feature_table.columns
        if col.startswith("[Vol] ") and parse_volume_label(col)["sample_type"] == "Sample"
    }
    codes |= set(sample_metadata["Donor_ID"].dropna().astype(str))
    codes -= {"NA", "POOL", "None", ""}
    alias = {code: f"D{i}" for i, code in enumerate(sorted(codes), start=1)}
    # longest codes first so that e.g. a three-letter code is not partially matched by a two-letter one
    pattern = re.compile(
        r"(?<![A-Za-z])(" + "|".join(sorted(map(re.escape, alias), key=len, reverse=True)) + r")(?=[_\s]|$)"
    )

    def sub(value):
        return pattern.sub(lambda m: alias[m.group(1)], value) if isinstance(value, str) else value

    feature_table = feature_table.rename(columns=sub)
    def sub_frame(frame: pd.DataFrame) -> pd.DataFrame:
        return frame.apply(lambda col: col.map(sub) if not pd.api.types.is_numeric_dtype(col) else col)

    sample_metadata, sequence = sub_frame(sample_metadata), sub_frame(sequence)
    return feature_table, sample_metadata, sequence


def to_numeric_frame(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(pd.to_numeric, errors="coerce")


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return q
    p_valid = p[valid]
    order = np.argsort(p_valid)
    ranked = p_valid[order]
    m = len(ranked)
    adjusted = ranked * m / (np.arange(m) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    q[valid] = restored
    return q


def impute_half_min_positive(matrix: pd.DataFrame) -> pd.DataFrame:
    """Replace missing/non-positive values by half of the minimum positive value per feature."""
    result = matrix.copy().astype(float)
    for idx in result.index:
        row = result.loc[idx]
        positives = row[row > 0]
        if len(positives) == 0:
            fill_value = 1.0
        else:
            fill_value = positives.min() / 2.0
        result.loc[idx] = row.where(row > 0, fill_value).fillna(fill_value)
    return result


def safe_wilcoxon(x: np.ndarray, y: np.ndarray) -> float:
    try:
        if np.allclose(x, y, equal_nan=False):
            return 1.0
        return float(wilcoxon(x, y, zero_method="wilcox", alternative="two-sided", mode="auto").pvalue)
    except Exception:
        return np.nan


def safe_mann_whitney(x: np.ndarray, y: np.ndarray) -> float:
    try:
        return float(mannwhitneyu(x, y, alternative="two-sided", method="auto").pvalue)
    except Exception:
        return np.nan


def prepare_paths(input_dir: Path | None, output_dir: Path | None) -> ProjectPaths:
    script_dir = Path(__file__).resolve().parent
    input_dir = input_dir or script_dir / "input_data"
    output_dir = output_dir or script_dir
    figures_dir = output_dir / "figures"
    processed_dir = output_dir / "processed_data"
    figures_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    return ProjectPaths(input_dir=input_dir, output_dir=output_dir, figures_dir=figures_dir, processed_dir=processed_dir)


def run_analysis(paths: ProjectPaths) -> Dict[str, object]:
    feature_path = paths.input_dir / "nontargeted table 5 000 threshold.csv"
    metadata_path = paths.input_dir / "metadata_M1_M12 2 05 2026.xlsx"

    if not feature_path.exists():
        raise FileNotFoundError(f"Missing feature table: {feature_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing metadata workbook: {metadata_path}")

    feature_table = pd.read_csv(feature_path)
    metadata_sheets = read_xlsx_sheets(metadata_path)
    sample_metadata = metadata_sheets["Sample_metadata"].copy()
    sequence = metadata_sheets["Sequence"].copy()
    # Donor_metadata and Timepoint_metadata (age, BMI, delivery mode, diet, health status, ...)
    # are deliberately not loaded: the analysis does not need them (data minimisation).

    # Drop entirely empty columns produced by the CSV export.
    feature_table = feature_table.loc[:, ~feature_table.columns.str.match(r"^Unnamed")]
    feature_table, sample_metadata, sequence = pseudonymise_donors(feature_table, sample_metadata, sequence)
    feature_table.insert(0, "feature_id", [f"F{i:04d}" for i in range(1, len(feature_table) + 1)])

    vol_columns = [col for col in feature_table.columns if col.startswith("[Vol] ")]
    vol_meta = pd.DataFrame([parse_volume_label(col) for col in vol_columns])

    # Attach sample metadata where possible. Biological samples are linked by inferred Sample_ID,
    # because the workbook contains several exact file-name inconsistencies for one donor's samples.
    sample_lookup_by_id = sample_metadata.set_index("Sample_ID", drop=False)
    sample_lookup_by_file = sample_metadata.set_index("File_name", drop=False)

    rows = []
    for _, row in vol_meta.iterrows():
        out = row.to_dict()
        if row["sample_type"] == "Sample" and row["sample_id_inferred"] in sample_lookup_by_id.index:
            md = sample_lookup_by_id.loc[row["sample_id_inferred"]].to_dict()
            out.update({f"metadata_{k}": v for k, v in md.items()})
            out["matched_to_metadata"] = True
            out["match_mode"] = "Sample_ID"
        elif row["file_name_inferred"] in sample_lookup_by_file.index:
            md = sample_lookup_by_file.loc[row["file_name_inferred"]].to_dict()
            out.update({f"metadata_{k}": v for k, v in md.items()})
            out["matched_to_metadata"] = True
            out["match_mode"] = "File_name"
        else:
            out["matched_to_metadata"] = False
            out["match_mode"] = "unmatched"
        rows.append(out)
    vol_meta = pd.DataFrame(rows)

    # Paired design: keep only donors with samples at both timepoints (one donor has no M1 samples).
    samples = vol_meta.query("sample_type == 'Sample'")
    paired_donors = [d for d, tps in samples.groupby("donor")["timepoint"] if {"M1", "M12"} <= set(tps)]
    used_sample_meta = samples[samples["donor"].isin(paired_donors)].copy()
    excluded_sample_meta = samples[~samples["donor"].isin(paired_donors)].copy()
    qc_meta = vol_meta.query("sample_type == 'QC'").copy()
    blank_meta = vol_meta.query("sample_type == 'Blank'").copy()

    sample_columns = used_sample_meta["column"].tolist()
    qc_columns = qc_meta["column"].tolist()
    blank_columns = blank_meta["column"].tolist()

    X_sample = to_numeric_frame(feature_table[sample_columns])
    X_qc = to_numeric_frame(feature_table[qc_columns])
    X_blank = to_numeric_frame(feature_table[blank_columns])

    raw_feature_count = len(feature_table)
    mean_sample = X_sample.fillna(0).mean(axis=1)
    mean_blank = X_blank.fillna(0).mean(axis=1)
    blank_ratio = pd.Series(np.where(mean_sample > 0, mean_blank / mean_sample, np.inf), index=feature_table.index)
    pass_blank = blank_ratio < 0.10

    qc_presence_count = (X_qc.fillna(0) > 0).sum(axis=1)
    qc_presence_fraction = qc_presence_count / len(qc_columns)
    pass_qc_100 = qc_presence_count >= len(qc_columns)
    pass_qc_70 = qc_presence_count >= math.ceil(0.70 * len(qc_columns))

    sample_presence_count = (X_sample.fillna(0) > 0).sum(axis=1)
    sample_presence_fraction = sample_presence_count / len(sample_columns)
    pass_sample_100 = sample_presence_count >= len(sample_columns)
    pass_sample_70 = sample_presence_count >= math.ceil(0.70 * len(sample_columns))

    sample_col_to_timepoint = used_sample_meta.set_index("column")["timepoint"].to_dict()
    m1_columns = [col for col in sample_columns if sample_col_to_timepoint[col] == "M1"]
    m12_columns = [col for col in sample_columns if sample_col_to_timepoint[col] == "M12"]
    m1_presence_count = (X_sample[m1_columns].fillna(0) > 0).sum(axis=1)
    m12_presence_count = (X_sample[m12_columns].fillna(0) > 0).sum(axis=1)
    pass_group_70 = (m1_presence_count >= math.ceil(0.70 * len(m1_columns))) | (
        m12_presence_count >= math.ceil(0.70 * len(m12_columns))
    )

    candidate_mask = pass_blank & pass_qc_70 & pass_group_70
    candidate_qc = X_qc.loc[candidate_mask].fillna(0)
    qc_mean = candidate_qc.mean(axis=1)
    qc_sd = candidate_qc.std(axis=1, ddof=1)
    candidate_rsd = pd.Series(np.where(qc_mean > 0, qc_sd / qc_mean * 100, np.nan), index=candidate_qc.index)
    qc_rsd_full = pd.Series(np.nan, index=feature_table.index, dtype=float)
    qc_rsd_full.loc[candidate_rsd.index] = candidate_rsd
    pass_rsd_30 = qc_rsd_full <= 30
    pass_rsd_20 = qc_rsd_full <= 20
    final_mask = candidate_mask & pass_rsd_30
    accepted_mask = candidate_mask & pass_rsd_20
    borderline_mask = candidate_mask & (qc_rsd_full > 20) & (qc_rsd_full <= 30)
    unstable_mask = candidate_mask & (qc_rsd_full > 30)

    feature_metrics = pd.DataFrame({
        "feature_id": feature_table["feature_id"],
        "mass_avg": pd.to_numeric(feature_table.get("Mass (avg)"), errors="coerce"),
        "rt_avg_min": pd.to_numeric(feature_table.get("RT (avg)"), errors="coerce"),
        "height_avg": pd.to_numeric(feature_table.get("Height (avg)"), errors="coerce"),
        "blank_to_sample_ratio": blank_ratio,
        "qc_presence_count": qc_presence_count,
        "qc_presence_fraction": qc_presence_fraction,
        "sample_presence_count": sample_presence_count,
        "sample_presence_fraction": sample_presence_fraction,
        "m1_presence_count": m1_presence_count,
        "m12_presence_count": m12_presence_count,
        "qc_rsd_percent": qc_rsd_full,
        "pass_blank_10pct": pass_blank,
        "pass_qc_100pct": pass_qc_100,
        "pass_qc_70pct": pass_qc_70,
        "pass_sample_100pct": pass_sample_100,
        "pass_sample_70pct_all": pass_sample_70,
        "pass_group_70pct_at_least_one_timepoint": pass_group_70,
        "candidate_before_rsd": candidate_mask,
        "final_rsd_le_30": final_mask,
        "accepted_rsd_le_20": accepted_mask,
        "borderline_rsd_20_30": borderline_mask,
        "unstable_rsd_gt_30": unstable_mask,
    })

    final_features = feature_table.loc[final_mask].reset_index(drop=True)
    final_metrics = feature_metrics.loc[final_mask].reset_index(drop=True)
    final_feature_ids = final_features["feature_id"].tolist()

    final_matrix_raw = pd.concat([
        final_features[["feature_id", "Mass (avg)", "RT (avg)"]].reset_index(drop=True),
        X_sample.loc[final_mask].reset_index(drop=True).rename(columns=dict(zip(sample_columns, used_sample_meta["sample_id_inferred"].tolist(), strict=True))),
        X_qc.loc[final_mask].reset_index(drop=True).rename(columns=dict(zip(qc_columns, qc_meta["sample_id_inferred"].tolist(), strict=True))),
    ], axis=1)

    # Preprocess for PCA and statistical analysis: zero/missing -> half minimum positive per feature; log2; autoscale for PCA.
    object_columns = sample_columns + qc_columns
    object_meta = pd.concat([used_sample_meta, qc_meta], ignore_index=True).copy()
    object_labels = object_meta["sample_id_inferred"].tolist()
    object_matrix = to_numeric_frame(feature_table.loc[final_mask, object_columns]).copy()
    object_matrix.columns = object_labels
    imputed_matrix = impute_half_min_positive(object_matrix)
    log2_matrix = np.log2(imputed_matrix)
    scaled_matrix = pd.DataFrame(
        StandardScaler(with_mean=True, with_std=True).fit_transform(log2_matrix.T),
        index=object_labels,
        columns=final_feature_ids,
    )

    pca = PCA(n_components=min(5, scaled_matrix.shape[0], scaled_matrix.shape[1]), random_state=0)
    pca_scores_array = pca.fit_transform(scaled_matrix)
    pca_scores = pd.DataFrame(
        pca_scores_array,
        columns=[f"PC{i+1}" for i in range(pca_scores_array.shape[1])],
        index=object_labels,
    ).reset_index(names="object_id")
    pca_scores = pca_scores.merge(
        object_meta[["sample_id_inferred", "sample_type", "donor", "timepoint", "measurement_day", "column"]],
        left_on="object_id",
        right_on="sample_id_inferred",
        how="left",
    )
    pca_loadings = pd.DataFrame(
        pca.components_.T,
        index=final_feature_ids,
        columns=[f"PC{i+1}" for i in range(pca.n_components_)],
    ).reset_index(names="feature_id")

    # Donor-level aggregation to avoid pseudoreplication.
    sample_object_matrix = imputed_matrix[used_sample_meta["sample_id_inferred"].tolist()]
    sample_log2 = np.log2(sample_object_matrix)
    sample_meta_for_columns = used_sample_meta.set_index("sample_id_inferred")

    donor_timepoint_values = []
    for feature_id in final_feature_ids:
        row = {"feature_id": feature_id}
        feature_values = sample_log2.loc[sample_log2.index[final_feature_ids.index(feature_id)]]
        for donor in sorted(used_sample_meta["donor"].dropna().unique()):
            for timepoint in ["M1", "M12"]:
                cols = sample_meta_for_columns.query("donor == @donor and timepoint == @timepoint").index.tolist()
                values = feature_values[cols].astype(float)
                row[f"{donor}_{timepoint}"] = float(values.mean()) if len(values) else np.nan
        donor_timepoint_values.append(row)
    donor_aggregated = pd.DataFrame(donor_timepoint_values)

    statistics_rows = []
    donors = sorted(used_sample_meta["donor"].dropna().unique())
    for _, frow in donor_aggregated.iterrows():
        m1 = np.array([frow.get(f"{donor}_M1", np.nan) for donor in donors], dtype=float)
        m12 = np.array([frow.get(f"{donor}_M12", np.nan) for donor in donors], dtype=float)
        valid = np.isfinite(m1) & np.isfinite(m12)
        m1 = m1[valid]
        m12 = m12[valid]
        diff = m12 - m1
        mean_log2fc = float(np.nanmean(diff)) if len(diff) else np.nan
        median_log2fc = float(np.nanmedian(diff)) if len(diff) else np.nan
        statistics_rows.append({
            "feature_id": frow["feature_id"],
            "n_paired_donors": int(len(diff)),
            "mean_log2fc_M12_minus_M1": mean_log2fc,
            "median_log2fc_M12_minus_M1": median_log2fc,
            "wilcoxon_p_paired": safe_wilcoxon(m12, m1) if len(diff) >= 3 else np.nan,
            "mann_whitney_p_unpaired": safe_mann_whitney(m12, m1) if len(diff) >= 3 else np.nan,
            "mean_log2_M1": float(np.nanmean(m1)) if len(m1) else np.nan,
            "mean_log2_M12": float(np.nanmean(m12)) if len(m12) else np.nan,
        })
    stats_df = pd.DataFrame(statistics_rows)
    stats_df["wilcoxon_q_bh"] = benjamini_hochberg(stats_df["wilcoxon_p_paired"].to_numpy())
    stats_df["mann_whitney_q_bh"] = benjamini_hochberg(stats_df["mann_whitney_p_unpaired"].to_numpy())
    stats_df["abs_mean_log2fc"] = stats_df["mean_log2fc_M12_minus_M1"].abs()
    stats_df = stats_df.merge(final_metrics[["feature_id", "mass_avg", "rt_avg_min", "qc_rsd_percent", "borderline_rsd_20_30"]], on="feature_id", how="left")
    stats_df = stats_df.sort_values(["wilcoxon_q_bh", "wilcoxon_p_paired", "abs_mean_log2fc"], ascending=[True, True, False])

    # Total signal summaries.
    total_signal_objects = object_matrix.fillna(0).sum(axis=0)
    total_signal = pd.DataFrame({"object_id": total_signal_objects.index, "total_signal": total_signal_objects.values})
    total_signal = total_signal.merge(
        object_meta[["sample_id_inferred", "sample_type", "donor", "timepoint", "measurement_day"]],
        left_on="object_id",
        right_on="sample_id_inferred",
        how="left",
    )

    # Write processed data.
    vol_meta.to_csv(paths.processed_dir / "volume_column_metadata_inferred.csv", index=False)
    feature_metrics.to_csv(paths.processed_dir / "feature_filter_metrics.csv", index=False)
    final_metrics.to_csv(paths.processed_dir / "final_feature_metrics.csv", index=False)
    final_matrix_raw.to_csv(paths.processed_dir / "final_matrix_raw_volumes.csv", index=False)
    log2_matrix.reset_index(names="feature_index").to_csv(paths.processed_dir / "final_matrix_log2_imputed.csv", index=False)
    pca_scores.to_csv(paths.processed_dir / "pca_scores.csv", index=False)
    pca_loadings.to_csv(paths.processed_dir / "pca_loadings.csv", index=False)
    donor_aggregated.to_csv(paths.processed_dir / "donor_aggregated_log2_matrix.csv", index=False)
    stats_df.to_csv(paths.processed_dir / "univariate_statistics_m1_vs_m12.csv", index=False)
    total_signal.to_csv(paths.processed_dir / "total_signal_by_object.csv", index=False)

    # Filtering summary table.
    filtering_summary = pd.DataFrame([
        {"stage": "Raw feature table", "features": raw_feature_count, "removed_vs_previous": 0},
        {"stage": "After blank filter <10%", "features": int(pass_blank.sum()), "removed_vs_previous": int(raw_feature_count - pass_blank.sum())},
        {"stage": "After blank + 100% QC presence", "features": int((pass_blank & pass_qc_100).sum()), "removed_vs_previous": None},
        {"stage": "After blank + 70% QC presence", "features": int((pass_blank & pass_qc_70).sum()), "removed_vs_previous": int(pass_blank.sum() - (pass_blank & pass_qc_70).sum())},
        {"stage": "After blank + 100% Sample presence", "features": int((pass_blank & pass_sample_100).sum()), "removed_vs_previous": None},
        {"stage": "After blank + 70% all Sample presence", "features": int((pass_blank & pass_sample_70).sum()), "removed_vs_previous": None},
        {"stage": "After blank + 70% in at least one timepoint", "features": int((pass_blank & pass_group_70).sum()), "removed_vs_previous": None},
        {"stage": "Candidate: blank + 70% QC + 70% group", "features": int(candidate_mask.sum()), "removed_vs_previous": int((pass_blank & pass_qc_70).sum() - candidate_mask.sum())},
        {"stage": "Final: candidate + QC RSD <=30%", "features": int(final_mask.sum()), "removed_vs_previous": int(candidate_mask.sum() - final_mask.sum())},
        {"stage": "Strict accepted subset: QC RSD <=20%", "features": int(accepted_mask.sum()), "removed_vs_previous": int(candidate_mask.sum() - accepted_mask.sum())},
    ])
    filtering_summary.to_csv(paths.processed_dir / "filtering_summary.csv", index=False)

    # Metadata audit summary.
    exact_feature_files = set(vol_meta["file_name_inferred"].dropna())
    metadata_files = set(sample_metadata["File_name"].dropna())
    sample_ids_feature = set(vol_meta.query("sample_type == 'Sample'")["sample_id_inferred"].dropna())
    sample_ids_metadata = set(sample_metadata.query("Sample_type == 'Sample'")["Sample_ID"].dropna())
    metadata_audit = {
        "feature_table_rows": int(raw_feature_count),
        "feature_table_columns_after_drop_unnamed": int(feature_table.shape[1]),
        "volume_columns_total": int(len(vol_columns)),
        "volume_columns_sample_total": int((vol_meta["sample_type"] == "Sample").sum()),
        "volume_columns_sample_used": int(len(sample_columns)),
        "volume_columns_sample_excluded_unpaired_donor": int(len(excluded_sample_meta)),
        "volume_columns_qc": int(len(qc_columns)),
        "volume_columns_blank": int(len(blank_columns)),
        "metadata_sample_rows": int((sample_metadata["Sample_type"] == "Sample").sum()),
        "metadata_qc_rows": int((sample_metadata["Sample_type"] == "QC").sum()),
        "metadata_blank_rows": int((sample_metadata["Sample_type"] == "Blank").sum()),
        "feature_files_not_in_metadata_exact": sorted(exact_feature_files - metadata_files),
        "metadata_files_not_in_feature_exact": sorted(metadata_files - exact_feature_files),
        "feature_sample_ids_not_in_metadata": sorted(sample_ids_feature - sample_ids_metadata),
        "metadata_sample_ids_not_in_feature": sorted(sample_ids_metadata - sample_ids_feature),
        "sequence_rows": int(len(sequence)),
        "sequence_sample_type_counts": sequence["Sample_type"].value_counts(dropna=False).to_dict(),
    }
    with open(paths.processed_dir / "metadata_audit_summary.json", "w", encoding="utf-8") as handle:
        json.dump(metadata_audit, handle, ensure_ascii=False, indent=2)

    # Figures.
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 160,
        "savefig.dpi": 300,
    })

    # 01 filtering waterfall
    wf_stages = [
        ("Raw", raw_feature_count),
        ("Blank\n<10%", int(pass_blank.sum())),
        ("+ QC\n70%", int((pass_blank & pass_qc_70).sum())),
        ("+ group\n70%", int(candidate_mask.sum())),
        ("+ RSD\n≤30%", int(final_mask.sum())),
        ("RSD\n≤20%", int(accepted_mask.sum())),
    ]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.bar([s for s, _ in wf_stages], [v for _, v in wf_stages], color=["#385E72", "#4C7C8A", "#5A9E97", "#73BFA1", "#E3A34B", "#9B7EBD"])
    for i, (_, val) in enumerate(wf_stages):
        ax.text(i, val + raw_feature_count * 0.015, f"{val}", ha="center", va="bottom", fontweight="bold")
    ax.set_ylabel("Liczba cech")
    ax.set_title("Filtracja cech non-targeted")
    ax.set_ylim(0, raw_feature_count * 1.15)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "01_filtering_waterfall.png", bbox_inches="tight")
    plt.close(fig)

    # 02 RSD histogram
    rsd_values = qc_rsd_full.loc[candidate_mask].dropna()
    bins = [0, 10, 20, 30, np.inf]
    labels = ["0-10", "10-20", "20-30", ">30"]
    rsd_counts = pd.cut(rsd_values, bins=bins, labels=labels, right=True, include_lowest=True).value_counts().reindex(labels).fillna(0)
    fig, ax = plt.subplots(figsize=(7, 4.4))
    ax.bar(labels, rsd_counts.values, color=["#4C7C8A", "#73BFA1", "#E3A34B", "#C15D4B"])
    for i, val in enumerate(rsd_counts.values):
        ax.text(i, val + max(rsd_counts.values) * 0.03, f"{int(val)}", ha="center", va="bottom", fontweight="bold")
    ax.set_xlabel("QC RSD [%]")
    ax.set_ylabel("Liczba cech")
    ax.set_title("Stabilność analityczna cech w QC")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "02_qc_rsd_histogram.png", bbox_inches="tight")
    plt.close(fig)

    # 03 RT vs m/z feature map
    all_rt = pd.to_numeric(feature_table["RT (avg)"], errors="coerce")
    all_mz = pd.to_numeric(feature_table["Mass (avg)"], errors="coerce")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(all_rt, all_mz, s=12, alpha=0.28, label="raw", color="#A0A6AD", edgecolors="none")
    ax.scatter(all_rt[candidate_mask], all_mz[candidate_mask], s=16, alpha=0.65, label="po filtrach obecności", color="#4C7C8A", edgecolors="none")
    ax.scatter(all_rt[final_mask], all_mz[final_mask], s=22, alpha=0.9, label="final RSD ≤30%", color="#C15D4B", edgecolors="none")
    ax.set_xlabel("RT [min]")
    ax.set_ylabel("m/z / masa średnia")
    ax.set_title("Mapa cech LC-MS: RT vs m/z")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "03_rt_mz_feature_map.png", bbox_inches="tight")
    plt.close(fig)

    # 04 PCA sample/QC
    pc1_var, pc2_var = pca.explained_variance_ratio_[0] * 100, pca.explained_variance_ratio_[1] * 100
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    groups = [
        ("QC", "QC", "#2D5F73", "s"),
        ("Sample M1", "M1", "#D27C2C", "o"),
        ("Sample M12", "M12", "#8E5C9E", "o"),
    ]
    for label, tp, color, marker in groups:
        if tp == "QC":
            mask = pca_scores["sample_type"] == "QC"
        else:
            mask = (pca_scores["sample_type"] == "Sample") & (pca_scores["timepoint"] == tp)
        ax.scatter(pca_scores.loc[mask, "PC1"], pca_scores.loc[mask, "PC2"], label=label, s=60, color=color, marker=marker, alpha=0.85, edgecolor="white", linewidth=0.6)
    for _, row in pca_scores.query("sample_type == 'QC'").iterrows():
        ax.text(row["PC1"], row["PC2"], str(row["measurement_day"]).replace("day", "d"), fontsize=7, ha="left", va="bottom")
    ax.axhline(0, color="#CCCCCC", linewidth=0.8)
    ax.axvline(0, color="#CCCCCC", linewidth=0.8)
    ax.set_xlabel(f"PC1 ({pc1_var:.1f}% wariancji)")
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}% wariancji)")
    ax.set_title("PCA po log2 i autoskalowaniu")
    ax.legend(frameon=False, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "04_pca_samples_qc.png", bbox_inches="tight")
    plt.close(fig)

    # 05 PCA biological samples by donor/timepoint
    fig, ax = plt.subplots(figsize=(7.2, 5.1))
    sample_scores = pca_scores.query("sample_type == 'Sample'").copy()
    colors = {"M1": "#D27C2C", "M12": "#8E5C9E"}
    for tp in ["M1", "M12"]:
        mask = sample_scores["timepoint"] == tp
        ax.scatter(sample_scores.loc[mask, "PC1"], sample_scores.loc[mask, "PC2"], s=58, label=tp, color=colors[tp], alpha=0.78, edgecolor="white", linewidth=0.5)
    for donor in sorted(sample_scores["donor"].dropna().unique()):
        donor_points = sample_scores.query("donor == @donor").groupby("timepoint")[["PC1", "PC2"]].mean()
        if {"M1", "M12"}.issubset(donor_points.index):
            ax.plot([donor_points.loc["M1", "PC1"], donor_points.loc["M12", "PC1"]], [donor_points.loc["M1", "PC2"], donor_points.loc["M12", "PC2"]], color="#666666", linewidth=0.8, alpha=0.55)
            ax.text(donor_points.mean()["PC1"], donor_points.mean()["PC2"], donor, fontsize=8, ha="center", va="center")
    ax.axhline(0, color="#CCCCCC", linewidth=0.8)
    ax.axvline(0, color="#CCCCCC", linewidth=0.8)
    ax.set_xlabel(f"PC1 ({pc1_var:.1f}%)")
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}%)")
    ax.set_title("PCA próbek biologicznych: M1 vs M12")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "05_pca_biological_timepoint.png", bbox_inches="tight")
    plt.close(fig)

    # 06 volcano
    volcano = stats_df.copy()
    volcano["neg_log10_p"] = -np.log10(volcano["wilcoxon_p_paired"].replace(0, np.nextafter(0, 1)))
    is_fdr = volcano["wilcoxon_q_bh"] < 0.05
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.scatter(volcano["mean_log2fc_M12_minus_M1"], volcano["neg_log10_p"], s=34, alpha=0.65, color="#4C7C8A", edgecolors="none", label="features")
    if is_fdr.any():
        ax.scatter(volcano.loc[is_fdr, "mean_log2fc_M12_minus_M1"], volcano.loc[is_fdr, "neg_log10_p"], s=45, color="#C15D4B", edgecolors="none", label="FDR < 0.05")
    ax.axvline(-1, color="#999999", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="#999999", linestyle="--", linewidth=0.8)
    ax.axhline(-math.log10(0.05), color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("średni log2FC M12 - M1")
    ax.set_ylabel("-log10(p), Wilcoxon parowany")
    ax.set_title("Analiza różnic M12 vs M1")
    ax.text(0.02, 0.96, f"FDR<0.05: {int(is_fdr.sum())}", transform=ax.transAxes, va="top", ha="left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "06_volcano_m12_vs_m1.png", bbox_inches="tight")
    plt.close(fig)

    # 07 heatmap top features by p-value/FC
    top_n = min(25, len(stats_df))
    top_ids = stats_df.sort_values(["wilcoxon_p_paired", "abs_mean_log2fc"], ascending=[True, False]).head(top_n)["feature_id"].tolist()
    donor_cols = [c for c in donor_aggregated.columns if re.match(r"^[A-Za-z][A-Za-z0-9]*_M(1|12)$", c)]
    heat = donor_aggregated.set_index("feature_id").loc[top_ids, donor_cols].astype(float)
    heat_z = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(9.2, max(4.8, 0.26 * top_n + 1.6)))
    im = ax.imshow(heat_z.values, aspect="auto", cmap="coolwarm", vmin=-2.2, vmax=2.2)
    ax.set_xticks(np.arange(len(donor_cols)))
    ax.set_xticklabels(donor_cols, rotation=45, ha="right")
    ax.set_yticks(np.arange(top_n))
    ax.set_yticklabels(top_ids)
    ax.set_title("Top cechy wg testu parowanego: log2 z-score na poziomie dawczyni")
    ax.set_xlabel("Dawczyni_timepoint")
    ax.set_ylabel("Feature ID")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("z-score")
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "07_heatmap_top_features.png", bbox_inches="tight")
    plt.close(fig)

    # 08 total signal distribution
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    groups_to_plot = [
        total_signal.query("sample_type == 'QC'")["total_signal"].values,
        total_signal.query("sample_type == 'Sample' and timepoint == 'M1'")["total_signal"].values,
        total_signal.query("sample_type == 'Sample' and timepoint == 'M12'")["total_signal"].values,
    ]
    ax.boxplot(groups_to_plot, tick_labels=["QC", "M1", "M12"], showmeans=True)
    for i, vals in enumerate(groups_to_plot, start=1):
        jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) else []
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=30, alpha=0.65, color="#4C7C8A", edgecolor="white", linewidth=0.4)
    ax.set_ylabel("Suma wolumenów cech finalnych")
    ax.set_title("Total signal po filtracji")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "08_total_signal_boxplot.png", bbox_inches="tight")
    plt.close(fig)

    # 09 QC total signal ordered by inferred day and suffix
    qc_total = total_signal.query("sample_type == 'QC'").copy()
    def qc_sort_key(obj_id: str) -> Tuple[int, int]:
        match = re.search(r"QC_day_([0-9]+)_([0-9]+)", obj_id)
        return (int(match.group(1)), int(match.group(2))) if match else (999, 999)
    qc_total["sort_key"] = qc_total["object_id"].apply(qc_sort_key)
    qc_total = qc_total.sort_values("sort_key")
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.plot(range(1, len(qc_total) + 1), qc_total["total_signal"].values, marker="o", linewidth=1.6, color="#2D5F73")
    ax.set_xticks(range(1, len(qc_total) + 1))
    ax.set_xticklabels(qc_total["object_id"].str.replace("QC_day_", "d", regex=False), rotation=45, ha="right")
    ax.set_xlabel("QC w kolejności dnia/sufiksu")
    ax.set_ylabel("Total signal")
    ax.set_title("Dryf łącznego sygnału w QC")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.figures_dir / "09_qc_total_signal_drift.png", bbox_inches="tight")
    plt.close(fig)

    # Create a compact JSON summary for document generation.
    qc_rsd_counts = {
        "le20": int((qc_rsd_full.loc[candidate_mask] <= 20).sum()),
        "20to30": int(((qc_rsd_full.loc[candidate_mask] > 20) & (qc_rsd_full.loc[candidate_mask] <= 30)).sum()),
        "gt30": int((qc_rsd_full.loc[candidate_mask] > 30).sum()),
    }
    top_candidates = stats_df.head(10)[[
        "feature_id", "mass_avg", "rt_avg_min", "mean_log2fc_M12_minus_M1", "wilcoxon_p_paired", "wilcoxon_q_bh", "qc_rsd_percent"
    ]].to_dict(orient="records")
    fdr_significant = int((stats_df["wilcoxon_q_bh"] < 0.05).sum())
    pca_var = [float(x) for x in pca.explained_variance_ratio_]

    summary = {
        "input_files_used": [feature_path.name, metadata_path.name],
        "targeted_analysis": "omitted_on_user_instruction",
        "raw_feature_count": raw_feature_count,
        "feature_table_column_count_after_drop_unnamed": int(feature_table.shape[1]),
        "volume_columns_total": int(len(vol_columns)),
        "sample_columns_used": int(len(sample_columns)),
        "qc_columns_used": int(len(qc_columns)),
        "blank_columns_used": int(len(blank_columns)),
        "excluded_unpaired_donor_sample_columns": int(len(excluded_sample_meta)),
        "missing_sample_cells_raw": int(X_sample.isna().sum().sum()),
        "missing_qc_cells_raw": int(X_qc.isna().sum().sum()),
        "filter_counts": {row["stage"]: int(row["features"]) for _, row in filtering_summary.iterrows()},
        "qc_rsd_counts_candidate": qc_rsd_counts,
        "final_feature_count": int(final_mask.sum()),
        "accepted_feature_count_rsd_le20": int(accepted_mask.sum()),
        "borderline_feature_count_rsd_20_30": int(borderline_mask.sum()),
        "unstable_candidate_count_rsd_gt30": int(unstable_mask.sum()),
        "pca_explained_variance_ratio": pca_var,
        "pca_pc1_pc2_percent": [float(pc1_var), float(pc2_var)],
        "fdr_significant_features_wilcoxon_q_lt_0_05": fdr_significant,
        "top_candidate_features": top_candidates,
        "metadata_audit": metadata_audit,
        "notes": [
            "Biological sample columns were linked by inferred Sample_ID because several exact file names differ between feature table headers and metadata.",
            "Columns of one donor without M1 samples were excluded (paired M1 vs M12 design).",
            "PCA was performed on log2-imputed values followed by feature-wise autoscaling (z-score).",
            "Primary M1 vs M12 statistics used donor-level aggregation and paired Wilcoxon tests to avoid pseudoreplication. Mann-Whitney U was also exported only as a worksheet-compatible sensitivity test.",
        ],
    }
    with open(paths.processed_dir / "results_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run non-targeted lipidomics analysis.")
    parser.add_argument("--input-dir", type=Path, default=None, help="Directory with input CSV/XLSX data.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for figures and processed_data.")
    args = parser.parse_args()
    paths = prepare_paths(args.input_dir, args.output_dir)
    run_analysis(paths)


if __name__ == "__main__":
    main()
