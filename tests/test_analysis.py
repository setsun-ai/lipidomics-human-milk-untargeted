"""Tests for analysis_non_targeted.py on synthetic data (the real input data are restricted)."""
import json
import re

import numpy as np
import pandas as pd
import pytest

import analysis_non_targeted as an
from tests.conftest import PAIRED, UNPAIRED


def test_parse_volume_label_variants():
    s = an.parse_volume_label("[Vol] AB 12_2B 1234")
    assert (s["sample_type"], s["donor"], s["timepoint"]) == ("Sample", "AB", "M12")
    assert (s["bio_rep"], s["extraction_rep"]) == (2, "B")
    assert an.parse_volume_label("[Vol] D3 1_1A 976")["donor"] == "D3"  # aliases contain digits
    assert an.parse_volume_label("[Vol] QC_day 8_9")["sample_type"] == "QC"
    assert an.parse_volume_label("[Vol] Blank extraction _day 9")["sample_type"] == "Blank"


def test_benjamini_hochberg_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.20])
    # step-up BH: sorted p * m / rank, then cumulative minimum from the top
    np.testing.assert_allclose(an.benjamini_hochberg(p), [0.04, 0.0533333, 0.0533333, 0.2], rtol=1e-5)


def test_half_min_imputation_uses_feature_minimum():
    m = pd.DataFrame({"a": [4.0, np.nan], "b": [2.0, 8.0]}, index=["f1", "f2"])
    out = an.impute_half_min_positive(m)
    assert out.notna().all().all()
    assert out.loc["f2", "a"] == pytest.approx(4.0)  # half of the f2 minimum (8.0)


def test_pseudonymisation_hides_original_codes():
    ft = pd.DataFrame(columns=["[Vol] AB 1_1A 1", "[Vol] EFG 12_1A 2", "[Vol] QC_day 6_1"])
    meta = pd.DataFrame({"Sample_ID": ["AB_M1_B1_EA"], "File_name": ["AB_1_1A_1.mzML"], "Donor_ID": ["AB"]})
    seq = pd.DataFrame({"Sample_ID": ["AB_M1_B1_EA"], "File_name": ["AB_1_1A_1.mzML"]})
    ft2, meta2, seq2 = an.pseudonymise_donors(ft, meta, seq)
    text = " ".join(ft2.columns) + meta2.to_string() + seq2.to_string()
    assert not re.search(r"(?<![A-Za-z])(AB|EFG)(?=[_\s])", text)
    assert "[Vol] D1 1_1A 1" in ft2.columns and "QC_day 6_1" in " ".join(ft2.columns)


@pytest.fixture(scope="module")
def run(synthetic_input, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    an.run_analysis(an.prepare_paths(synthetic_input, out))
    return out


def test_end_to_end_outputs(run):
    for name in ["01_filtering_waterfall.png", "05_pca_biological_timepoint.png", "07_heatmap_top_features.png"]:
        assert (run / "figures" / name).stat().st_size > 0
    summary = json.loads((run / "processed_data" / "results_summary.json").read_text(encoding="utf-8"))
    assert summary["sample_columns_used"] == 2 * 4 * len(PAIRED)  # 2 timepoints x 4 replicates
    assert summary["excluded_unpaired_donor_sample_columns"] == 4
    assert summary["qc_columns_used"] == 10
    assert 0 < summary["final_feature_count"] <= 60


def test_no_original_donor_codes_in_any_output(run):
    codes = "|".join(PAIRED + [UNPAIRED])
    pattern = re.compile(rf"(?<![A-Za-z])({codes})(?=[_\s\"])")
    for f in run.rglob("*"):
        if f.suffix in {".csv", ".json", ".txt"}:
            assert not pattern.search(f.read_text(encoding="utf-8")), f.name


def test_donor_level_statistics_are_paired(run):
    stats = pd.read_csv(run / "processed_data" / "univariate_statistics_m1_vs_m12.csv")
    assert (stats["n_paired_donors"] == len(PAIRED)).all()
    assert stats["wilcoxon_q_bh"].between(0, 1).all()
