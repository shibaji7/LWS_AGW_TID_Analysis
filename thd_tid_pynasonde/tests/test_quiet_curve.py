"""Offline checks for quiet_curve.py (bin/median logic; no network, no OMNI)."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import quiet_curve as qcurve

C = dict(qcurve.DEFAULTS, bin_min=15.0, min_days=3)


def day(offset, hours=24, step=5):
    t = pd.date_range("2025-06-01", periods=hours * 60 // step, freq=f"{step}min")
    return pd.DataFrame({"datetime": t, "foF2": 5.0 + offset + (t.hour + t.minute / 60) * 0.1})


class BinTests(unittest.TestCase):
    def test_bin_day_averages_within_bin(self):
        s = qcurve.bin_day(day(0), "foF2", 15.0)
        self.assertEqual(len(s), 96)
        self.assertEqual(s.index[1], 0.25)
        self.assertAlmostEqual(s.iloc[0], 5.0 + 0.1 * (5 / 60), places=6)  # mean of 00:00,05:00,10:00

    def test_median_is_robust_to_one_outlier_day(self):
        series = {i: qcurve.bin_day(day(o), "foF2", 15.0) for i, o in enumerate([0.0, 0.1, -0.1, 0.05, 5.0])}
        cv = qcurve.quiet_curve(series, dict(C, min_days=5))
        base = qcurve.bin_day(day(0), "foF2", 15.0).to_numpy()
        np.testing.assert_allclose(cv["median"].to_numpy(), base + 0.05, atol=1e-9)
        self.assertTrue((cv.n_days == 5).all())

    def test_bins_with_too_few_days_are_nan(self):
        short = {i: qcurve.bin_day(day(0, hours=6 if i else 24), "foF2", 15.0) for i in range(4)}
        cv = qcurve.quiet_curve({0: short[0], 1: short[1]}, C)  # 2 days < min_days=3
        self.assertTrue(cv["median"].isna().all())
        cv = qcurve.quiet_curve(short, C)  # 4 days early, 1 day late
        self.assertTrue(cv.loc[cv.ut_hour < 6, "median"].notna().all())
        self.assertTrue(cv.loc[cv.ut_hour > 7, "median"].isna().all())


if __name__ == "__main__":
    unittest.main()
