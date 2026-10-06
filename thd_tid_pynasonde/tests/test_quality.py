"""Offline checks for quality.py (QC + band-pass)."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import quality


def make_day(cadence=7.5, hours=24, fn=None, seed=0):
    t = pd.date_range("2025-06-01", periods=int(hours * 60 / cadence), freq=f"{cadence}min")
    x = np.arange(len(t)) * cadence
    y = (8 + np.random.default_rng(seed).normal(0, 0.03, len(t))) if fn is None else fn(x)
    return pd.DataFrame({"datetime": t, "foF2": y})


class CleanDayTests(unittest.TestCase):
    def test_clean_day_passes(self):
        seg, m = quality.clean_day(make_day(), "foF2")
        self.assertEqual(quality.evaluate(m, m["noise"]), "")
        self.assertEqual(m["cadence_min"], 7.5)
        self.assertGreater(m["seg_hours"], 23)

    def test_spike_removed_and_out_of_range_dropped(self):
        df = make_day()
        df.loc[50, "foF2"] = 15.0   # in range, but a spike
        df.loc[60, "foF2"] = 99.0   # out of physical range
        seg, m = quality.clean_day(df, "foF2")
        self.assertEqual(m["n_spikes"], 1)
        self.assertLess(seg.foF2.max(), 9)

    def test_keeps_longest_segment_across_gap(self):
        df = make_day()
        df = df.drop(index=range(40, 60))  # 2.5 h hole; segments of ~5 h and ~16.5 h remain
        seg, m = quality.clean_day(df, "foF2")
        self.assertGreater(m["seg_hours"], 16)
        self.assertEqual(seg.datetime.iloc[0], pd.Timestamp("2025-06-01 07:30"))  # first sample after the hole

    def test_short_day_fails(self):
        seg, m = quality.clean_day(make_day(hours=3), "foF2")
        self.assertIn("segment", quality.evaluate(m, m["noise"]))

    def test_noisy_day_fails_relative_to_station(self):
        seg, m = quality.clean_day(make_day(), "foF2")
        noisy = dict(m, noise=m["noise"] * 10)
        self.assertIn("noise", quality.evaluate(noisy, m["noise"]))


class BandpassTests(unittest.TestCase):
    def test_passes_30min_wave_rejects_trend_and_fast_noise(self):
        wave = lambda x: 8 + 0.002 * x + 0.3 * np.sin(2 * np.pi * x / 30) + 0.3 * np.sin(2 * np.pi * x / 360)
        seg = make_day(cadence=5, fn=wave)
        f, band = quality.bandpass(seg, "foF2", 5.0)
        self.assertEqual(band, (15.0, 60.0))
        mid = f.iloc[len(f) // 4: -len(f) // 4]
        self.assertAlmostEqual(mid.std() * np.sqrt(2), 0.3, delta=0.05)  # only the 30-min wave remains

    def test_band_respects_nyquist(self):
        self.assertEqual(quality.effective_band(15.0), (33.0, 60.0))
        self.assertEqual(quality.effective_band(7.5), (16.5, 60.0))

    def test_masks_output_inside_long_gap_only(self):
        seg = make_day(cadence=5, fn=lambda x: 8 + 0.3 * np.sin(2 * np.pi * x / 30))
        seg = seg.drop(index=range(100, 112))  # 60-min hole (> 30-min mask limit)
        f, _ = quality.bandpass(seg, "foF2", 5.0)
        hole = f[(f.index > seg.datetime.iloc[99]) & (f.index < seg.datetime.iloc[100])]
        self.assertTrue(len(hole) > 0 and hole.isna().all())
        self.assertTrue(f[f.index <= seg.datetime.iloc[99]].notna().all())
        self.assertTrue(f[f.index >= seg.datetime.iloc[100]].notna().all())

    def test_short_gap_not_masked(self):
        seg = make_day(cadence=5, fn=lambda x: 8 + 0.3 * np.sin(2 * np.pi * x / 30))
        f, _ = quality.bandpass(seg.drop(index=[100, 101]), "foF2", 5.0)  # 15-min hole
        self.assertTrue(f.notna().all())

    def test_too_short_segment_returns_none(self):
        f, _ = quality.bandpass(make_day(cadence=15, hours=3), "foF2", 15.0)
        self.assertIsNone(f)


if __name__ == "__main__":
    unittest.main()
