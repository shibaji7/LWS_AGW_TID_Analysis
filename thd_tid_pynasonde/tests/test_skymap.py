"""Offline checks for skymap.py (wavelet, filters, tilt)."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import skymap

C = dict(skymap.scfg())


class WaveletTests(unittest.TestCase):
    def test_ridge_at_sinusoid_period(self):
        dt = 7.5
        t = np.arange(0, 12 * 60, dt)
        y = np.sin(2 * np.pi * t / 30.0)
        periods = np.geomspace(10, 120, 60)
        W, coi = skymap.morse_cwt(y, dt, periods, 3.0, 20.0)
        mid = np.abs(W[:, len(t) // 3: 2 * len(t) // 3]) ** 2
        self.assertAlmostEqual(periods[int(np.argmax(mid.mean(axis=1)))], 30.0, delta=3.0)
        self.assertTrue((np.diff(coi) > 0).all())  # cone widens with period

    def test_ridge_in_window_picks_burst_period(self):
        idx = pd.date_range("2025-06-01", periods=192, freq="7.5min")
        t = np.arange(192) * 7.5
        y = np.where((t > 400) & (t < 900), np.sin(2 * np.pi * t / 45.0), 0.0) + 0.2 * np.sin(2 * np.pi * t / 20.0)
        grid, periods, power, coi, _ = skymap.wavelet_of(pd.Series(y, index=idx), C)
        P, _ = skymap.ridge_in_window(grid, periods, power, coi, idx[70], idx[110])
        self.assertAlmostEqual(P, 45.0, delta=6.0)


class FilterTests(unittest.TestCase):
    def echoes(self):
        rng = np.random.default_rng(0)
        n = 200
        d = pd.DataFrame(dict(file_time=pd.Timestamp("2025-06-01 12:00"), group_range=300.0,
                              x_coord=rng.normal(2, 3, n), y_coord=rng.normal(-1, 3, n),
                              spect_dop=rng.normal(5, 1, n)))
        extra = pd.DataFrame(dict(file_time=pd.Timestamp("2025-06-01 12:00"),
                                  group_range=[100.0, 600.0, 300.0, 300.0], x_coord=[0, 0, 24.0, 1.0],
                                  y_coord=[0, 0, 24.0, 1.0], spect_dop=[5, 5, 5, 40.0]))
        return pd.concat([d, extra], ignore_index=True)

    def test_range_zenith_and_doppler_filters(self):
        kept = skymap.filter_echoes(self.echoes(), C)
        self.assertFalse(kept.group_range.lt(150).any() or kept.group_range.gt(450).any())
        self.assertTrue((kept.zenith <= 25).all())
        self.assertFalse((kept.spect_dop > 20).any())
        self.assertGreater(len(kept), 150)

    def test_tilt_is_echo_centre_and_needs_min_echoes(self):
        kept = skymap.filter_echoes(self.echoes(), C)
        t = skymap.tilt_series(kept, C)
        self.assertAlmostEqual(t.ew.iloc[0], 2.0, delta=0.8)
        self.assertAlmostEqual(t.ns.iloc[0], -1.0, delta=0.8)
        t2 = skymap.tilt_series(kept.head(10), C)
        self.assertTrue(np.isnan(t2.ew.iloc[0]))

    def test_azimuth_from_north_clockwise(self):
        g = skymap.add_geometry(pd.DataFrame(dict(x_coord=[0.0, 5.0, 0.0], y_coord=[5.0, 0.0, -5.0])))
        np.testing.assert_allclose(g.azimuth, [0.0, 90.0, 180.0])


if __name__ == "__main__":
    unittest.main()
