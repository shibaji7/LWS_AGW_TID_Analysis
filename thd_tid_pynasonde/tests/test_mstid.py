"""Offline checks for mstid.py (baseline interpolation, event detection)."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mstid

C = dict(min_event_min=60.0, edge_min=30.0)


def curve(values):
    x = np.arange(len(values)) * 24.0 / len(values) + 12.0 / len(values)
    return pd.DataFrame(dict(ut_hour=x, median=values, q25=values, q75=values))


class BaselineTests(unittest.TestCase):
    def test_interpolates_and_wraps_midnight(self):
        c = curve(np.array([0.0, 6.0, 12.0, 18.0]))  # bin centres 3, 9, 15, 21 h
        t = pd.Series(pd.to_datetime(["2025-06-01 06:00", "2025-06-01 00:00"]))
        out = mstid.baseline_at(c, t)["median"]
        self.assertAlmostEqual(out[0], 3.0)           # halfway between 0 (3 h) and 6 (9 h)
        self.assertAlmostEqual(out[1], 18.0 + (0.0 - 18.0) * 0.5)  # wrap: 21 h -> 3 h next day

    def test_nan_bin_propagates(self):
        c = curve(np.array([0.0, np.nan, 12.0, 18.0]))
        out = mstid.baseline_at(c, pd.Series(pd.to_datetime(["2025-06-01 07:00", "2025-06-01 16:00"])))["median"]
        self.assertTrue(np.isnan(out[0]))
        self.assertFalse(np.isnan(out[1]))


class EventTests(unittest.TestCase):
    def packet(self, amp, period=30.0, cadence=5.0, hours=12.0, start_h=4.0, dur_h=3.0):
        idx = pd.date_range("2025-06-01", periods=int(hours * 60 / cadence), freq=f"{cadence}min")
        m = (idx - idx[0]).total_seconds().to_numpy() / 60.0
        env = ((m >= start_h * 60) & (m < (start_h + dur_h) * 60)).astype(float)
        return pd.Series(amp * env * np.sin(2 * np.pi * m / period), index=idx)

    def test_detects_wave_packet_and_period(self):
        ev = mstid.detect_events(self.packet(1.0), amp_thr=0.5, cadence=5.0, c=C)
        self.assertEqual(len(ev), 1)
        self.assertAlmostEqual(ev[0]["period_min"], 30.0, delta=2.0)
        self.assertAlmostEqual(ev[0]["duration_min"], 180.0, delta=30.0)
        self.assertGreater(ev[0]["amplitude"], 0.8)

    def test_below_threshold_or_too_short_gives_nothing(self):
        self.assertEqual(mstid.detect_events(self.packet(0.3), 0.5, 5.0, C), [])
        self.assertEqual(mstid.detect_events(self.packet(1.0, dur_h=0.5), 0.5, 5.0, C), [])

    def test_nan_gap_splits_events_and_edges_are_ignored(self):
        s = self.packet(1.0, hours=14.0, start_h=2.0, dur_h=10.0)
        self.assertEqual(mstid.detect_events(s, np.nan, 5.0, C), [])  # NaN threshold is safe
        s.iloc[72:84] = np.nan  # blank 06:00-07:00 in the middle of the packet
        ev = mstid.detect_events(s, 0.5, 5.0, C)
        self.assertEqual(len(ev), 2)  # one event each side of the blanked hour
        self.assertTrue(all(e["end"] <= pd.Timestamp("2025-06-01 06:00") or e["start"] >= pd.Timestamp("2025-06-01 07:00") for e in ev))


class TerminatorTests(unittest.TestCase):
    def test_flags_sunrise_window_not_midday(self):
        lat, lon = 45.0, -83.5  # northern Michigan, June: sunrise ~09:50 UT
        rise = pd.Timestamp("2025-06-01 09:30")
        self.assertTrue(mstid.near_terminator(rise, rise + pd.Timedelta(hours=1), lat, lon, 60.0))
        noon = pd.Timestamp("2025-06-01 17:00")  # local solar noon ~17:34 UT
        self.assertFalse(mstid.near_terminator(noon, noon + pd.Timedelta(hours=1), lat, lon, 60.0))


if __name__ == "__main__":
    unittest.main()
