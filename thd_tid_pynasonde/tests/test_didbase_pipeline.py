"""Offline checks: python -m unittest discover -s tests -v (run from thd_tid_pynasonde/)."""
import datetime as dt
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import didbase_pipeline as pipeline


class DownloadTests(unittest.TestCase):
    def test_seasons_and_future_dates(self):
        days = list(pipeline.requested_dates([2024, 2025, 2026, 2026], dt.date(2027, 1, 1)))
        self.assertEqual(len(days), 549)
        self.assertEqual(days[0], dt.date(2024, 4, 1))
        self.assertEqual(days[-1], dt.date(2026, 9, 30))
        self.assertTrue(all(4 <= day.month <= 9 for day in days))
        self.assertEqual(list(pipeline.requested_dates([2026], dt.date(2026, 3, 31))), [])
        self.assertEqual(len(list(pipeline.requested_dates([2026], dt.date(2026, 4, 2)))), 2)

    def test_geography_and_catalogue_parser(self):
        rows = [['URSI', 'Country', 'Lat', 'Long'],
                ['AU930', 'USA', 30.4, 262.3], ['AU930', 'USA', 30.4, 262.3],
                ['CA001', 'Canada', 49, -100], ['MX001', 'Mexico', 30, -110],
                ['AK001', 'USA', 65, -150], ['HI001', 'USA', 21, -157],
                ['PR001', 'USA', 18, -66], ['NO001', 'USA', 40, -100]]
        catalog = pipeline.parse_station_html('var jsonData = ' + json.dumps(rows) + ';')
        available = pd.DataFrame({'URSI': ['AU930', 'CA001', 'MX001', 'AK001', 'HI001', 'PR001']})
        selected = pipeline.select_conus_stations(available, catalog)
        self.assertEqual(selected.URSI.tolist(), ['AU930'])
        self.assertAlmostEqual(selected.LONG.iloc[0], -97.7)

    def test_resume_empty_partial_corrupt_and_failed_days(self):
        now = dt.datetime(2026, 9, 29, 12)
        day = dt.date(2024, 4, 1)
        db = Mock()
        db.fetch_scaled_params.return_value = pd.DataFrame({
            'datetime': [dt.datetime.combine(day, dt.time.max)],
            'foF2': [8.123456789], 'scaler_id': [1], 'foF2_qletter': ['A']})
        with TemporaryDirectory() as directory:
            root = Path(directory)
            def download(date=day, **kwargs):
                return pipeline.download_day(db, 'AU930', date, ['foF2'], root, now=now, **kwargs)
            self.assertEqual(download(), 'saved')
            self.assertEqual(db.fetch_scaled_params.call_args.kwargs['end'], dt.datetime.combine(day, dt.time.max))
            path = root / '2024/AU930/2024-04-01.csv'
            frame = pd.read_csv(path)
            self.assertEqual(frame.foF2.iloc[0], 8.123456789)
            self.assertEqual(frame.foF2_qletter.iloc[0], 'A')
            self.assertEqual(download(), 'cached')
            self.assertEqual(db.fetch_scaled_params.call_count, 1)
            self.assertEqual(download(overwrite=True), 'saved')
            path.write_text('truncated')
            self.assertEqual(download(), 'saved')
            db.fetch_scaled_params.return_value = pd.DataFrame()
            empty_day = dt.date(2024, 4, 2)
            self.assertEqual(download(empty_day), 'empty')
            self.assertEqual(download(empty_day), 'cached')
            self.assertEqual(download(now.date()), 'empty')
            self.assertEqual(download(now.date()), 'empty')
            self.assertEqual(db.fetch_scaled_params.call_args.kwargs['end'], now)
            db.fetch_scaled_params.side_effect = RuntimeError('connection lost')
            with self.assertRaises(RuntimeError):
                download(dt.date(2024, 4, 3))
            self.assertFalse((root / '2024/AU930/2024-04-03.json').exists())


if __name__ == '__main__':
    unittest.main()
