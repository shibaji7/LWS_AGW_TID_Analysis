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

    def test_sky_download_cache_and_event_days(self):
        now = dt.datetime(2026, 10, 6, 12)
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            db = Mock()
            db.fetch_sky_files.return_value = [root / 'a.SKY', root / 'b.SKY']
            day = dt.date(2024, 4, 2)
            self.assertEqual(pipeline.download_day_sky(db, 'AL945', day, root, now=now), 'saved')
            self.assertEqual(pipeline.download_day_sky(db, 'AL945', day, root, now=now), 'cached')
            self.assertEqual(db.fetch_sky_files.call_count, 1)
            self.assertEqual(json.loads((root / '2024/AL945/sky/2024-04-02.json').read_text())['files'], 2)
            db.fetch_sky_files.return_value = []
            self.assertEqual(pipeline.download_day_sky(db, 'AL945', dt.date(2024, 4, 3), root, now=now), 'empty')

            csv = root / 'days.csv'
            pd.DataFrame(dict(station=['AL945', 'AL945', 'AL945', 'MHJ45', 'AL945'],
                              date=['2024-04-02', '2024-04-02', '2024-04-05', '2024-05-01', '2024-11-01'],
                              enhanced=[True, True, False, True, True])).to_csv(csv, index=False)
            self.assertEqual(pipeline.event_station_days(csv), [('AL945', dt.date(2024, 4, 2)), ('MHJ45', dt.date(2024, 5, 1))])
            self.assertEqual(pipeline.event_station_days(csv, stations=['AL945']), [('AL945', dt.date(2024, 4, 2))])

    def test_sky_explicit_dates_ignore_season_and_years(self):
        calls = []
        class FakeDb:
            def __init__(self, *a): pass
            def __enter__(self): return self
            def __exit__(self, *a): return False
        import sys, types
        mod = types.ModuleType('pynasonde.digisonde.driftbase')
        mod.DriftBaseConnector = FakeDb
        orig = pipeline.download_day_sky
        pipeline.download_day_sky = lambda db, st, day, out, **k: (calls.append((st, day)) or 'saved')
        old = {k: sys.modules.get(k) for k in ('pynasonde.digisonde.driftbase',)}
        sys.modules['pynasonde.digisonde.driftbase'] = mod
        try:
            with TemporaryDirectory() as tmp:
                pipeline.run(sky=True, stations=['mhj45'], output_dir=tmp, dates=['2023-01-15', '2023-04-24'], retries=1)
        finally:
            pipeline.download_day_sky = orig
            sys.modules['pynasonde.digisonde.driftbase'] = old['pynasonde.digisonde.driftbase']
        self.assertEqual(calls, [('MHJ45', dt.date(2023, 1, 15)), ('MHJ45', dt.date(2023, 4, 24))])

    def test_months_override_and_combine(self):
        days = list(pipeline.requested_dates([2026], dt.date(2026, 10, 25), 1, 12))
        self.assertEqual((days[0], days[-1]), (dt.date(2026, 1, 1), dt.date(2026, 10, 25)))
        with TemporaryDirectory() as tmp:
            root = Path(tmp) / '2026' / 'MHJ45'
            root.mkdir(parents=True)
            pd.DataFrame(dict(station='MHJ45', datetime=['2026-01-01 00:00', '2026-01-01 00:07'], foF2=[5.0, 5.1])).to_csv(root / '2026-01-01.csv', index=False)
            pd.DataFrame(dict(station=[], datetime=[], foF2=[])).to_csv(root / '2026-01-02.csv', index=False)  # empty day
            pd.DataFrame(dict(station='MHJ45', datetime=['2026-01-03 00:00'], foF2=[6.0])).to_csv(root / '2026-01-03.csv', index=False)
            out = pipeline.combine_scaled(tmp, ['MHJ45', 'XX000'], [2026])
            self.assertEqual([p.name for p in out], ['MHJ45_2026_scaled.csv'])
            self.assertEqual(len(pd.read_csv(out[0])), 3)


if __name__ == '__main__':
    unittest.main()
