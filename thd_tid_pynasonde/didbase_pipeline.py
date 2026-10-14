"""Download CONUS DIDBASE observations for April-September 2024-2026.

Outputs (relative to this script), default settings from ``config.toml``:

    Data/didbase_giro_scaled/<year>/<station>/<YYYY-MM-DD>.csv          (scaled params)
    Data/didbase_giro_scaled/<year>/<station>/<YYYY-MM-DD>.json         (completion record)
    Data/didbase_giro_scaled/<year>/<station>/raw/<YYYY-MM-DD>/*.RSF    (--raw only)
    Data/didbase_giro_scaled/<year>/<station>/raw/<YYYY-MM-DD>.json     (--raw only)
    Data/didbase_giro_scaled/<year>/<station>/sky/<YYYY-MM-DD>/*.SKY    (--sky only)
    Data/didbase_giro_scaled/<year>/<station>/sky/<YYYY-MM-DD>.json     (--sky only)

Scaled downloads request every catalogued parameter, quality flag and
provenance field. The connector selects one scaler per measurement,
preferring human review; these files do not archive all alternative scaler
revisions. Raw downloads (``--raw``) fetch the underlying ionogram files
(RSF, falling back to MMM) that ``plot_viterbi.py`` needs to draw traces;
DIDBASE's scaled-parameter query alone cannot supply these. All dates are
UTC. Historical downloads resume using checksummed JSON completion records.
Empty days are recorded too. --overwrite refreshes late arrivals/revisions.
Today's partial data are always refreshed; future dates are never queried.

Geography is discovered from the Digisonde instrument catalogue and
intersected with DIDBASE's live station list. --station-catalog supplies an
alternative CSV with URSI, Country, Lat, Long columns. Unknown stations are
reported for review.

``--sky`` downloads raw DDM skymap files (SKY) from GIRO DriftBase instead (via
``pynasonde.digisonde.driftbase``; env ``DRIFTBASE_URL/USER/PASSWORD/ROLE``).
Combine with ``--stations`` and ``--days-csv`` (e.g. ``Data/analysis/mstid_days.csv``,
using its ``station``/``date`` columns and, if present, only rows with
``enhanced == True``) to fetch skymaps just for event days.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import time
from itertools import product
from urllib.request import Request, urlopen

import sys

import pandas as pd
from loguru import logger
from tqdm import tqdm

from utils import load_config

# loguru's default handler level is DEBUG, so logger.debug(...) calls (used
# for the per-station-day success lines that would otherwise tear the tqdm
# progress bar apart) print by default just like INFO. Raise the default
# handler to INFO; set DIDBASE_LOG_LEVEL=DEBUG to see per-day detail again.
logger.remove()
logger.add(sys.stderr, level=os.environ.get('DIDBASE_LOG_LEVEL', 'INFO'))

_CFG = load_config()
_DOWNLOAD_CFG = _CFG["download"]
_STATIONS_CFG = _CFG["stations"]

DATA_DIR = Path(__file__).resolve().parent / _DOWNLOAD_CFG["output_dir"]
STATION_URL = 'https://www.digisonde.com/stationlist.php'
DEFAULT_YEARS = tuple(_DOWNLOAD_CFG["years"])
DEFAULT_MONTH_START = _DOWNLOAD_CFG["month_start"]
DEFAULT_MONTH_END = _DOWNLOAD_CFG["month_end"]
DEFAULT_RETRIES = _DOWNLOAD_CFG["retries"]
RAW_FORMATS = tuple(_DOWNLOAD_CFG["raw_formats"])


class _TableRows(HTMLParser):
    """Read the catalogue without optional HTML parsing dependencies."""
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.cell = [], [], None

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.row = []
        elif tag in ('td', 'th'):
            self.cell = []

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.cell is not None:
            self.row.append(' '.join(''.join(self.cell).split()))
            self.cell = None
        elif tag == 'tr' and self.row:
            self.rows.append(self.row)


def parse_station_html(html):
    # The live site embeds the complete table as JSON for client-side rendering.
    match = re.search(r'\bvar\s+jsonData\s*=\s*', html)
    if match:
        rows, _ = json.JSONDecoder().raw_decode(html[match.end():])
        if rows and {'URSI', 'Country', 'Lat', 'Long'}.issubset(rows[0]):
            return pd.DataFrame(rows[1:], columns=rows[0])
    parser = _TableRows()
    parser.feed(html)
    for index, row in enumerate(parser.rows):
        if {'URSI', 'Country', 'Lat', 'Long'}.issubset(row):
            return pd.DataFrame([r for r in parser.rows[index + 1:] if len(r) == len(row)], columns=row)
    raise ValueError('Station catalogue has no URSI/Country/Lat/Long table')


def load_station_catalog(path=None):
    if path is not None:
        return pd.read_csv(path)
    request = Request(STATION_URL, headers={'User-Agent': 'DIDBASE-CONUS-downloader/1.0'})
    with urlopen(request, timeout=60) as response:
        catalog = parse_station_html(response.read().decode('utf-8'))
    # Stations present in DIDBASE but absent from the public instrument
    # catalogue are supplemented from config.toml's [[stations.extra]].
    extra = [row for row in _STATIONS_CFG.get("extra", []) if row["ursi"] not in set(catalog.URSI)]
    if extra:
        catalog = pd.concat([catalog, pd.DataFrame([
            dict(URSI=row["ursi"], Country=row["country"], Lat=row["lat"],
                 Long=row["long"], Location=row["location"])
            for row in extra
        ])], ignore_index=True)
    return catalog


def select_conus_stations(available, catalog):
    """Country plus coordinates excludes Canada/Mexico, AK/HI and territories."""
    catalog = catalog.copy()
    catalog.columns = catalog.columns.str.strip().str.upper()
    if not {'URSI', 'COUNTRY', 'LAT', 'LONG'}.issubset(catalog.columns):
        raise ValueError('Station catalogue requires URSI, Country, Lat, Long')
    catalog['URSI'] = catalog.URSI.fillna('').astype(str).str.strip().str.upper()
    catalog['LAT'] = pd.to_numeric(catalog.LAT, errors='coerce')
    catalog['LONG'] = (pd.to_numeric(catalog.LONG, errors='coerce') + 180) % 360 - 180
    available = available.copy()
    available['URSI'] = available.URSI.astype(str).str.strip().str.upper()
    unknown = sorted(set(available.URSI) - set(catalog.URSI))
    if unknown:
        logger.warning('Stations without catalogue geography: {}. Review --station-catalog for omissions.', ', '.join(unknown))
    country = catalog.COUNTRY.astype(str).str.strip().str.upper()
    selected = catalog[
        country.isin({c.upper() for c in _STATIONS_CFG["country_names"]})
        & catalog.LAT.between(_STATIONS_CFG["lat_min"], _STATIONS_CFG["lat_max"])
        & catalog.LONG.between(_STATIONS_CFG["lon_min"], _STATIONS_CFG["lon_max"])
        & catalog.URSI.str.fullmatch(r'[A-Z0-9]{5}')
    ].drop_duplicates('URSI')
    return available.merge(selected[['URSI', 'COUNTRY', 'LAT', 'LONG']], on='URSI',
                           how='inner', suffixes=('_DIDBASE', '')).sort_values('URSI').reset_index(drop=True)


def requested_dates(years, today, month_start=DEFAULT_MONTH_START, month_end=DEFAULT_MONTH_END):
    for year in sorted(set(years)):
        day = dt.date(year, month_start, 1)
        last = (dt.date(year, month_end + 1, 1) - dt.timedelta(days=1)) if month_end < 12 else dt.date(year, 12, 31)
        while day <= min(last, today):
            yield day
            day += dt.timedelta(days=1)


def _atomic_csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    frame.to_csv(temp, index=False)
    temp.replace(path)


def _atomic_json(payload, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(payload, indent=2) + '\n')
    temp.replace(path)


def download_day(db, station, day, params, output_dir, *, now, overwrite=False):
    path = output_dir / str(day.year) / station / f'{day.isoformat()}.csv'
    status_path = path.with_suffix('.json')
    if not overwrite and day < now.date() and path.is_file() and status_path.is_file():
        try:
            status = json.loads(status_path.read_text())
            if (status.get('complete') and status.get('parameters') == params
                    and status.get('station') == station and status.get('date') == day.isoformat()
                    and status.get('sha256') == hashlib.sha256(path.read_bytes()).hexdigest()):
                return 'cached'
        except (OSError, ValueError):
            pass
    start = dt.datetime.combine(day, dt.time.min)
    end = min(dt.datetime.combine(day, dt.time.max), now)
    frame = db.fetch_scaled_params(station=station, start=start, end=end, params=params,
                                   prefer_human=True, include_quality=True, include_provenance=True).copy()
    if frame.empty:
        frame = pd.DataFrame(columns=['datetime', 'scaler_id'])
    else:
        times = pd.to_datetime(frame.datetime, utc=True, errors='raise')
        if times.isna().any() or not times.dt.date.eq(day).all():
            raise ValueError(f'{station} {day}: invalid or out-of-day timestamps')
        frame['datetime'] = times.dt.tz_localize(None)
        frame = frame.sort_values('datetime').reset_index(drop=True)
    columns = [c for p in params for c in (p, f'{p}_qletter', f'{p}_dletter')]
    columns += ['scaling_quality', 'scaling_confidence', 'equipment_id', 'program_id']
    frame = frame.reindex(columns=list(frame.columns) + [c for c in columns if c not in frame])
    frame.insert(0, 'station', station)
    _atomic_csv(frame, path)
    result = 'empty' if frame.empty else 'saved'
    _atomic_json(dict(station=station, date=day.isoformat(), status=result,
        rows=len(frame), parameters=params, complete=day < now.date(), queried_through_utc=end.isoformat(),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest()), status_path)
    logger.debug('{} {}: {} ({} rows)', station, day, result, len(frame))
    return result


def download_day_raw(db, station, day, output_dir, *, now, overwrite=False, formats=RAW_FORMATS):
    """Fetch raw ionogram files for one station-day; feeds ``plot_viterbi.py``."""
    day_dir = output_dir / str(day.year) / station / 'raw' / day.isoformat()
    status_path = day_dir.with_suffix('.json')
    if not overwrite and day < now.date() and status_path.is_file():
        try:
            status = json.loads(status_path.read_text())
            if (status.get('complete') and status.get('station') == station
                    and status.get('date') == day.isoformat()):
                return 'cached'
        except (OSError, ValueError):
            pass
    start = dt.datetime.combine(day, dt.time.min)
    end = min(dt.datetime.combine(day, dt.time.max), now)
    day_dir.mkdir(parents=True, exist_ok=True)
    used_format, paths = None, []
    for fmt in formats:
        paths = db.fetch_files(station=station, start=start, end=end, fmt=fmt, out_dir=day_dir)
        if paths:
            used_format = fmt
            break
    result = 'empty' if not paths else 'saved'
    _atomic_json(dict(station=station, date=day.isoformat(), status=result, files=len(paths),
        format=used_format, complete=day < now.date(), queried_through_utc=end.isoformat()), status_path)
    logger.debug('{} {} (raw): {} ({} {} files)', station, day, result, len(paths), used_format or '-')
    return result


def download_day_sky(db, station, day, output_dir, *, now, overwrite=False):
    """Fetch raw skymap (SKY) files for one station-day from DriftBase."""
    day_dir = output_dir / str(day.year) / station / 'sky' / day.isoformat()
    status_path = day_dir.with_suffix('.json')
    if not overwrite and day < now.date() and status_path.is_file():
        try:
            status = json.loads(status_path.read_text())
            if (status.get('complete') and status.get('station') == station
                    and status.get('date') == day.isoformat()):
                return 'cached'
        except (OSError, ValueError):
            pass
    start = dt.datetime.combine(day, dt.time.min)
    end = min(dt.datetime.combine(day, dt.time.max), now)
    day_dir.mkdir(parents=True, exist_ok=True)
    paths = db.fetch_sky_files(station=station, start=start, end=end, out_dir=day_dir)
    result = 'empty' if not paths else 'saved'
    _atomic_json(dict(station=station, date=day.isoformat(), status=result, files=len(paths),
        complete=day < now.date(), queried_through_utc=end.isoformat()), status_path)
    logger.debug('{} {} (sky): {} ({} files)', station, day, result, len(paths))
    return result


def event_station_days(csv_path, stations=None, years=None, month_start=DEFAULT_MONTH_START, month_end=DEFAULT_MONTH_END):
    """(station, date) pairs from a CSV with ``date`` (and optionally ``station``/``enhanced``) columns."""
    frame = pd.read_csv(csv_path, parse_dates=['date'])
    if 'enhanced' in frame:
        frame = frame[frame.enhanced.astype(bool)]
    if 'station' in frame and stations:
        frame = frame[frame.station.isin(stations)]
    if years:
        frame = frame[frame.date.dt.year.isin(years)]
    frame = frame[frame.date.dt.month.between(month_start, month_end)]
    cols = ['station', 'date'] if 'station' in frame else ['date']
    pairs = frame[cols].drop_duplicates().sort_values(cols)
    if 'station' not in pairs:
        if not stations:
            raise ValueError('days CSV has no station column; pass --stations')
        return [(s, d.date()) for s in stations for d in pairs.date]
    return [(r.station, r.date.date()) for r in pairs.itertuples()]


def run_sky(years=None, stations=None, output_dir=DATA_DIR, days_csv=None,
            overwrite=False, retries=DEFAULT_RETRIES, dates=None):
    from pynasonde.digisonde.driftbase import DriftBaseConnector

    if not stations and not days_csv and not dates:
        raise ValueError('--sky needs --stations and/or --days-csv/--dates (skymaps are large)')
    if dates and not stations and not days_csv:
        raise ValueError('--dates needs --stations')
    now = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None)
    stations = [s.strip().upper() for s in stations] if stations else None
    if dates:
        days = sorted({dt.date.fromisoformat(str(d)) for d in dates})
        pairs = [(s, d) for s in stations for d in days]  # explicit dates: no season/year filtering
    elif days_csv:
        pairs = event_station_days(days_csv, stations, years)  # years=None: use the CSV dates as given
    else:
        pairs = [(s, d) for s in stations for d in requested_dates(years or DEFAULT_YEARS, now.date())]
    pairs = [(s, d) for s, d in pairs if d <= now.date()]
    if not pairs:
        logger.info('No station-days to download.')
        return
    output_dir = Path(output_dir)
    logger.info('Skymap download: {} station-days ({})', len(pairs), ', '.join(sorted({s for s, _ in pairs})))
    with DriftBaseConnector(
        os.environ.get('DRIFTBASE_URL', 'jdbc:firebirdsql://129.63.134.45/drift'),
        os.environ.get('DRIFTBASE_USER', os.environ.get('DIDBASE_USER', 'ERAU')),
        os.environ.get('DRIFTBASE_PASSWORD', os.environ.get('DIDBASE_PASSWORD', 'eraucar')),
        os.environ.get('DRIFTBASE_ROLE', os.environ.get('DIDBASE_ROLE', 'COMMON')),
    ) as db:
        counts, failures = _download_with_progress(
            None, None, retries, desc=f'{output_dir.name} (sky)', pairs=pairs,
            fetch_one=lambda station, day: download_day_sky(db, station, day, output_dir, now=now, overwrite=overwrite),
            label=lambda station, day: f'{station} {day} (sky)',
        )
    _atomic_csv(pd.DataFrame(failures, columns=['station', 'date', 'error']), output_dir / 'sky_failures.csv')
    logger.info('Skymap download summary: {}', counts)
    if failures:
        raise RuntimeError(f"{len(failures)} skymap station-days failed; see {output_dir / 'sky_failures.csv'}")


def _select_stations(db, station_catalog, stations):
    catalog = load_station_catalog(station_catalog)
    available = db.stations()
    selected = select_conus_stations(available, catalog)
    if stations:
        wanted = {s.strip().upper() for s in stations}
        missing = wanted - set(selected.URSI)
        if missing:
            raise ValueError(f'Stations not found in CONUS DIDBASE selection: {sorted(missing)}')
        selected = selected[selected.URSI.isin(wanted)]
    if selected.empty:
        raise RuntimeError('No CONUS stations found; check the station catalogue')
    return catalog, available, selected


def run(years=None, stations=None, output_dir=DATA_DIR,
        station_catalog=None, overwrite=False, retries=DEFAULT_RETRIES,
        list_stations=False, raw=False, sky=False, days_csv=None, dates=None):
    if sky:
        return run_sky(years=years, stations=stations, output_dir=output_dir,
                       days_csv=days_csv, overwrite=overwrite, retries=retries, dates=dates)
    years = years or DEFAULT_YEARS
    from pynasonde.digisonde.didbase import DidBaseConnector

    now = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None)
    days = list(requested_dates(years, now.date()))
    if not days:
        logger.info('No requested dates have occurred yet.')
        return
    output_dir = Path(output_dir)
    with DidBaseConnector(
        os.environ.get('DIDBASE_URL', 'jdbc:firebirdsql://didbase.giro.uml.edu/didb'),
        os.environ.get('DIDBASE_USER', 'ERAU'), os.environ.get('DIDBASE_PASSWORD', 'eraucar'),
        os.environ.get('DIDBASE_ROLE', 'COMMON'),
    ) as db:
        catalog, available, selected = _select_stations(db, station_catalog, stations)
        logger.info('Selected {} stations: {}', len(selected), ', '.join(selected.URSI))
        if list_stations:
            print(selected.to_string(index=False))
            return
        _atomic_csv(catalog, output_dir / 'station_catalog.csv')
        _atomic_csv(available, output_dir / 'didbase_stations.csv')
        _atomic_csv(selected, output_dir / 'conus_stations.csv')

        if raw:
            counts, failures = _download_with_progress(
                selected.URSI, days, retries, desc=f'{output_dir.name} (raw)',
                fetch_one=lambda station, day: download_day_raw(db, station, day, output_dir, now=now, overwrite=overwrite),
                label=lambda station, day: f'{station} {day} (raw)',
            )
            _atomic_csv(pd.DataFrame(failures, columns=['station', 'date', 'error']), output_dir / 'raw_failures.csv')
            logger.info('Raw download summary: {}', counts)
            if failures:
                raise RuntimeError(f"{len(failures)} raw station-days failed; see {output_dir / 'raw_failures.csv'}")
            return

        param_catalog = db.scaled_params_catalog()
        params = sorted(set(param_catalog.NAME.dropna().astype(str).str.strip()) - {''})
        if not params or any(len(p) > 15 for p in params):
            raise ValueError('Empty or invalid DIDBASE parameter catalogue')
        _atomic_csv(param_catalog, output_dir / 'parameter_catalog.csv')
        counts, failures = _download_with_progress(
            selected.URSI, days, retries, desc=output_dir.name,
            fetch_one=lambda station, day: download_day(db, station, day, params, output_dir, now=now, overwrite=overwrite),
            label=lambda station, day: f'{station} {day}',
        )
        _atomic_csv(pd.DataFrame(failures, columns=['station', 'date', 'error']), output_dir / 'failures.csv')
        logger.info('Download summary: {}', counts)
        if failures:
            raise RuntimeError(f"{len(failures)} station-days failed; see {output_dir / 'failures.csv'}")


def _download_with_progress(stations, days, retries, *, desc, fetch_one, label, pairs=None):
    """Run (station, day) downloads behind a single updating tqdm bar.

    ``fetch_one(station, day)`` should call the actual downloader and return
    its result string ('saved'/'empty'/'cached'); this handles retries,
    per-day progress, and the running saved/empty/cached/failed postfix.
    """
    pairs = list(pairs) if pairs is not None else list(product(list(stations), days))
    counts = dict(saved=0, empty=0, cached=0, failed=0)
    failures = []
    # pynasonde's own DidBaseConnector logs an INFO line per station-day
    # query (plus a full parameter-list dump); left enabled it tears the
    # single-line tqdm bar apart on every iteration. Silence it for the
    # duration of the loop only -- our own warnings/errors still show.
    logger.disable('pynasonde')
    try:
        with tqdm(total=len(pairs), desc=desc, unit='day',
                  dynamic_ncols=True, mininterval=0.2) as bar:
            for station, day in pairs:
                for attempt in range(retries):
                    try:
                        counts[fetch_one(station, day)] += 1
                        break
                    except Exception as exc:
                        logger.warning('{} attempt {}/{}: {}', label(station, day), attempt + 1, retries, exc)
                        if attempt + 1 == retries:
                            counts['failed'] += 1
                            failures.append(dict(station=station, date=str(day), error=str(exc)))
                        else:
                            time.sleep(min(2 ** attempt, 10))
                bar.set_postfix(counts, refresh=False)
                bar.update(1)
    finally:
        logger.enable('pynasonde')
    return counts, failures


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--years', type=int, nargs='+', default=None,
                         help='Years to download (default: config.toml [download] years). With --sky --days-csv the CSV '
                              'dates are used as given unless --years is passed.')
    parser.add_argument('--stations', nargs='+', help='Optional subset of discovered CONUS URSI codes')
    parser.add_argument('--output-dir', type=Path, default=DATA_DIR)
    parser.add_argument('--station-catalog', type=Path, help='Alternate URSI/Country/Lat/Long CSV')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--list-stations', action='store_true', help='Discover stations without downloading observations')
    parser.add_argument('--raw', action='store_true',
                         help='Fetch raw ionogram files (RSF/MMM) instead of scaled parameters; feeds plot_viterbi.py')
    parser.add_argument('--sky', action='store_true',
                         help='Fetch raw DDM skymap (SKY) files from DriftBase; use with --stations/--days-csv')
    parser.add_argument('--dates', nargs='+', metavar='YYYY-MM-DD',
                         help='With --sky: download exactly these dates for --stations (e.g. one event day)')
    parser.add_argument('--days-csv', type=Path,
                         help='CSV with date (and optional station/enhanced) columns restricting --sky downloads, '
                              'e.g. Data/analysis/mstid_days.csv (enhanced rows only)')
    parser.add_argument('--retries', type=int, default=DEFAULT_RETRIES)
    args = parser.parse_args()
    if args.retries < 1:
        parser.error('--retries must be at least 1')
    return args


if __name__ == '__main__':
    run(**vars(parse_args()))
