import datetime as dt
from pathlib import Path

import numpy as np


def _normalize_date(date_value: str | dt.date | dt.datetime) -> dt.date:
    if isinstance(date_value, dt.datetime):
        return date_value.date()
    if isinstance(date_value, dt.date):
        return date_value
    if isinstance(date_value, str):
        for fmt in ("%Y-%m-%d", "%Y_%m_%d", "%Y%m%d"):
            try:
                return dt.datetime.strptime(date_value, fmt).date()
            except ValueError:
                continue
    raise ValueError(
        "Invalid date format. Use datetime/date or string in YYYY-MM-DD, YYYY_MM_DD, or YYYYMMDD."
    )


def _normalize_datetime(
    time_value: str | dt.datetime | None, day: dt.date
) -> dt.datetime | None:
    if time_value is None:
        return None
    if isinstance(time_value, dt.datetime):
        return time_value
    if isinstance(time_value, str):
        dt_formats = (
            "%Y-%m-%d %H:%M",
            "%Y-%m-%d %H:%M:%S",
            "%Y-%m-%dT%H:%M",
            "%Y-%m-%dT%H:%M:%S",
        )
        for fmt in dt_formats:
            try:
                return dt.datetime.strptime(time_value, fmt)
            except ValueError:
                continue
        t_formats = ("%H:%M", "%H:%M:%S")
        for fmt in t_formats:
            try:
                parsed_t = dt.datetime.strptime(time_value, fmt).time()
                return dt.datetime.combine(day, parsed_t)
            except ValueError:
                continue
    raise ValueError(
        "Invalid specific_time format. Use datetime or string in "
        "YYYY-MM-DD HH:MM[:SS], YYYY-MM-DDTHH:MM[:SS], or HH:MM[:SS]."
    )


def _matlab_datenum_to_datetime(dn: float) -> dt.datetime:
    day = dt.datetime.fromordinal(int(dn))
    frac = dt.timedelta(days=float(dn) % 1)
    return day + frac - dt.timedelta(days=366)


def _round_datetime_to_minute(value: dt.datetime) -> dt.datetime:
    rounded = value.replace(second=0, microsecond=0)
    if value.second >= 30 or value.microsecond >= 500000:
        rounded += dt.timedelta(minutes=1)
    return rounded


def _decode_time_array(time_arr):
    arr = np.asarray(time_arr).squeeze()
    decoded = []
    for val in arr:
        x = float(val)
        if 7.0e5 <= x <= 9.0e5:
            decoded.append(_round_datetime_to_minute(_matlab_datenum_to_datetime(x)))
        elif x >= 1.0e12:
            decoded.append(
                _round_datetime_to_minute(dt.datetime.utcfromtimestamp(x / 1000.0))
            )
        else:
            decoded.append(_round_datetime_to_minute(dt.datetime.utcfromtimestamp(x)))
    return np.array(decoded, dtype=object)


def _load_mat(mat_path: Path):
    try:
        from scipy.io import loadmat

        data = loadmat(mat_path)
        return {k: v for k, v in data.items() if not k.startswith("__")}
    except NotImplementedError:
        pass

    try:
        import h5py
    except Exception as exc:
        raise RuntimeError(
            "MATLAB v7.3 reader unavailable. Install a compatible h5py build in this environment."
        ) from exc

    with h5py.File(mat_path, "r") as h5f:
        data = {k: np.array(h5f[k]) for k in h5f.keys()}
    return data


def _subset_by_time_axis(arr: np.ndarray | None, mask: np.ndarray, time_len: int):
    if arr is None:
        return None
    arr = np.asarray(arr)
    if arr.ndim == 0:
        return arr
    for axis, size in enumerate(arr.shape):
        if size-1 == time_len:
            return np.compress(mask, arr, axis=axis)
    return arr


def _slice_at_time_index(arr: np.ndarray | None, time_index: int, time_len: int):
    if arr is None:
        return None
    arr = np.asarray(arr)
    if arr.ndim == 0:
        return arr
    for axis, size in enumerate(arr.shape):
        if size-1 == time_len:
            return np.take(arr, indices=time_index, axis=axis)
    return arr


def load_nexrad_data(
    date: dt.datetime | dt.date,
    mat_dir: str | Path = "../data",
    downsample_step: int = 4,
):
    day = _normalize_date(date)

    mat_dir = Path(mat_dir)
    mat_file = mat_dir / f"nexrad_{day:%Y_%m_%d}.mat"
    if not mat_file.exists():
        raise FileNotFoundError(f"NEXRAD file not found: {mat_file}")

    data = _load_mat(mat_file)
    required = {"lat", "lon", "timefull"}
    missing = required.difference(data.keys())
    if missing:
        raise KeyError(f"Missing variables in {mat_file.name}: {sorted(missing)}")

    lat = np.asarray(data["lat"]).squeeze()
    lon = np.asarray(data["lon"]).squeeze()
    timefull = np.asarray(data["timefull"]).squeeze()
    time = _decode_time_array(timefull)
    precipfull = data.get("precipfull")
    precip_raw = data.get("precip")

    # MATLAB equivalent: nexradlat = lat(1:4:end); nexradlon = lon(1:4:end);
    lat_ds = lat[::downsample_step]
    lon_ds = lon[::downsample_step]

    return {
        "file": mat_file,
        "lat": lat_ds,
        "lon": lon_ds,
        "timefull": timefull,
        "time": time,
        "precip": precip_raw,
        "precipfull": precipfull,
        "raw": data,
    }


def get_nexrad_data_by_date(
    date: str | dt.date | dt.datetime,
    mat_dir: str | Path = "../data",
    downsample_step: int = 4,
    start_time: dt.datetime | None = None,
    end_time: dt.datetime | None = None,
):
    day = _normalize_date(date)
    data = load_nexrad_data(
        date=day,
        mat_dir=mat_dir,
        downsample_step=downsample_step,
    )
    time = data["time"]
    mask = np.array([t == date for t in time], dtype=bool)
    if start_time is not None:
        mask &= np.array([t >= start_time for t in time], dtype=bool)
    if end_time is not None:
        mask &= np.array([t <= end_time for t in time], dtype=bool)

    time_indices = np.flatnonzero(mask)
    print(f"Selected date {date} at indices {time_indices}")
    precip_window = _subset_by_time_axis(data["precipfull"], mask, len(time))
    precip_window[precip_window<=0] = np.nan

    return {
        "file": data["file"],
        "time": time[mask],
        "time_index": time_indices,
        "lat": data["lat"],
        "lon": data["lon"],
        "precip": precip_window,
    }
