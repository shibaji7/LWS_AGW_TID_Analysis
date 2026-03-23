import datetime as dt
from pathlib import Path
import pandas as pd

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


def _day_ordinal(day: int) -> str:
    if 11 <= day <= 13:
        return f"{day}th"
    return f"{day}{['th','st','nd','rd','th','th','th','th','th','th'][day % 10]}"


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


def load_fulltimedata(
    date: str | dt.date | dt.datetime,
    beam: int = None,
    mat_dir: str | Path = "data",
    start_time: str | dt.datetime | None = None,
    end_time: str | dt.datetime | None = None,
    elevation_cutoff: float = 0.0,
) -> dict:
    """Load a fulltimedata_<date>_WithinBeam<beam>.mat file (MATLAB v7.3 / HDF5).

    Each row in the file is one TEC measurement. The 12 columns are:

        0  – Time (MATLAB Serial Date Number / datenum)
        1  – Elevation angle (0°–90°)
        2  – Azimuth angle  (-180° to 180°, 0° = geographic north)
        3  – IPP latitude
        4  – IPP longitude
        5  – Absolute vTEC  (not fully implemented, unreliable)
        6  – 10-min high-pass filtered vTEC
        7  – 5–40 min band-pass filtered vTEC
        8  – 60-min high-pass filtered vTEC
        9  – 120-min high-pass filtered vTEC
       10  – Station latitude
       11  – Station longitude

    Parameters
    ----------
    date       : date of the data file
    beam       : beam number used in the filename
    mat_dir    : directory that contains the .mat file
    start_time        : optional lower bound for time filtering
    end_time          : optional upper bound for time filtering
    elevation_cutoff  : discard measurements with elevation < this value (degrees, default 0)

    Returns
    -------
    dict with keys:
        file            – Path to the loaded file
        beam            – beam number
        utt             – raw UTT MATLAB datenum array (one per UTT timestep)
        time            – decoded datetime for each measurement (flat, all timesteps concatenated)
        elevation       – elevation angle per measurement (degrees)
        azimuth         – azimuth angle per measurement (degrees)
        ipp_lat         – IPP latitude per measurement
        ipp_lon         – IPP longitude per measurement
        vtec_abs        – absolute vTEC (unreliable)
        vtec_hp10       – 10-min high-pass filtered vTEC
        vtec_bp5_40     – 5–40 min band-pass filtered vTEC
        vtec_hp60       – 60-min high-pass filtered vTEC
        vtec_hp120      – 120-min high-pass filtered vTEC
        station_lat     – station latitude per measurement
        station_lon     – station longitude per measurement
    """
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError(
            "h5py is required to read fulltimedata files. "
            "Install it in this environment."
        ) from exc

    day = _normalize_date(date)
    mat_dir = Path(mat_dir)

    date_str = f"{_day_ordinal(day.day)}{day:%b}{day.year}"
    if beam is not None:
        mat_file = mat_dir / f"fulltimedata_{date_str}_WithinBeam{beam}.mat"
    else:
        mat_file = mat_dir / f"fulltimedata_05_27_2017-GRE.mat"
    if not mat_file.exists():
        raise FileNotFoundError(f"fulltimedata file not found: {mat_file}")

    print(f"Matfile: {mat_file}")
    start_dt = _normalize_datetime(start_time, day) if start_time is not None else None
    end_dt = _normalize_datetime(end_time, day) if end_time is not None else None

    with h5py.File(mat_file, "r") as h5f:
        utt = np.array(h5f["UTT"]).squeeze()
        times = _decode_time_array(utt)

        fd = h5f["fulltimedata"]
        n_rows = fd.shape[0]

        chunks = {
            "time": [], "elevation": [], "azimuth": [],
            "ipp_lat": [], "ipp_lon": [], "vtec_abs": [],
            "vtec_hp10": [], "vtec_bp5_40": [], "vtec_hp60": [], "vtec_hp120": [],
            "station_lat": [], "station_lon": [],
        }
        _fields = list(chunks.keys())

        for i in range(n_rows):
            t = times[i]

            if start_dt is not None and t < start_dt:
                continue
            if end_dt is not None and t > end_dt:
                continue

            ref = fd[i, 0]
            arr = np.array(h5f[ref])

            # Skip null sentinel rows (last rows are empty uint64 zeros)
            if arr.dtype.kind != "f" or arr.ndim != 2 or arr.shape[0] != 12:
                continue

            n_meas = arr.shape[1]
            chunks["time"].append(_decode_time_array(arr[0, :]))
            chunks["elevation"].append(arr[1, :])
            chunks["azimuth"].append(arr[2, :])
            chunks["ipp_lat"].append(arr[3, :])
            chunks["ipp_lon"].append(arr[4, :])
            chunks["vtec_abs"].append(arr[5, :])
            chunks["vtec_hp10"].append(arr[6, :])
            chunks["vtec_bp5_40"].append(arr[7, :])
            chunks["vtec_hp60"].append(arr[8, :])
            chunks["vtec_hp120"].append(arr[9, :])
            chunks["station_lat"].append(arr[10, :])
            chunks["station_lon"].append(arr[11, :])

    result = {"file": mat_file, "beam": beam, "utt": utt}
    result["time"] = np.concatenate(chunks["time"]).astype(object)
    for field in _fields[1:]:
        result[field] = np.concatenate(chunks[field]).astype(float)
    df = pd.DataFrame({
        'time':        result['time'],
        'elevation':   result['elevation'],
        'azimuth':     result['azimuth'],
        'ipp_lat':     result['ipp_lat'],
        'ipp_lon':     result['ipp_lon'],
        'vtec_abs':    result['vtec_abs'],
        'vtec_hp10':   result['vtec_hp10'],
        'vtec_bp5_40': result['vtec_bp5_40'],
        'vtec_hp60':   result['vtec_hp60'],
        'vtec_hp120':  result['vtec_hp120'],
        'station_lat': result['station_lat'],
        'station_lon': result['station_lon'],
    })
    df['time'] = pd.to_datetime(df['time'])
    if elevation_cutoff > 0:
        print(f"Applying elevation cutoff: {elevation_cutoff}°: len before={len(df)}")
        df = df[df['elevation'] >= elevation_cutoff].reset_index(drop=True)
        print(f"len after={len(df)}")
    return df

if __name__ == "__main__":
    import pandas as pd

    df = load_fulltimedata(
        '2017-05-27', beam=11,
        mat_dir='/home/chakras4/Research/Individual_Studies/LWS_AGW_TID_Analysis/data',
        start_time='16:00', end_time='22:00',
        elevation_cutoff=20.0,
    )

    

    print(df)
    print()
    print(df.dtypes)
    print()
    print(df.describe())