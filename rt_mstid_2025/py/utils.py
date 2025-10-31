import ephem
import numpy as np
from loguru import logger
import pandas as pd
import datetime as dt
import xarray as xr
from tqdm import tqdm

def snap_to_nearest_five(dx: dt.datetime) -> dt.datetime:
    # total minutes from start of day
    minute_total = dx.hour * 60 + dx.minute + dx.second / 60 + dx.microsecond / 60_000_000
    nearest = int(round(minute_total / 5.0) * 5)
    snapped = dx.replace(hour=0, minute=0, second=0, microsecond=0) + dt.timedelta(minutes=nearest)
    return snapped

def get_fov_pyEclipse(
    dates, lats, lons, wl="X", 
    file_ext="_150km_alleof.nc",
    data_folder="/home/chakras4/OneDrive/Chakras4/Projects/Chakraborty.Projects/byProjects/2024 Eclipse Project/Datasets/2024_Apr_Eclipse/mask/"
):
    p = np.nan * np.zeros((len(dates), lats.shape[0], lats.shape[1]))
    for i, date in tqdm(enumerate(dates)):
        date = snap_to_nearest_five(date)
        file = f"{data_folder}{date.strftime('%Y%m%d%H%M%S')}{file_ext}"
        of_data = xr.open_dataset(file)
        sza = of_data["sza"].values
        for j in tqdm(range(lats.shape[0])):
            for k in tqdm(range(lats.shape[1])):
                of = of_data[str(wl)].values
                of[sza > 90] = np.nan
                glats = of_data["glat"].values
                glons = of_data["glon"].values
                ix, jx = np.argmin(np.abs(glats-lats[j,k])), np.argmin(np.abs(glons-lons[j,k]))
                p[i,j,k] = of[ix, jx]
        of_data.close()
    return p

def get_eclipse_contours_pyEclipse(
    dates, lat, lon, wl="X", 
    file_ext="_150km_alleof.nc",
    data_folder="/home/chakras4/OneDrive/Chakras4/Projects/Chakraborty.Projects/byProjects/2024 Eclipse Project/Datasets/2024_Apr_Eclipse/mask/"
):
    Of = []
    for date in dates:
        date = snap_to_nearest_five(date)
        file = f"{data_folder}{date.strftime('%Y%m%d%H%M%S')}{file_ext}"
        of_data = xr.open_dataset(file)
        sza = of_data["sza"].values
        if wl is None:
            o = []
            for wl in [
                    "X", "131", "1600", 
                    "1700", "171", "193",
                    "211", "284", "304", 
                    "335", "94"
                ]:
                of = of_data[str(wl)].values
                of[sza > 90] = np.nan
                lats = of_data["glat"].values
                lons = of_data["glon"].values
                i, j = np.argmin(np.abs(lats-lat)), np.argmin(np.abs(lons-lon))
                o.append(of[i, j])
            Of.append(o)
        else:
            of = of_data[str(wl)].values
            of[sza > 90] = np.nan
            lats = of_data["glat"].values
            lons = of_data["glon"].values
            i, j = np.argmin(np.abs(lats-lat)), np.argmin(np.abs(lons-lon))
            Of.append(of[i, j])
        of_data.close()
    Of = np.array(Of)
    return Of

def get_eclipse_of_pyEclipse(
    dates: list = [dt.datetime(2024, 4, 8, 16)], wl="X",
    file_ext="_150km_alleof.nc",
    data_folder="/home/chakras4/OneDrive/Chakras4/Projects/Chakraborty.Projects/byProjects/2024 Eclipse Project/Datasets/2024_Apr_Eclipse/mask/"
):
    Of = []
    for date in dates:
        file = f"{data_folder}{date.strftime('%Y%m%d%H%M%S')}{file_ext}"
        of_data = xr.open_dataset(file)
        sza = of_data["sza"].values
        lat = of_data["glat"].values
        lon = of_data["glon"].values
        of = of_data[str(wl)].values
        of[sza > 90] = np.nan
        Of.append(of)
        of_data.close()
    return lat, lon, np.array(Of)



class Eclipse(object):
    def __init__(self):
        return

    def calculate_w2naf_shadow(self, d, lat, lon, alt=300.0):
        obsc = eclipse_calc.calculate_obscuration(d, lat, lon, alt)
        return obsc

    def intersection(slef, r0, r1, d, n_s=100):
        A1 = np.zeros([n_s, n_s])
        A2 = np.zeros([n_s, n_s])
        I = np.zeros([n_s, n_s])
        x = np.linspace(-2.0 * r0, 2.0 * r0, num=n_s)
        y = np.linspace(-2.0 * r0, 2.0 * r0, num=n_s)
        xx, yy = np.meshgrid(x, y)
        A1[np.sqrt((xx + d) ** 2.0 + yy**2.0) < r0] = 1.0
        n_sun = np.sum(A1)
        A2[np.sqrt(xx**2.0 + yy**2.0) < r1] = 1.0
        S = A1 + A2
        I[S > 1] = 1.0
        eclipse = np.sum(I) / n_sun
        return eclipse

    def create_eclipse_shadow(self, d, lat, lon, alt):
        obs = ephem.Observer()
        t0 = ephem.date(
            (
                d.year,
                d.month,
                d.day,
                d.hour,
                d.minute,
                d.second,
            )
        )
        obs.lon, obs.lat = "%1.2f" % (lon), "%1.2f" % (lat)  # ESR
        obs.elevation = alt
        obs.date = t0
        sun, moon = ephem.Sun(), ephem.Moon()
        sun.compute(obs)
        moon.compute(obs)
        r_sun = (sun.size / 2.0) / 3600.0
        r_moon = (moon.size / 2.0) / 3600.0
        s = np.degrees(ephem.separation((sun.az, sun.alt), (moon.az, moon.alt)))
        percent_eclipse = 0.0

        if s < (r_moon + r_sun):
            if s < 1e-3:
                percent_eclipse = 1.0
            else:
                percent_eclipse = self.intersection(r_moon, r_sun, s, n_s=100)
        if np.degrees(sun.alt) <= r_sun:
            if np.degrees(sun.alt) <= -r_sun:
                percent_eclipse = 2
            else:
                percent_eclipse = 1.0 - (
                    (np.degrees(sun.alt) + r_sun) / (2.0 * r_sun)
                ) * (1.0 - percent_eclipse)
        return percent_eclipse


def create_eclipse_path_local(dates, lat, lon, alt=300, limit=None):
    e = Eclipse()
    p = np.nan * np.zeros((len(dates)))
    for i, d in enumerate(tqdm(dates)):
        shadow = e.create_eclipse_shadow(d, lat, lon, alt)
        if limit is not None:
            shadow[shadow < limit] = np.float64(0.)
        p[i] = shadow
    return p


def get_fov_eclipse(
    dates, lats, lons, alt=300, limit=None
):
    from tqdm import tqdm
    e = Eclipse()
    p = np.nan * np.zeros((len(dates), lats.shape[0], lats.shape[1]))
    for i, date in tqdm(enumerate(dates)):
        for j in tqdm(range(lats.shape[0])):
            for k in tqdm(range(lats.shape[1])):
                shadow = e.create_eclipse_shadow(date, lats[j,k], lons[j,k], alt)
                if limit is not None:
                    shadow = shadow if shadow >= limit else 0.
                p[i, j, k] = shadow
    return p


def extract_all_ionosonde_datasets_from_inogram_image(folder, date, code):
    import glob
    from image_parser import IonogramTableExtractor
    files = glob.glob(f"{folder}/*.png")
    logger.info(f"Found {len(files)} files in {folder}")
    files.sort()
    records = []
    for i, f in enumerate(files):
        logger.info(f"Extracting {f}")
        hours, minutes, seconds = (
            int(f.split("/")[-1].split("_")[-1].split(".")[0][:2]),
            int(f.split("/")[-1].split("_")[-1].split(".")[0][2:4]),
            int(f.split("/")[-1].split("_")[-1].split(".")[0][4:]),
        )
        record = IonogramTableExtractor.extract_ionogram_table(f, verbose=False)
        record["date"] = date.replace(
            hour=hours, minute=minutes, second=seconds
        )
        record["image_path"] = f
        records.append(record)
    records = pd.concat(records, ignore_index=True)

    records.to_csv(f"data/{date.year}/{code}.csv", index=False, header=True, float_format="%g")
    return

def interpolate_missing_values(
        df, xparam, yparam, intv=60., 
        cutoff_norm = 0.0006, order=2
):
    import datetime as dt
    from scipy import signal

    xvalues, yvalues = df[xparam].tolist(), np.array(df[yparam])
    dx = np.array([(x-xvalues[0]).total_seconds() for x in xvalues])
    dx = dx[~np.isnan(yvalues)]
    yvalues = yvalues[~np.isnan(yvalues)]
    n_minutes = (xvalues[-1] - xvalues[0]).total_seconds() / intv
    newdx = np.arange(0, n_minutes*intv, intv)
    yvalues_new = np.interp(newdx, dx, yvalues)
    dates = [xvalues[0] + dt.timedelta(seconds=x) for x in newdx]

    fs = 1/intv
    nyquist = 0.5 * fs
    cutoff_norm = cutoff_norm / nyquist
    b, a = signal.butter(order, cutoff_norm, btype='low', analog=False)
    yvalues_new = signal.lfilter(b, a, yvalues_new)

    return dates, yvalues_new

def plasma_freq_mhz(ne_per_m3: float) -> float:
    return 8.98e-6 * np.sqrt(ne_per_m3)