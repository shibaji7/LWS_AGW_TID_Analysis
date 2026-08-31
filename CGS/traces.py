"""
3-D NVIS raytrace experiment with synthetic TID/bump structure.

This is a forward-modeling diagnostic for comparing modeled turning-point
locations with WI937 VIPIR XL/YL echo-location skymaps.
"""

import argparse
import datetime as dt
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.extend([
    str(Path(__file__).resolve().parents[1]),
    str(Path(__file__).resolve().parents[2]),
])

TRACE_ROOT = Path("/home/chakras4/Research/CodeBase/trace")
if str(TRACE_ROOT) not in sys.path:
    sys.path.insert(0, str(TRACE_ROOT))

from hfpytrace.model.rt3d import RT3D, RT3DProfile

# from apep.jgr.style import apply_jgr_style


OUTPUT_DIR = Path("manuscript_figures/jgr")
STATION = "WI937"
STATION_LAT = 37.94
STATION_LON = -75.58
EVENT_TIME = dt.datetime(2024, 4, 8, 19, 20, 0)

REGIONS = {
    "e": {"title": "E-region", "freqs": np.arange(1.5, 3.01, 0.5)},
    "f": {"title": "F-region", "freqs": np.concatenate([np.arange(3.0, 6.01, 1.0), np.arange(8.0, 14.01, 2.0)])},
}


def plasma_freq_mhz(ne_m3):
    """Electron density [m^-3] to plasma frequency [MHz]."""
    return 8.98e-6 * np.sqrt(np.asarray(ne_m3, dtype=float))


def ne_from_fp_mhz(fp_mhz):
    """Plasma frequency [MHz] to electron density [m^-3]."""
    return (np.asarray(fp_mhz, dtype=float) / 8.98e-6) ** 2


def local_axes_to_latlon(x_km, y_km, lat0=STATION_LAT, lon0=STATION_LON):
    km_per_deg_lat = 111.32
    km_per_deg_lon = km_per_deg_lat * np.cos(np.deg2rad(lat0))
    lats = lat0 + np.asarray(y_km) / km_per_deg_lat
    lons = lon0 + np.asarray(x_km) / km_per_deg_lon
    return lats, lons


def latlon_to_local_axes(lats, lons, lat0=STATION_LAT, lon0=STATION_LON):
    km_per_deg_lat = 111.32
    km_per_deg_lon = km_per_deg_lat * np.cos(np.deg2rad(lat0))
    y_km = (np.asarray(lats, dtype=float) - lat0) * km_per_deg_lat
    x_km = (np.asarray(lons, dtype=float) - lon0) * km_per_deg_lon
    return x_km, y_km


def build_synthetic_profile(dense=False):
    if dense:
        x_km = np.linspace(-180.0, 180.0, 361)
        y_km = np.linspace(-180.0, 180.0, 361)
        alts = np.linspace(0.0, 500.0, 151)
    else:
        x_km = np.linspace(-180.0, 180.0, 61)
        y_km = np.linspace(-180.0, 180.0, 61)
        alts = np.linspace(0.0, 500.0, 126)
    lats, lons = local_axes_to_latlon(x_km, y_km)

    yy, xx, zz = np.meshgrid(y_km, x_km, alts, indexing="ij")

    layers = [
        {"h": 112.0, "scale": 12.0, "fp": 3.0},
        {"h": 190.0, "scale": 28.0, "fp": 6.0},
        {"h": 300.0, "scale": 55.0, "fp": 11.5},
    ]
    ne = np.zeros_like(zz, dtype=float)
    for layer in layers:
        layer_height = float(layer["h"])
        if layer_height < 130.0:
            e_tilt = 13.5 * np.tanh((0.45 * xx - 0.75 * yy) / 30.0)
            e_tilt += 4.0 * np.sin(2.0 * np.pi * (0.65 * xx + 0.35 * yy) / 55.0)
            layer_height = layer_height + e_tilt
        ne += ne_from_fp_mhz(layer["fp"]) * np.exp(-((zz - layer_height) / layer["scale"]) ** 2)

    # Weak E-region structure: small-amplitude tilted waves plus an asymmetric
    # packet/depletion pair to move E-region returns a few km off center.
    e_env = np.exp(-((zz - 112.0) / 13.0) ** 2)
    e_phase_1 = 2.0 * np.pi * (0.70 * xx + 0.45 * yy) / 36.0 + 0.3
    e_phase_2 = 2.0 * np.pi * (-0.40 * xx + 0.90 * yy) / 28.0 - 0.5
    e_packet = np.exp(-((xx - 12.0) / 24.0) ** 2 - ((yy + 4.0) / 18.0) ** 2)
    e_depletion = np.exp(-((xx + 16.0) / 22.0) ** 2 - ((yy - 8.0) / 20.0) ** 2)
    e_gradient = np.tanh((0.55 * xx - 0.83 * yy) / 28.0)
    ne *= 1.0 + e_env * (
        0.11 * np.sin(e_phase_1)
        + 0.075 * np.sin(e_phase_2)
        + 0.075 * e_gradient
    )
    ne *= 1.0 + 0.16 * e_packet * e_env
    ne *= 1.0 - 0.13 * e_depletion * e_env

    # Weak, broad F2 structure and stronger confined F1 structure. The F1
    # packets are intentionally smaller scale to mimic eclipse-time fine
    # structure seen in the VIPIR skymaps.
    f2_env = np.exp(-((zz - 300.0) / 75.0) ** 2)
    f2_phase_1 = 2.0 * np.pi * (0.30 * xx + 0.95 * yy) / 135.0
    f2_phase_2 = 2.0 * np.pi * (-0.70 * xx + 0.42 * yy) / 105.0 + 0.8
    f2_phase_3 = 2.0 * np.pi * (0.82 * xx - 0.55 * yy) / 72.0 - 1.2
    ne *= 1.0 + f2_env * (
        0.08 * np.sin(f2_phase_1)
        + 0.06 * np.sin(f2_phase_2)
        + 0.035 * np.sin(f2_phase_3)
    )

    f1_env = np.exp(-((zz - 190.0) / 32.0) ** 2)
    f1_phase_1 = 2.0 * np.pi * (0.35 * xx + 1.05 * yy) / 58.0 + 0.4
    f1_phase_2 = 2.0 * np.pi * (-0.95 * xx + 0.35 * yy) / 42.0 - 0.7
    f1_phase_3 = 2.0 * np.pi * (0.85 * xx - 0.55 * yy) / 30.0 + 1.1
    f1_phase_4 = 2.0 * np.pi * (-0.20 * xx - 1.10 * yy) / 24.0 + 2.0
    f1_tid = (
        0.22 * np.sin(f1_phase_1)
        + 0.16 * np.sin(f1_phase_2)
        + 0.10 * np.sin(f1_phase_3)
        + 0.055 * np.sin(f1_phase_4)
    )
    ne *= 1.0 + f1_env * f1_tid

    u1 = 0.45 * (xx - 12.0) + 0.89 * (yy - 22.0)
    v1 = -0.89 * (xx - 12.0) + 0.45 * (yy - 22.0)
    u2 = -0.72 * (xx + 38.0) + 0.69 * (yy - 8.0)
    v2 = -0.69 * (xx + 38.0) - 0.72 * (yy - 8.0)
    u3 = 0.28 * (xx - 4.0) - 0.96 * (yy + 48.0)
    v3 = 0.96 * (xx - 4.0) + 0.28 * (yy + 48.0)
    packet_1 = np.exp(-(u1 / 58.0) ** 2 - (v1 / 16.0) ** 2)
    packet_2 = np.exp(-(u2 / 45.0) ** 2 - (v2 / 14.0) ** 2)
    packet_3 = np.exp(-(u3 / 36.0) ** 2 - (v3 / 12.0) ** 2)
    ne *= 1.0 + 0.26 * packet_1 * np.exp(-((zz - 180.0) / 25.0) ** 2)
    ne *= 1.0 - 0.20 * packet_2 * np.exp(-((zz - 205.0) / 30.0) ** 2)
    ne *= 1.0 + 0.15 * packet_3 * np.exp(-((zz - 165.0) / 22.0) ** 2)
    ne = np.clip(ne, 0.0, None)

    profile = RT3DProfile(
        lats=np.asarray(lats, dtype=float),
        lons=np.asarray(lons, dtype=float),
        alts_km=alts,
        time=EVENT_TIME,
        ne_m3=ne,
        source="synthetic_tid_bump",
    )
    return RT3D(profile=profile)


def azimuth_fan(centers, width_deg, step_deg):
    offsets = np.arange(-width_deg, width_deg + 0.001, step_deg)
    azimuths = []
    for center in centers:
        azimuths.extend(center + offsets)
    return np.mod(np.asarray(azimuths, dtype=float), 360.0)


def trace_region(rt, region, mode="quick", ray_paths=None, path_stride=1):
    freqs = REGIONS[region]["freqs"]
    if mode == "quick":
        representative = {"e": [2.2], "f": [4.5, 9.0, 12.5]}
        freqs = np.asarray(representative[region], dtype=float)
        azimuths = np.arange(0.0, 360.0, 30.0)
        elevations = np.asarray([78.0, 82.0, 85.0, 87.0, 89.0])
        max_step_km = 4.0
    elif mode == "dense":
        representative = {"e": [2.2], "f": [4.0, 6.0, 9.0, 12.0]}
        freqs = np.asarray(representative[region], dtype=float)
        azimuths = np.arange(0.0, 360.0, 15.0)
        elevations = np.asarray([77.0, 79.0, 81.0, 83.0, 85.0, 87.0, 88.5])
        max_step_km = 3.0
    else:
        azimuths = np.arange(0.0, 360.0, 10.0)
        elevations = np.arange(74.0, 90.0, 1.0)
        max_step_km = 4.0

    rows = []
    path_counter = 0
    for freq_mhz in freqs:
        for elev in elevations:
            for az in azimuths:
                try:
                    ray = rt.oblique_trace(
                        freq_hz=float(freq_mhz) * 1e6,
                        elevation_deg=float(elev),
                        azimuth_deg=float(az),
                        coordinate_system="cartesian",
                        solver="hamiltonian",
                        x0_km=0.0,
                        y0_km=0.0,
                        z0_km=0.0,
                        s_max_km=700.0,
                        max_step_km=max_step_km,
                        mode="O",
                    )
                except Exception as exc:
                    rows.append(
                        {
                            "region": region,
                            "freq_mhz": freq_mhz,
                            "elevation_deg": elev,
                            "azimuth_deg": az,
                            "status": f"error:{exc}",
                            "x_turn_km": np.nan,
                            "y_turn_km": np.nan,
                            "z_turn_km": np.nan,
                            "x_ground_km": np.nan,
                            "y_ground_km": np.nan,
                        }
                    )
                    continue

                z = np.asarray(ray.z_km, dtype=float)
                x = np.asarray(ray.x_km, dtype=float)
                y = np.asarray(ray.y_km, dtype=float)
                if z.size:
                    k = int(np.nanargmax(z))
                    x_turn, y_turn, z_turn = float(x[k]), float(y[k]), float(z[k])
                    x_ground, y_ground = float(x[-1]), float(y[-1])
                else:
                    x_turn = y_turn = z_turn = x_ground = y_ground = np.nan

                rows.append(
                    {
                        "region": region,
                        "freq_mhz": freq_mhz,
                        "elevation_deg": elev,
                        "azimuth_deg": az,
                        "status": str(ray.status),
                        "x_turn_km": x_turn,
                        "y_turn_km": y_turn,
                        "z_turn_km": z_turn,
                        "x_ground_km": x_ground,
                        "y_ground_km": y_ground,
                        "group_path_km": float(ray.group_path_km),
                    }
                )
                if ray_paths is not None and path_counter % int(path_stride) == 0:
                    ray_paths.append(
                        {
                            "region": region,
                            "freq_mhz": float(freq_mhz),
                            "elevation_deg": float(elev),
                            "azimuth_deg": float(az),
                            "status": str(ray.status),
                            "x_km": np.asarray(ray.x_km, dtype=float),
                            "y_km": np.asarray(ray.y_km, dtype=float),
                            "z_km": np.asarray(ray.z_km, dtype=float),
                        }
                    )
                path_counter += 1
    return pd.DataFrame(rows)


def plot_modeled_skymaps(results, out_path):
    # apply_jgr_style(fontsize=9, tick_size=8, line_width=0.8)
    fig, axes = plt.subplots(1, 2, figsize=(5.35, 2.65), sharex=True, sharey=True)

    sc = None
    for ax, region in zip(axes, ["e", "f"]):
        df = results[(results["region"] == region) & np.isfinite(results["x_turn_km"])]
        returned = df["status"].eq("ground")
        if np.any(~returned):
            ax.scatter(
                df.loc[~returned, "x_turn_km"],
                df.loc[~returned, "y_turn_km"],
                s=8,
                color="0.65",
                alpha=0.35,
                linewidths=0,
                label="not returned",
            )
        if np.any(returned):
            sc = ax.scatter(
                df.loc[returned, "x_turn_km"],
                df.loc[returned, "y_turn_km"],
                c=df.loc[returned, "freq_mhz"],
                s=12,
                cmap="viridis",
                vmin=1.5,
                vmax=14.0,
                alpha=0.8,
                linewidths=0,
                label="returned",
            )
        ax.axhline(0, color="0.35", lw=0.6, ls="--")
        ax.axvline(0, color="0.35", lw=0.6, ls="--")
        ax.set_title(REGIONS[region]["title"], fontsize=9)
        ax.set_xlim(-50, 50)
        ax.set_ylim(-50, 50)
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(direction="in", top=True, right=True, labelsize=8)
        ax.minorticks_on()
        ax.set_xlabel("XL (east, km)")
    axes[0].set_ylabel("YL (north, km)")

    if sc is not None:
        cax = fig.add_axes([0.91, 0.28, 0.018, 0.46])
        cb = fig.colorbar(sc, cax=cax)
        cb.set_label("Frequency (MHz)", fontsize=9)
        cb.ax.tick_params(labelsize=8)

    fig.text(
        0.5,
        0.98,
        "Synthetic 3-D NVIS ray turning locations: E versus F-region structure",
        ha="center",
        va="top",
        fontsize=10,
    )
    fig.subplots_adjust(left=0.10, right=0.87, bottom=0.18, top=0.84, wspace=0.16)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def ray_curve_score(path):
    """Score ray-path curvature/loopiness for side-view thinning."""
    x = np.asarray(path["x_km"], dtype=float)
    y = np.asarray(path["y_km"], dtype=float)
    z = np.asarray(path["z_km"], dtype=float)
    if x.size < 4 or y.size < 4 or z.size < 4:
        return 0.0
    ds = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2 + np.diff(z) ** 2)
    path_len = float(np.nansum(ds))
    chord = float(np.hypot(np.hypot(x[-1] - x[0], y[-1] - y[0]), z[-1] - z[0]))
    vx, vy, vz = np.diff(x), np.diff(y), np.diff(z)
    norm = np.sqrt(vx**2 + vy**2 + vz**2)
    valid = norm > 0
    if np.count_nonzero(valid) < 3:
        turn = 0.0
    else:
        dirs = np.column_stack([vx[valid] / norm[valid], vy[valid] / norm[valid], vz[valid] / norm[valid]])
        dots = np.sum(dirs[:-1] * dirs[1:], axis=1)
        turn = float(np.nansum(np.arccos(np.clip(dots, -1.0, 1.0))))
    excess = path_len / max(chord, 1.0)
    return excess + 0.35 * turn


def is_escaping_ray(path):
    """True for non-ground rays that do not turn back down toward ground."""
    if path["status"] == "ground":
        return False
    z = np.asarray(path["z_km"], dtype=float)
    x = np.asarray(path["x_km"], dtype=float)
    y = np.asarray(path["y_km"], dtype=float)
    if z.size < 4 or not np.all(np.isfinite(z)):
        return False
    apex = float(np.nanmax(z))
    apex_idx = int(np.nanargmax(z))
    z_end = float(z[-1])
    tail = np.diff(z[-min(8, z.size):])
    final_slope = float(np.nanmedian(tail)) if tail.size else 0.0
    if apex_idx < z.size - 2 and z_end < apex - 5.0:
        return False
    if final_slope < 0.0:
        return False
    if z_end > 430.0:
        return True
    if x.size and y.size and (abs(float(x[-1])) > 145.0 or abs(float(y[-1])) > 145.0):
        return z_end > 260.0
    return False


def ray_paths_to_dataframe(ray_paths):
    rows = []
    for path_id, path in enumerate(ray_paths):
        x = np.asarray(path["x_km"], dtype=float)
        y = np.asarray(path["y_km"], dtype=float)
        z = np.asarray(path["z_km"], dtype=float)
        if x.size == 0 or y.size == 0 or z.size == 0:
            continue
        near_station = path["status"] == "ground" and np.hypot(float(x[-1]), float(y[-1])) <= 10.0
        curve_score = ray_curve_score(path)
        for point_id, (xp, yp, zp) in enumerate(zip(x, y, z)):
            rows.append(
                {
                    "path_id": path_id,
                    "point_id": point_id,
                    "region": path["region"],
                    "freq_mhz": path["freq_mhz"],
                    "elevation_deg": path["elevation_deg"],
                    "azimuth_deg": path["azimuth_deg"],
                    "status": path["status"],
                    "near_station": near_station,
                    "curve_score": curve_score,
                    "x_km": xp,
                    "y_km": yp,
                    "z_km": zp,
                }
            )
    return pd.DataFrame(rows)


def dataframe_to_ray_paths(paths_df):
    ray_paths = []
    for _, group in paths_df.sort_values(["path_id", "point_id"]).groupby("path_id", sort=True):
        first = group.iloc[0]
        ray_paths.append(
            {
                "region": first["region"],
                "freq_mhz": float(first["freq_mhz"]),
                "elevation_deg": float(first["elevation_deg"]),
                "azimuth_deg": float(first["azimuth_deg"]),
                "status": str(first["status"]),
                "near_station": bool(first["near_station"]),
                "curve_score": float(first["curve_score"]),
                "x_km": group["x_km"].to_numpy(dtype=float),
                "y_km": group["y_km"].to_numpy(dtype=float),
                "z_km": group["z_km"].to_numpy(dtype=float),
            }
        )
    return ray_paths


def select_sideview_paths(ray_paths, max_ground_paths=110, max_escape_paths=60):
    f_paths = [path for path in ray_paths if path["region"] == "f"]
    green_paths = []
    ground_candidates = []
    escape_candidates = []
    for path in f_paths:
        x = np.asarray(path["x_km"], dtype=float)
        y = np.asarray(path["y_km"], dtype=float)
        near_station = bool(path.get("near_station", False))
        if x.size and y.size and path["status"] == "ground":
            near_station = near_station or np.hypot(float(x[-1]), float(y[-1])) <= 10.0
        path["curve_score"] = float(path.get("curve_score", ray_curve_score(path)))
        if path["status"] == "ground" and near_station:
            green_paths.append(path)
        elif path["status"] == "ground":
            ground_candidates.append(path)
        elif is_escaping_ray(path):
            escape_candidates.append(path)
        else:
            continue

    ground_candidates = sorted(ground_candidates, key=lambda item: item["curve_score"], reverse=True)[::2]
    escape_candidates = sorted(escape_candidates, key=lambda item: float(np.nanmax(item["z_km"])), reverse=True)[::2]
    black_paths = ground_candidates[:max_ground_paths] + escape_candidates[:max_escape_paths]
    return black_paths, green_paths


def plot_side_views(rt, ray_paths, out_path):
    apply_jgr_style(fontsize=9, tick_size=8, line_width=0.8)
    x_km, y_km = latlon_to_local_axes(rt.profile.lats, rt.profile.lons)
    alts = np.asarray(rt.profile.alts_km, dtype=float)
    fp = plasma_freq_mhz(rt.profile.ne_m3)
    ix0 = int(np.nanargmin(np.abs(x_km)))
    iy0 = int(np.nanargmin(np.abs(y_km)))
    fp_xz = fp[iy0, :, :].T
    fp_yz = fp[:, ix0, :].T

    fig, axes = plt.subplots(1, 2, figsize=(6.35, 3.05), sharey=True)
    panels = [
        (axes[0], x_km, fp_xz, "XL side view", "XL (east, km)", "x_km"),
        (axes[1], y_km, fp_yz, "YL side view", "YL (north, km)", "y_km"),
    ]
    black_paths, green_paths = select_sideview_paths(ray_paths)

    mesh = None
    for ax, axis_km, fp_slice, title, xlabel, coord_key in panels:
        mesh = ax.pcolormesh(
            axis_km,
            alts,
            fp_slice,
            cmap="magma_r",
            vmin=1.0,
            vmax=12.0,
            shading="auto",
        )
        for path in black_paths:
            coord = np.asarray(path[coord_key], dtype=float)
            z = np.asarray(path["z_km"], dtype=float)
            if coord.size < 2 or z.size < 2:
                continue
            ax.plot(coord, z, color="black", lw=0.45, alpha=0.50, ls="-")
        for path in green_paths:
            coord = np.asarray(path[coord_key], dtype=float)
            z = np.asarray(path["z_km"], dtype=float)
            if coord.size < 2 or z.size < 2:
                continue
            ax.plot(coord, z, color="#00a33a", lw=0.45, alpha=0.95, ls="-", zorder=5)
        ax.axvline(0, color="white", lw=0.7, ls="--", alpha=0.8)
        ax.set_title(title, fontsize=9)
        ax.set_xlim(-120, 120)
        ax.set_ylim(0, 500)
        ax.set_xlabel(xlabel)
        ax.tick_params(direction="in", top=True, right=True, labelsize=8)
        ax.minorticks_on()
    axes[0].set_ylabel("Height (km)")

    cax = fig.add_axes([0.90, 0.25, 0.018, 0.52])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label(r"$f_0$ (MHz)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.text(
        0.5,
        0.98,
        "F-region ray paths through altered electron density",
        ha="center",
        va="top",
        fontsize=10,
    )
    fig.subplots_adjust(left=0.10, right=0.87, bottom=0.17, top=0.84, wspace=0.16)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Run synthetic 3-D NVIS raytrace skymaps.")
    parser.add_argument(
        "--region",
        choices=["all", "e", "f"],
        default="all",
        help="Region to run. Default runs E and combined F regions.",
    )
    parser.add_argument("--quick", action="store_true", help="Use a small angular sweep.")
    parser.add_argument("--dense", action="store_true", help="Use denser plasma grid and ray sweep.")
    parser.add_argument("--sideview", action="store_true", help="Save XL/YL side-view ray projections.")
    parser.add_argument(
        "--sideview-from-csv",
        type=Path,
        default=None,
        help="Draw side-view figure from a saved ray-path CSV instead of tracing rays.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mode = "dense" if args.dense else ("quick" if args.quick else "full")
    rt = build_synthetic_profile(dense=args.dense)
    if args.sideview_from_csv is not None:
        ray_paths = dataframe_to_ray_paths(pd.read_csv(args.sideview_from_csv))
        side_path = OUTPUT_DIR / f"FigureRT3D_NVIS_TID_Bump_sideview_{mode}.png"
        plot_side_views(rt, ray_paths, side_path)
        print(f"Saved: {side_path}")
        return

    regions = ["e", "f"] if args.region == "all" else [args.region]
    ray_paths = [] if args.sideview else None
    path_stride = 1
    frames = [
        trace_region(
            rt,
            region,
            mode=mode,
            ray_paths=ray_paths,
            path_stride=path_stride,
        )
        for region in regions
    ]
    results = pd.concat(frames, ignore_index=True)

    suffix = mode
    csv_path = OUTPUT_DIR / f"FigureRT3D_NVIS_TID_Bump_{args.region}_{suffix}.csv"
    results.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")
    print(results.groupby(["region", "status"]).size().to_string())

    if args.region == "all":
        fig_path = OUTPUT_DIR / f"FigureRT3D_NVIS_TID_Bump_{suffix}.png"
        plot_modeled_skymaps(results, fig_path)
        print(f"Saved: {fig_path}")
        if args.sideview and ray_paths:
            ray_path_csv = OUTPUT_DIR / f"FigureRT3D_NVIS_TID_Bump_raypaths_{suffix}.csv"
            ray_paths_to_dataframe(ray_paths).to_csv(ray_path_csv, index=False)
            print(f"Saved: {ray_path_csv}")
            side_path = OUTPUT_DIR / f"FigureRT3D_NVIS_TID_Bump_sideview_{suffix}.png"
            saved_paths = dataframe_to_ray_paths(pd.read_csv(ray_path_csv))
            plot_side_views(rt, saved_paths, side_path)
            print(f"Saved: {side_path}")


if __name__ == "__main__":
    main()
