import datetime as dt
from pathlib import Path

import numpy as np

def _load_mat(mat_path: Path| str = "data/dvTEC/vTECdata_1800.mat") -> dict:
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
        data = {"tec": np.array(h5f[k]) for k in h5f.keys()}
    data["lat"],data["lon"] = np.arange(27,48,0.125), np.arange(-103,-75,0.125)
    return data


if __name__ == "__main__":
    mat_file = Path("data/dvTEC/vTECdata_1800.mat")
    snap = _load_mat(mat_file)
    print(f"Loaded data keys: {list(snap.keys())}")
    print(
        "Shape of precip_at_time:", snap["tec"].shape, snap["lat"].shape, snap["lon"].shape)