#!/usr/bin/env python

"""to_netCDF.py: utility module to stored as a netCDF file."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import sys

from loguru import logger

sys.path.extend(["py/txUtils/", "py/tid/", "py/davitPy/"])
import tidUtils
from fetchUtils import FetchData
from rtiUtils import RTI, plot_FFT
from filterUtils import fft_analysis
from tidFilters import filterBasedOnLimits, fftBasedOnLimits

# Initializations
CASES = [
    "filter",
]
case = "filter"
generateRTI = False
rads = ["fhw"]
dates = [
    dt.datetime(2022, 12, 20),
]
fdMap = {}
bm_gate_cells = [(7, 1170), (7, 1215), (7, 1260), (7, 1305)]


def plot_RTI(
    rad,
    dat,
    date_range,
    beams=[],
    nGates=100,
    angle_th=100.0,
    vhm=None,
    tec_mat_file=None,
    tec_param="cdvTECgrid2",
    power_vlim=[-5, 5],
    tec_vlim=[-0.3, 0.3],
):
    """
    Plot RTI plots by beams
    """
    logger.info("Inside RTI plots")
    beams = beams if beams and len(beams) > 0 else dat.bmnum.unique()

    if tec_mat_file:
        tec, tec_times = tidUtils.read_tec_mat_files(tec_mat_file)

    for b in beams:
        logger.info(f" Beam:{b}")
        file = tidUtils.get_folder(date_range[0]) + f"/{rad}-filter-{'%02d'%b}.png"
        rt = RTI(
            100,
            date_range,
            None,
            date_range,
            f"{date_range[0].strftime('%Y-%m-%d')}/{rad}/{b}",
            num_subplots=3,
            angle_th=angle_th,
            vhm=None,
            ylim=[180, 2500],
        )
        rt.addParamPlot(
            dat,
            b,
            "",
            xlabel="",
            cbar=True,
            plot_fov=False,
            vlim=power_vlim,
        )
        if tec_mat_file:
            rt.ovearlay_TEC(
                tec,
                tec_times,
                b,
                plot_fov=False,
                xlabel="",
                tec_param=tec_param,
                vlim=tec_vlim,
                cmap="magma",
            )
            axis = rt.ovearlay_TEC(
                tec,
                tec_times,
                b,
                cbar_xOff=0.15,
                tec_param=tec_param,
                plot_fov=False,
                vlim=tec_vlim,
                cmap="magma",
            )
            rt.addParamPlot(
                dat,
                b,
                "",
                xlabel="",
                plot_fov=False,
                ax=axis,
                vlim=power_vlim,
            )
        rt.save(file)
        rt.close()
    return


# Conduct RTI analysis with filter
if case == "filter":
    for d in dates:
        for rad in rads:
            fd = FetchData.fetch(
                rad,
                [d, d + dt.timedelta(1)],
            )
            rsep, frang = (fd.frame.rsep.tolist()[0], fd.frame.frang.tolist()[0])
            if generateRTI:
                fd.plot_RTI(
                    tec_mat_file=f"data/{d.strftime('%Y-%m-%d')}/{rad}geom.mat",
                    tec_param="cdvTECgrid2",
                    power_vlim=[20, 50],
                    tec_vlim=[-0.3, 0.3],
                )
            print(sorted(fd.frame.srange.unique()))
            for bm_gate_cell in bm_gate_cells:
                o = fftBasedOnLimits(
                    fd,
                    bm_gate_cell[0],
                    bm_gate_cell[1],
                    [d + dt.timedelta(hours=17), d + dt.timedelta(hours=21)],
                    timeRes=60,
                )
    #             dat = filterBasedOnLimits(
    #                 fd,
    #                 [d, d + dt.timedelta(1)],
    #                 [20, 60],
    #                 [d + dt.timedelta(hours=17), d + dt.timedelta(hours=21)],
    #                 timeRes=60,
    #                 numtaps=101,
    #                 cutoffs_min=[15, 40],
    #                 pfilt_thresh=0.3,
    #             )
    #            dat["rsep"], dat["frang"], dat["srange"] = rsep, frang, rsep*dat.slist+frang
                #dat = dat.dropna()
    #             plot_RTI(
    #                 rad,
    #                 dat,
    #                 [d + dt.timedelta(hours=17), d + dt.timedelta(hours=21)],
    #                 tec_mat_file=f"data/{d.strftime('%Y-%m-%d')}/{rad}geom.mat",
    #             )
    #             o = fft_analysis(
    #                 dat,
    #                 bm_gate_cell[1],
    #                 bm_gate_cell[0],
    #                 T=60,
    #             )
    #             print(o)
                fname = f"data/{d.strftime('%Y-%m-%d')}/{rad}-{bm_gate_cell[0]}-{bm_gate_cell[1]}-fft.png"
                plot_FFT(o, bm_gate_cell[0], bm_gate_cell[1], fname)