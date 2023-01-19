#!/usr/bin/env python

"""
    music.py: music processing module
    A module for running the MUltiple SIgnal Classification (MUSIC) algorithm for the detection of
    MSTIDs and wave-like structures in SuperDARN data.
    For usage examples, please see the iPython notebooks included in the docs folder of the DaViTPy distribution.
    References
    ----------
    See Samson et al. [1990] and Bristow et al. [1994] for details regarding the MUSIC algorithm and SuperDARN-observed MSTIDs.
    Bristow, W. A., R. A. Greenwald, and J. C. Samson (1994), Identification of high-latitude acoustic gravity wave sources
        using the Goose Bay HF Radar, J. Geophys. Res., 99(A1), 319-331, doi:10.1029/93JA01470.
    Samson, J. C., R. A. Greenwald, J. M. Ruohoniemi, A. Frey, and K. B. Baker (1990), Goose Bay radar observations of Earth-reflected,
        atmospheric gravity waves in the high-latitude ionosphere, J. Geophys. Res., 95(A6), 7693-7709, doi:10.1029/JA095iA06p07693.
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import datetime as dt
import pandas as pd
from scipy.interpolate import interp1d
from scipy import signal, fftpack
from types import SimpleNamespace
import json
from tqdm import tqdm
import numpy as np
from loguru import logger

import pydarn

import sys
sys.path.extend(["agw/", "agw/tid/"])
from plots import RTI, PlotKarr, plot_Dlm
from fetch_fit_data import FetchData


def greatCircleAzm(lat1, lon1, lat2, lon2):
    """
    Calculates the azimuth from the coordinates of a start point to and end
    point along a great circle path.
    """
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    dlon = lon2 - lon1
    y = np.sin(dlon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    azm = np.degrees(np.arctan2(y,x))
    return azm

def greatCircleDist(lat1, lon1, lat2, lon2):
    """
    Calculates the distance in radians along a great circle path between two
    points.
    """
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlat = (lat2 - lat1) / 2.0
    dlon = (lon2 - lon1) / 2.0
    a = np.sin(dlat)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon)**2
    radDist = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    return radDist

class Dataset(object):
    """
    This class is the basic container for holding MUSIC data sets.
    """

    def __init__(
        self,
        rad,
        dates,
        frame,
        param="p_l"
    ):
        self.rad = rad
        self.dates = dates
        self.param = param
        self.load_params()
        
        self.frames = {"DS.raw": frame if frame else self._fetch()}
        return
    
    def load_params(self):
        """
        Load parameters from params.json
        """
        hdw_data = pydarn.read_hdw_file(self.rad)
        self.lats, self.lons = pydarn.Coords.GEOGRAPHIC(hdw_data.stid)
        self.beams, self.gates = np.arange(1, self.lats.shape[1]+1), np.arange(1, self.lats.shape[0]+1)
        with open("config/params.json") as f:
            self.conf = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
        self.base = self.conf.files.base.format(
            date=self.dates[0].strftime("%Y-%m-%d"), run_id=self.conf.run_id
        )
        if not os.path.exists(self.base):
            os.makedirs(self.base, exist_ok=True)
        return
    
    def _fetch(self):
        """
        Fetch radar data from repository
        """
        logger.info(" Fetching data...")
        if not os.path.exists(self.fetch_file("DS.raw")):
            logger.info(
                f" Read radar file {self.rad} for {[d.strftime('%Y.%m.%dT%H.%M') for d in self.dates]}"
            )
            self.fd = FetchData(self.rad, self.dates, ftype=self.conf.ftype)
            _, scans, self.data_exists = self.fd.fetch_data(by="scan")
            if self.data_exists:
                scan_time = (scans[0].etime-scans[0].stime).total_seconds()
                frame = self.fd.scans_to_pandas(scans)
                if len(frame) > 0:
                    frame["srange"] = frame.frang + (
                        frame.slist * frame.rsep
                    )
                    frame["intt"] = (
                        frame["intt.sc"] + 1.0e-6 * frame["intt.us"]
                    )
                    frame = frame.progress_apply(
                        self.__get_location__, axis=1
                    )
                    frame["scan_time"] = scan_time
                    self._save(frame, "raw")
                else:
                    logger.info(f" Radar file does not exist!")
        else:
            frame = pd.read_csv(self.fetch_file("DS.raw"), parse_dates=["time"])
        return frame
    
    def _load(self, p):
        """
        Load all other files
        """
        if os.path.exists(self.fetch_file(p)):
            self.frames[p] = pd.read_csv(self.fetch_file(p), parse_dates=["time"])
        return
    
    def _save(self, df, p):
        """
        Print data into parquet/gzip file for later processing.
        """
        df.to_csv(self.fetch_file(p), index=False, header=True)
        return
    
    def fetch_file(self, p):
        """
        Return filename
        """
        stime, etime = self.dates[0].strftime("%H%M"), self.dates[-1].strftime("%H%M")
        file = self.base + self.conf.files.csv.format(
            rad=self.rad, stime=stime, etime=etime, kind=p
        )
        return file

    def __get_location__(self, row):
        """
        Get AACGMV2 magnetic location
        """
        lat, lon, dtime = (
            self.lats[row["slist"], row["bmnum"]],
            self.lons[row["slist"], row["bmnum"]],
            row["time"],
        )
        row["glat"], row["glon"] = lat, lon
        return row
    
    
class MUSIC(object):
    """
    This class is the basic container for holding MUSIC analysis.

    References
    ----------
    Bristow, W. A., R. A. Greenwald, and J. C. Samson (1994), Identification of high-latitude acoustic gravity wave sources
        using the Goose Bay HF Radar, J. Geophys. Res., 99(A1), 319-331, doi:10.1029/93JA01470.
    """
    
    def __init__(
        self, 
        ds, 
        timeres=120, 
        numtaps=101,
        beams=[]
    ):
        self.dataset = ds
        self.timeres = timeres
        self.numtaps = numtaps
        self.beams = self.dataset.frames["DS.raw"].bmnum.unique()
        return
    
    def applyLimits(
        self, 
        gateLimits, 
        timeLimits
    ):
        """
        Filter by gate and time limits
        """
        td = dt.timedelta(seconds=(int(self.numtaps*self.timeres/2.)))
        timeLimits = [
            timeLimits[0] - td,
            timeLimits[1] + td
        ]
        raw = self.dataset.frames["DS.raw"].copy()
        self.dataset.frames["DS.sel"] = raw[
            (raw.slist>=gateLimits[0])
            & (raw.slist<=gateLimits[1])
            & (raw.time>=timeLimits[0])
            & (raw.time<=timeLimits[1])
        ]
        self.gateLimits = gateLimits
        self.timeLimits = timeLimits
        self.dtime = [ self.timeLimits[0] + dt.timedelta(seconds=self.timeres*i)
            for i in range(int((self.timeLimits[1]-self.timeLimits[0]).total_seconds()/self.timeres))
        ]
        return
    
    def beamInterpolation(self):
        """
        Interolate by each beam along gates
        """
        logger.info("Running beam interpolation..")
        iGates = np.arange(self.gateLimits[0], self.gateLimits[1])
        sel = self.dataset.frames["DS.sel"].copy()
        o = pd.DataFrame()
        for t in sel.time.unique():
            for b in sel.bmnum.unique():
                x = sel[
                    (sel.bmnum == b)
                    & (sel.time == t)
                ]
                if len(x) > 2:
                    gates = np.array(x.slist)
                    val = np.array(x[self.dataset.param])
                    intFn = interp1d(
                        gates, val, 
                        bounds_error=False, 
                        fill_value=0
                    )
                    x = pd.DataFrame()
                    x["slist"], x[self.dataset.param], x["time"], x["bmnum"] = (
                        iGates,
                        intFn(iGates),
                        t,
                        b
                    )
                    o = pd.concat([x, o])
        self.dataset.frames["DS.intp_by_beam"] = o
        return
    
    def timeInterpolation(self):
        """
        Interolate along each beam by time
        """
        logger.info("Running time interpolation..")
        dtime = np.array([
            t.hour * 3600 + t.minute * 60 + t.second
            for t in self.dtime
        ])
        intp = self.dataset.frames["DS.intp_by_beam"].copy()
        o = pd.DataFrame()
        for g in intp.slist.unique():
            for b in intp.bmnum.unique():
                xp = intp[
                    (intp.bmnum == b)
                    & (intp.slist == g)
                ]
                if len(xp) > 0:
                    time = np.array(
                        xp.time.apply(
                            lambda t: t.hour * 3600 + t.minute * 60 + t.second
                        )
                    )
                    val = np.array(xp[self.dataset.param])
                    intFn = interp1d(
                        time, val, 
                        bounds_error=False,
                        fill_value=0
                    )
                    x = pd.DataFrame()
                    x[self.dataset.param], x["time"] = (
                        intFn(dtime),
                        self.dtime
                    )
                    x["slist"], x["bmnum"] = (g,b)
                    o = pd.concat([x, o])
        self.dataset.frames["DS.intp_by_time"] = o
        return
    
    def addfilters(
        self,  
        cutoff_low, 
        cutoff_high, 
        width=None, 
        window="blackman",
        pass_zero=True, 
        scale=True
    ):
        """
        Create filters for band-pass filtering
        """
        logger.info("Filtering dataset..")
        intp = self.dataset.frames["DS.intp_by_time"].copy()
        nyq = self.nyquistFrequency()
        lp = signal.firwin(
            numtaps=self.numtaps, 
            cutoff=cutoff_high, 
            width=width, 
            window=window, 
            pass_zero=pass_zero, 
            scale=scale, 
            nyq=nyq
        )
        hp = -signal.firwin(
            numtaps=self.numtaps, 
            cutoff=cutoff_low, 
            width=width, 
            window=window, 
            pass_zero=pass_zero, 
            scale=scale, 
            nyq=nyq
        )
        d = -(lp+hp)
        d[int(self.numtaps/2)] = d[int(self.numtaps/2)] + 1
        d = -1.*d
        self.cutoff_low = cutoff_low
        self.cutoff_high = cutoff_high
        self.nyq = nyq
        self.ir = d
        
        o = pd.DataFrame()
        shift = np.int32(-np.floor(len(self.ir)/2.))
        for b in intp.bmnum.unique():
            for g in intp.slist.unique():
                xp = intp[
                    (intp.bmnum == b)
                    & (intp.slist == g)
                ]
                if len(xp) > 0:
                    xp[self.dataset.param] = xp[self.dataset.param].interpolate(method="nearest")
                    tmp = signal.lfilter(self.ir,[1.0],xp[self.dataset.param])
                    tmp = np.roll(tmp, shift)
                    x = pd.DataFrame()
                    x[self.dataset.param], x["time"] = (
                        tmp[:],
                        xp.time,
                    )
                    x["slist"], x["bmnum"] = (g,b)
                    o = pd.concat([x, o])
        self.dataset.frames["DS.fil"] = o
        return
    
    def nyquistFrequency(self):
        """
        Calculate nyquist frequency
        """
        nyq = float(1. / (2*self.timeres))
        return nyq
    
    def calculateFFT(self):
        """
        Calculate FFT
        """
        logger.info("Running FFT..")
        nyq = self.nyquistFrequency()
        freq_ax = np.arange(len(self.dtime), dtype="f8")
        freq_ax = (freq_ax / max(freq_ax)) - 0.5
        freq_ax = freq_ax * 2. * nyq
        
        fil = self.dataset.frames["DS.fil"].copy()
        o = pd.DataFrame()
        for b in fil.bmnum.unique():
            for g in fil.slist.unique():
                xp = fil[
                    (fil.bmnum == b)
                    & (fil.slist == g)
                ]
                if len(xp) > 0:
                    val = fftpack.fftshift(
                        fftpack.fft(np.array(xp[self.dataset.param]))
                    ) / len(xp)
                    x = pd.DataFrame()
                    x[self.dataset.param], x["time"], x["fft"], x["frq"] = (
                        xp[self.dataset.param],
                        xp.time,
                        val,
                        freq_ax
                    )
                    x["slist"], x["bmnum"] = (g, b)
                    o = pd.concat([x, o])
        self.dataset.frames["DS.fft"] = o
        self.freqVec = freq_ax
        return
    
    def determineRelativePosition(self, altitude=250.):
        """
        Finds the center cell of the field-of-view of a musicArray data object.
        The range, azimuth, x-range, and y-range from the center to each cell in the FOV
        is calculated and saved to the FOV object.
        """
        logger.info("Compiling relative positions..")
        Re = 6378.1
        ctrBeamInx = int(self.dataset.lats.shape[1]/2)
        ctrGateInx = int(self.dataset.lats.shape[0]/2)
        lat1 = np.zeros_like(self.dataset.lats) + self.dataset.lats[ctrGateInx,ctrBeamInx]
        lon1 = np.zeros_like(self.dataset.lons) + self.dataset.lons[ctrGateInx,ctrBeamInx]
        lat2, lon2 = self.dataset.lats, self.dataset.lons
        self.azm = greatCircleAzm(lat1,lon1,lat2,lon2)
        self.dist = (Re + altitude) * greatCircleDist(lat1,lon1,lat2,lon2)
        self.relative_x = self.dist * np.sin(np.radians(self.azm))
        self.relative_y = self.dist * np.cos(np.radians(self.azm))
        return
    
    def calculateDlm(self):
        """
        Calculate the cross-spectral matrix of a musicaArray object. FFT must already have been calculated.
        """
        logger.info("Compiling Dlm matrix..")
        fil = self.dataset.frames["DS.fil"].copy()
        nCells = len(fil.bmnum.unique()) * len(fil.slist.unique())
        posInx = np.where(self.freqVec > 0)[0]
        
        llList = []
        for b in fil.bmnum.unique():
            for g in fil.slist.unique():
                llList.append((g,b))
        
        self.llLookupTable = np.zeros([5,nCells])
        self.Dlm = np.zeros([nCells,nCells], dtype=np.complex128)
        for ll in range(nCells):
            llAI  = llList[ll]
            ew_dist = self.relative_x[llAI]
            ns_dist = self.relative_y[llAI]
            self.llLookupTable[:,ll] = [
                ll, 
                self.dataset.beams[llAI[1]], 
                self.dataset.gates[llAI[0]],
                ns_dist,
                ew_dist
            ]            
            spectL = self.get_fft_spectrum(posInx, llAI)
            for mm in range(nCells):
                mmAI = llList[mm]
                spectM = self.get_fft_spectrum(posInx,mmAI)
                self.Dlm[ll,mm] = np.sum(spectL * np.conj(spectM))
        return
    
    def get_fft_spectrum(self, posInx, llAI):
        """
        Fetch spectrum
        """
        fft = self.dataset.frames["DS.fft"]
        f = fft[
            (fft.bmnum==llAI[1])
            & (fft.slist==llAI[0])
        ]
        if len(f) > 0:
            spect = np.array(f["fft"])[posInx]
        else:
            spect = np.zeros_like(posInx)
        return spect
    
    def calculateKarr(
        self,
        kxMax=0.05,
        kyMax=0.05,
        dkx=0.001,
        dky=0.001,
        threshold=0.15
    ):
        logger.info("Compiling K Martrix..")
        # Calculate eigenvalues, eigenvectors
        eVals, eVecs = np.linalg.eig(np.transpose(self.Dlm))
        nkx = np.ceil(2*kxMax/dkx)
        if (nkx % 2) == 0: 
            nkx = nkx+1
        kxVec = kxMax * (2*np.arange(nkx)/(nkx-1) - 1)

        nky = np.ceil(2*kyMax/dky)
        if (nky % 2) == 0: nky = nky+1
        kyVec = kyMax * (2*np.arange(nky)/(nky-1) - 1)

        nkx = int(nkx)
        nky = int(nky)
        
        xm = self.llLookupTable[4,:] #x is in the E-W direction.
        ym = self.llLookupTable[3,:] #y is in the N-S direction.
        
        threshold = 0.15
        maxEval = np.max(np.abs(eVals))
        
        minEvalsInx = np.where(eVals <= threshold*maxEval)[0]
        cnt = np.size(minEvalsInx)
        maxEvalsInx = np.where(eVals > threshold*maxEval)[0]
        nSigs = np.size(maxEvalsInx)
        
        def vCalc(um,v):
            return np.dot( np.conj(um), v) * np.dot( np.conj(v), um)
        
        vList = [eVecs[:,minEvalsInx[ee]] for ee in range(cnt)]
        kArr  = np.zeros((nkx,nky),dtype=np.complex64)
        for kk_kx in range(nkx):
            kx  = kxVec[kk_kx]
            for kk_ky in range(nky):
                ky  = kyVec[kk_ky]
                um  = np.exp(1j*(kx*xm + ky*ym))
                kArr[kk_kx,kk_ky]= 1. / np.sum(list(map(lambda v: vCalc(um,v), vList)))
                
        self.karr = kArr
        self.kxVec = kxVec
        self.kyVec = kyVec
        return
    
    def close(self):
        """
        Closing the analysis by saving datasets , plots
        """
        rti = RTI(75, self.dataset.dates, fig_title="", num_subplots=3)
        rti.addParamPlot(
            self.dataset.frames["DS.raw"], 7,
            title=f"{self.dataset.rad.upper()} / Beam:{7} / Raw",
            p_step=6,
            vlines=self.timeLimits,
            hlines=self.gateLimits 
        )
        rti.save(self.dataset.base + "/RTI.png")
        rti.close()
        
        for key in self.dataset.frames.keys():
            o = self.dataset.frames[key]
            self.dataset._save(o, key)
        return
        
        
def music_analysis(
    rad, dates, frame,
    analysis_prop,
    param="p_l", timeres=120,
    numtaps=101, gateLimits=[30,45],
    
):
    tqdm.pandas()
    ds = Dataset(rad, dates, frame, param)
    music = MUSIC(ds, timeres, numtaps)
    # Create music analysis
    music.determineRelativePosition()
    music.applyLimits(
        gateLimits=analysis_prop["gateLimits"],
        timeLimits=analysis_prop["timeLimits"]
    )
    music.beamInterpolation()
    music.timeInterpolation()
    music.addfilters(
        cutoff_low=analysis_prop["cutoff_low"], 
        cutoff_high=analysis_prop["cutoff_high"]
    )
    music.calculateFFT()
    music.calculateDlm()
#     music.calculateKarr()
    # Closing analysis
    music.close()
    plot_Dlm(music.Dlm)
    #plot = PlotKarr(music.karr, music.kxVec, music.kyVec)
    #plot.save("out.png")
    return
    
    
if __name__ == "__main__":
    rad = "wal"
    dates = [
        dt.datetime(2011,5,9,8),
        dt.datetime(2011,5,9,19)
    ]
    analysis_prop = {
        "gateLimits": [30,45],
        "timeLimits": [
            dt.datetime(2011,5,9,14),
            dt.datetime(2011,5,9,16),
        ],
        "cutoff_low": 0.0003, 
        "cutoff_high": 0.0012,
    }
    music_analysis(
        rad, dates, None, analysis_prop
    )