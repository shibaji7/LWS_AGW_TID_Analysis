#!/usr/bin/env python

"""filterUtils.py: module is dedicated to run 3X3X3 boxcar filtering."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import numpy as np
from dataUtils import Beam, Gate, Scan
from loguru import logger


class Boxcar(object):
    """Class to filter data - Boxcar median filter."""

    def __init__(
        self,
        thresh=0.7,
        w=None,
        nrange_scatter=None,
    ):
        """
        initialize variables

        thresh: Threshold of the weight matrix
        w: Weight matrix
        pbnd: Lower and upper bounds of IS / GS probability
        pth: Probability of the threshold
        """
        self.thresh = thresh
        if w is None:
            w = np.array(
                [
                    [[1, 2, 1], [2, 3, 2], [1, 2, 1]],
                    [[2, 3, 2], [3, 5, 3], [2, 3, 2]],
                    [[1, 2, 1], [2, 3, 2], [1, 2, 1]],
                ]
            )
        self.w = w
        self.nrange_scatter = nrange_scatter
        return

    def __discard_repeting_beams__(self, scan, ch=True):
        """
        Discard all more than one repeting beams
        scan: SuperDARN scan
        """
        oscan = Scan()
        if ch:
            scan.beams = sorted(scan.beams, key=lambda bm: (bm.bmnum))
        else:
            scan.beams = sorted(scan.beams, key=lambda bm: (bm.bmnum, bm.time))
        bmnums = []
        for bm in scan.beams:
            if bm.bmnum not in bmnums:
                if hasattr(bm, "slist") and len(getattr(bm, "slist")) > 0:
                    oscan.beams.append(bm)
                    bmnums.append(bm.bmnum)
        if len(oscan.beams) > 0:
            oscan.update_time()
        oscan.beams = sorted(oscan.beams, key=lambda bm: bm.bmnum)
        return oscan

    def do_filter(self, scans, params_to_run_filter=["v", "w_l", "p_l", "gflg"]):
        """
        Median filter based on the weight given by matrix (3X3X3) weight,
        and threshold based on thresh

        params_to_run_filter: List of parameters to run boxcar filter
        """
        scans = [self.__discard_repeting_beams__(s) for s in scans]
        if len(scans) == 3:
            logger.info(f"Running filter for scan: {scans[1].stime}")
            w = self.w
            oscan = Scan()
            l_bmnum, r_bmnum = scans[1].beams[0].bmnum, scans[1].beams[-1].bmnum

            for b in scans[1].beams:
                bmnum = b.bmnum
                beam = Beam(self.nrange_scatter)
                beam.copy(b)

                for key in params_to_run_filter + ["slist", "srange"]:
                    if type(getattr(beam, key)) == np.ndarray:
                        setattr(beam, key, [])

                for r in range(0, b.nrang):
                    box = [
                        [[None for j in range(3)] for k in range(3)] for n in range(3)
                    ]

                    for j in range(0, 3):  # iterate through time
                        for k in range(-1, 2):  # iterate through beam
                            for n in range(-1, 2):  # iterate through gate
                                # get the scan we are working on
                                s = scans[j]
                                if s == None:
                                    continue
                                # get the beam we are working on
                                if s == None:
                                    continue
                                # get the beam we are working on
                                tbm = None
                                for bm in s.beams:
                                    if bm.bmnum == bmnum + k:
                                        tbm = bm
                                if tbm == None:
                                    continue
                                # check if target gate number is in the beam
                                if r + n in tbm.slist:
                                    ind = np.array(tbm.slist).tolist().index(r + n)
                                    box[j][k + 1][n + 1] = Gate(
                                        tbm, ind, params=params_to_run_filter
                                    )
                                else:
                                    box[j][k + 1][n + 1] = 0
                    pts = 0.0
                    tot = 0.0
                    params = {}
                    for pm in params_to_run_filter:
                        params[pm] = list()

                    for j in range(0, 3):  # iterate through time
                        for k in range(0, 3):  # iterate through beam
                            for n in range(0, 3):  # iterate through gate
                                bx = box[j][k][n]
                                if bx == None:
                                    continue
                                wt = w[j][k][n]
                                tot += wt
                                if bx != 0:
                                    pts += wt
                                    for m in range(0, wt):
                                        for pm in params_to_run_filter:
                                            params[pm].append(getattr(bx, pm))
                    if pts / tot >= self.thresh:  # check if we meet the threshold
                        for pm in params_to_run_filter:
                            getattr(beam, pm).append(np.median(params[pm]))
                        beam.slist.append(r)
                        beam.srange.append(r * beam.rsep + beam.frang)
                for k in beam.__dict__.keys():
                    if type(getattr(beam, key)) == list:
                        setattr(beam, key, np.asarray(getattr(beam, key)))
                oscan.beams.append(beam)

            if len(oscan.beams) > 0:
                oscan.update_time()
                sorted(oscan.beams, key=lambda bm: bm.bmnum)
        else:
            oscan = Scan()
        return oscan

    
def fft_analysis(dat, srange, beam, T=60):
    """
    FFT analysis for each range cell
    """
    dat = dat[
        (dat.bmnum==beam) &
        (dat.srange==srange)
    ]
    import numpy as np
    import pandas as pd
    from numpy.fft import fft, ifft
    o = pd.DataFrame()
    o["fft"], n = np.array(fft(dat.p_l)), np.arange(len(dat))
    o["freq"] = n/T
    return o