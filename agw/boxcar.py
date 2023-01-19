#!/usr/bin/env python

"""boxcar.py: module is dedicated to run all analysis and filtering."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import copy

import numpy as np
from fetch_fit_data import Beam, Gate, Scan
from loguru import logger
from scipy import signal
from scipy import stats as st
from scipy.stats import beta


def create_gaussian_weights(mu, sigma, _kernel=3, base_w=5):
    """
    Method used to create gaussian weights
    mu: n 1D list of all mean values
    sigma: n 1D list of sigma matrix
    """
    _k = (_kernel - 1) / 2
    _kNd = np.zeros((3, 3, 3))
    for i in range(_kernel):
        for j in range(_kernel):
            for k in range(_kernel):
                _kNd[i, j, k] = np.exp(
                    -(
                        (float(i - _k) ** 2 / (2 * sigma[0] ** 2))
                        + (float(j - _k) ** 2 / (2 * sigma[1] ** 2))
                        + (float(k - _k) ** 2 / (2 * sigma[2] ** 2))
                    )
                )
    _kNd = np.floor(_kNd * base_w).astype(int)
    return _kNd


class BxFilter(object):
    """Class to filter data - Boxcar median filter."""

    def __init__(self, thresh=0.7, w=None, verbose=True):
        """
        initialize variables

        thresh: Threshold of the weight matrix
        w: Weight matrix
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
        self.verbose = verbose
        if self.verbose:
            logger.info(f"Initialize median filter with default weigth matrix")
        return

    def _discard_repeting_beams(self, scan, ch=True):
        """
        Discard all more than one repeting beams
        scan: SuperDARN scan
        """
        oscan = Scan(scan.stime, scan.etime, scan.s_mode)
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

    def doFilter(self, i_scans, comb=False, gflg_type=-1):
        """
        Median filter based on the weight given by matrix (3X3X3) w, and threshold based on thresh

        i_scans: 3 consecutive radar scans
        comb: combine beams
        """
        if comb:
            scans = []
            for s in i_scans:
                scans.append(self._discard_repeting_beams(s))
        else:
            scans = i_scans
        logger.info(f"Sample time {scans[1].stime}")
        self.scans = scans
        if len(scans) == 3:
            w = self.w
            oscan = Scan(scans[1].stime, scans[1].etime, scans[1].s_mode)
            if w is None:
                w = np.array(
                    [
                        [[1, 2, 1], [2, 3, 2], [1, 2, 1]],
                        [[2, 3, 2], [3, 5, 3], [2, 3, 2]],
                        [[1, 2, 1], [2, 3, 2], [1, 2, 1]],
                    ]
                )
            l_bmnum, r_bmnum = scans[1].beams[0].bmnum, scans[1].beams[-1].bmnum

            for b in scans[1].beams:
                bmnum = b.bmnum
                beam = Beam()
                beam.copy(b)

                for key in beam.__dict__.keys():
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
                                        tbm, ind, gflg_type=gflg_type
                                    )
                                else:
                                    box[j][k + 1][n + 1] = 0
                    pts = 0.0
                    tot = 0.0
                    v, w_l, p_l, gfx = list(), list(), list(), list()

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
                                        v.append(bx.v)
                                        w_l.append(bx.w_l)
                                        p_l.append(bx.p_l)
                                        gfx.append(bx.gflg)
                    if pts / tot >= self.thresh:  # check if we meet the threshold
                        beam.slist.append(r)
                        beam.v.append(np.median(v))
                        beam.w_l.append(np.median(w_l))
                        beam.p_l.append(np.median(p_l))
                oscan.beams.append(beam)
            if len(oscan.beams) > 0:
                oscan.update_time()
                sorted(oscan.beams, key=lambda bm: bm.bmnum)
        else:
            oscan = Scan(None, None, "")
        return oscan


if __name__ == "__main__":
    create_gaussian_weights(mu=[0, 0, 0], sigma=[3, 3, 3], base_w=7)
    filter = Filter()
