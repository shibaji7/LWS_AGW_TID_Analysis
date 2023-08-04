#!/usr/bin/env python3

"""utils.py: simulate python program utility"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import numpy as np
from scipy import array
from scipy.interpolate import interp1d


def extrap1d(x, y, kind="linear"):
    """This method is used to extrapolate 1D paramteres"""
    interpolator = interp1d(x, y, kind=kind)
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        elif x > xs[-1]:
            return ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(list(map(pointwise, array(xs))))

    return ufunclike


def interpolate_by_altitude(h, hx, param, scale="log", kind="cubic", method="intp"):
    if scale == "linear":
        pnew = (
            interp1d(h, param, kind=kind)(hx)
            if method == "intp"
            else extrap1d(h, param, kind=kind)(hx)
        )
    if scale == "log":
        pnew = (
            10 ** interp1d(h, np.log10(param), kind=kind)(hx)
            if method == "intp"
            else 10 ** extrap1d(h, np.log10(param), kind=kind)(hx)
        )
    return pnew
