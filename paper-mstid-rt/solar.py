
import os
import pandas as pd
import numpy as np
from loguru import logger
import datetime as dt



from sunpy import timeseries as ts # type: ignore
from sunpy.net import Fido # type: ignore
from sunpy.net import attrs as a # type: ignore

os.environ["OMNIDATA_PATH"] = "/home/chakras4/OMNI/"

from rays import StackPlots

class SolarDataset(object):
    """
    This class is dedicated to plot GOES, FISM, and OMNI data
    from the repo using SunPy
    """

    def __init__(self, dates):
        """
        Parameters
        ----------
        dates: list of datetime object for start and end of TS
        """
        self.dates = dates
        self.dfs = {}
        self.__load_FISM__()
        self._load_omni_()
        self.__loadGOES__()
        return
    
    def _load_omni_(self, res=1):
        import pyomnidata
        logger.info(f"OMNIDATA_PATH: {os.environ['OMNIDATA_PATH']}")
        pyomnidata.UpdateLocalData()
        self.omni = pd.DataFrame(
            pyomnidata.GetOMNI(self.dates[0].year,Res=res)
        )
        self.omni["time"] = self.omni.apply(
            lambda r: (
                dt.datetime(
                    int(str(r.Date)[:4]), 
                    int(str(r.Date)[4:6]),
                    int(str(r.Date)[6:].replace(".0","")) 
                ) 
                + dt.timedelta(hours=r.ut)
            ), 
            axis=1
        )
        self.omni = self.omni[
            (self.omni.time>=self.dates[0])
            & (self.omni.time<=self.dates[1])
        ]
        return
    
    def __load_FISM__(self):
        year, doy = (
            self.dates[0].year,
            self.dates[0].timetuple().tm_yday
        )
        url = f"https://lasp.colorado.edu/eve/data_access/eve_data/fism/flare_hr_data/{year}/"
        fname = f"FISM_60sec_{year}{doy}_v02_01.sav"
        link = url + fname
        filepath = "paper-mstid-rt/figures/"+fname
        if not os.path.exists(filepath):
            os.system(f"wget -O {filepath} {link}")
        from scipy.io import readsav
        self.fism = readsav(filepath)
        return
    
    def get_fism_spectrum_by_time(self, date):
        o = pd.DataFrame()
        i = int((date-self.dates[0]).total_seconds()/60)
        o["wv"], o["ir"] = (
            self.fism["wavelength"],
            self.fism["irradiance"][i,:]
        )
        return o

    def __loadGOES__(self):
        """
        Load GOES data from remote/local repository
        """
        self.flare = {}
        self.dfs["goes"], self.goes, self.flareHEK = pd.DataFrame(), [], None
        result = Fido.search(
            a.Time(
                self.dates[0].strftime("%Y-%m-%d %H:%M"),
                self.dates[1].strftime("%Y-%m-%d %H:%M"),
            ),
            a.Instrument("XRS") | a.hek.FL & (a.hek.FRM.Name == "SWPC"),
        )
        if len(result) > 0:
            logger.info(f"Fetching GOES ...")
            tmpfiles = Fido.fetch(result, progress=False)
            for tf in tmpfiles:
                if "avg1m" in tf:
                    self.goes.append(ts.TimeSeries(tf))
                    self.dfs["goes"] = pd.concat(
                        [self.dfs["goes"], self.goes[-1].to_dataframe()]
                    )
            if len(self.dfs["goes"]) > 0:
                self.dfs["goes"].index.name = "time"
                self.dfs["goes"] = self.dfs["goes"].reset_index()
                self.dfs["goes"] = self.dfs["goes"][
                    (self.dfs["goes"].time >= self.dates[0])
                    & (self.dfs["goes"].time <= self.dates[1])
                ]
            # Retrieve HEKTable from the Fido result and then load
            hek_results = result["hek"]
            if len(hek_results) > 0:
                self.flare = hek_results[
                    "event_starttime",
                    "event_peaktime",
                    "event_endtime",
                    "fl_goescls",
                    "ar_noaanum",
                ]
        self.dfs["goes"].drop_duplicates(subset="time", inplace=True)
        return
    
    def create_stackplots(self, fname):
        dates = [dt.datetime(2017, 5, 27, 12), dt.datetime(2017, 5, 28)]
        sp = StackPlots(
            4, 1, datetime=True
        )
        _, ax = sp.plot_stack_plots(
            np.array(self.dfs["goes"].time), self.dfs["goes"].xrsa, color="b", label=r"$\lambda$ (0.05-0.4 nm)",
            xlabel="", ylabel="", ylim=[1e-7, 1e-3], xlim=dates,
            title="Geospace Condition / 27 May 2017"
        )
        sp.plot_stack_plots(
            np.array(self.dfs["goes"].time), self.dfs["goes"].xrsb, color="r", label=r"$\lambda$ (0.1-0.8 nm)",
            xlabel="", ylabel=r"GOES X-rays, $W/m^2/nm$", ylim=[1e-9, 1e-3], ax=ax,
            xlim=dates,
        )
        ax.set_yscale("log")
        ax.legend(loc=1)
        _, ax = sp.plot_stack_plots(
            np.array(self.omni.time), self.omni.BzGSE, color="b",
            xlabel="", ylabel=r"IMF ($B_i$), nT", ylim=[-20, 20], label=r"$B_z$",
            xlim=dates,
        )
        sp.plot_stack_plots(
            np.array(self.omni.time), self.omni.ByGSE, color="k",
            xlabel="", ylabel="", ylim=[-20, 20], label=r"$B_y$",
            ax=ax, xlim=dates,
        )
        ax.legend(loc=1)
        _, ax = sp.plot_stack_plots(
            np.array(self.omni.time), self.omni.FlowSpeed, color="k",
            xlabel="", ylabel="SW Speed, km/s", ylim=[250, 500], xlim=dates,
        )
        sp.plot_stack_plots(
            np.array(self.omni.time), self.omni.ProtonDensity, color="r", 
            xlabel="", ylabel="Proton Density, /cc", ylim=[0, 100], ax=ax.twinx(), ylabel_color="r",
            xlim=dates,
        )
        _, ax = sp.plot_stack_plots(
            np.array(self.omni.time), self.omni.AsyH, color="k", 
            xlabel="Time, UT", ylabel="AsyH, nT", ylim=[], xlim=dates,
        )
        sp.plot_stack_plots(
            np.array(self.omni.time), self.omni.AE, color="r", 
            xlabel="Time, UT", ylabel="AE, nT", ylim=[], ax=ax.twinx(), ylabel_color="r",
            xlim=dates,
        )
        print(self.omni.columns)
        sp.save_fig(fname)
        sp.close()
        return