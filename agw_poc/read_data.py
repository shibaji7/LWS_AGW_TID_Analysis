import h5py as hf
import datetime as dt
import numpy as np
import pandas as pd

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

def round_time(dx: dt.datetime, unit=dt.timedelta(seconds=30)):
    seconds = (dx - dt.datetime.min).total_seconds()
    unit_seconds = unit.total_seconds()
    half_over = seconds + unit_seconds / 2
    rounded_seconds = half_over - half_over % unit_seconds
    return dt.datetime.min + dt.timedelta(seconds=rounded_seconds)

def read_tec_paul(ftec,filt=5,timed=dt.datetime(2017,5,27,17,30,)):
    '''
    Parameters
    ----------
    ftec : full path to .mat file containing fulltimemat type scatter data
        .mat file.
    filt : 4/5/6/7
        4: acoustics. 10 min high pass
        5: 5-40min filter
        6: 60min high pass filter
        7: 120min highpass filter
        Default is 5
    Returns
    -------
    timedt : datetimes of recorded times
    tecs : vertical TECs at every timestamp
    ipplats : Latitude of observations at every timestamp
    ipplons : Longitude of observations at every timestamp
    '''
    alltec=hf.File(ftec,'r')
    timestamps=alltec['UTT']
    times=np.concatenate(timestamps,axis=0)
    timedt=np.array([])
    ipplats=[]
    ipplons=[]
    tecs=[]
    # arrs=alltec['fulltimedata']
    for tidx,ts in enumerate(times):
        tsdt = round_time(
            dt.datetime.fromordinal(int(ts)) + dt.timedelta(days=ts%1)-dt.timedelta(days=366)
        )
        ############################
        if timed is not None:
            if timed==tsdt:
                arrs=alltec['fulltimedata']
                ref=arrs[tidx][0]
                dset=np.array(alltec[ref])
                try:
                    ntimes,nstat=np.shape(dset)
                except:
                    continue
                ipplats.append(dset[2])
                ipplons.append(dset[3])
                tecs.append(dset[filt])
                timedt=np.append(timedt,tsdt)
        else:
            arrs=alltec['fulltimedata']
            ref=arrs[tidx][0]
            dset=np.array(alltec[ref])
            try:
                ntimes,nstat=np.shape(dset)
            except:
                continue
            #################
            ipplats.append(dset[2])
            ipplons.append(dset[3])
            tecs.append(dset[filt])
            timedt=np.append(timedt,tsdt)
            #################
            # if tidx==100: break
    return timedt,tecs,ipplats,ipplons

def get_gridded_parameters(
    q, xparam="lon", yparam="lat", zparam="tec", r=1, rounding=True
):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[[xparam, yparam, zparam]]
    if rounding:
        plotParamDF.loc[:, xparam] = np.round(plotParamDF[xparam].tolist(), r)
        plotParamDF.loc[:, yparam] = np.round(plotParamDF[yparam].tolist(), r)
    plotParamDF = plotParamDF.groupby([xparam, yparam]).mean().reset_index()
    plotParamDF = plotParamDF[[xparam, yparam, zparam]].pivot(index=xparam, columns=yparam)
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y = np.meshgrid(x, y)
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
        np.isnan(plotParamDF[zparam].values), plotParamDF[zparam].values
    )
    return X, Y, Z

def plot_datetime(timedt,tecs,ipplats,ipplons):
    import matplotlib.pyplot as plt
    time, tecs,ipplats,ipplons = timedt[0], tecs[0],ipplats[0],ipplons[0]
    o = pd.DataFrame()
    o["lat"], o["lon"], o["tec"] = ipplats, ipplons, tecs
    X, Y, Z = get_gridded_parameters(o)
    plt.pcolormesh(Y, X, Z.T, cmap="plasma")
    plt.ylim(20, 60)
    plt.xlim(-150, -60)
    plt.colorbar(label="dTEC")
    plt.clim(-.5,.5)
    plt.xlabel("Lon")
    plt.ylabel("Lat")
    plt.title(f"UT: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    plt.savefig("agw_poc/example_out.png")
    return

if __name__ == "__main__":
    # timedt,tecs,ipplats,ipplons = read_tec_paul("fulltimedata_05_27_2017-GRE.mat")
    # plot_datetime(timedt,tecs,ipplats,ipplons)
    timedt,tecs,ipplats,ipplons=read_tec_paul("fulltimedata_05_27_2017-GRE.mat", timed=None)
    uthr=17.5
    idx=int((uthr%24)*60*2)
    maplons=ipplons[idx]
    maplats=ipplats[idx]
    tecs=tecs[idx]
    tstamp=timedt[idx]
    print(maplons, maplats, tecs)
    # elvals=el[idx]
    # mskel=elvals>30
    ### test plot data
    fig, ax = plt.subplots(figsize=(15, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    # Add map features
    ax.set_extent([-130, -60, 20, 55], crs=ccrs.PlateCarree())  # U.S. bounds
    ax.add_feature(cfeature.BORDERS, linewidth=1)
    # ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)
    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, alpha=0.7)
    gl.top_labels = False  # Remove top labels
    gl.right_labels = False  # Remove right labels
    # sc = ax.scatter(maplons[mskel], maplats[mskel], c=tecs[mskel], cmap='RdBu_r',
    #                 s=25, vmin=-0.15, vmax=0.15)
    sc = ax.scatter(maplats, maplons, c=tecs, cmap='RdBu_r',
                    s=25, vmin=-0.15, vmax=0.15, transform=ccrs.PlateCarree(), alpha=0.5)
    # Add colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal',
                        pad=0.07, aspect=40, shrink=0.8, label='dTEC (TECu)')
    ax.set_title('%s'%tstamp.strftime('%D %H:%M:%S'))
    fig.savefig("agw_poc/example_out.png")