import logging


def plot_freq(ax, times, freq, xlim=None, xticks=None):
    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator

    # Format the yaxis.
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction="out")
    ax.set_ylim(bottom=8, top=20)
    ax.yaxis.set_minor_locator(MultipleLocator())
    ax.yaxis.set_tick_params(direction="out", which="minor")

    # Plot the TX frequency.
    ax.plot_date(
        date2num(times), freq, fmt="k-", tz=None, xdate=True, ydate=False, markersize=2
    )

    # Format the xaxis.
    if xlim is not None:
        ax.set_xlim(xlim)
    if xticks is not None:
        ax.set_xticks(xticks)

    # Add labels to identify the frequency axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]
    fig.text(pos[0] - 0.01, pos[1] + 0.005, "10", ha="right", va="bottom", size=8)
    fig.text(pos[0] - 0.01, pos[1] + pos[3] - 0.015, "16", ha="right", va="top", size=8)
    fig.text(
        pos[0] - 0.07,
        pos[1] + pos[3] / 2.0,
        "Freq",
        ha="center",
        va="center",
        size=9,
        rotation="vertical",
    )
    fig.text(
        pos[0] - 0.05,
        pos[1] + pos[3] / 2.0,
        "[MHz]",
        ha="center",
        va="center",
        size=7,
        rotation="vertical",
    )
    l = Line2D(
        [pos[0] - 0.04, pos[0] - 0.04],
        [pos[1] + 0.01, pos[1] + pos[3] - 0.01],
        transform=fig.transFigure,
        clip_on=False,
        ls="-",
        color="k",
        lw=1.5,
    )
    ax.add_line(l)
    ax.set_xticklabels([])
    # use only 2 major yticks
    ax.set_yticks([10, 16])
    ax.set_yticklabels([" ", " "])
    return


def plot_nave(ax, times, nave, xlim=None, xticks=None, ytickside="right"):

    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator

    # Format the yaxis
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction="out")
    ax.set_ylim(bottom=0, top=80)
    ax.yaxis.set_minor_locator(MultipleLocator(base=5))
    ax.yaxis.set_tick_params(direction="out", which="minor")

    # Plot the number of averages.
    ax.plot_date(
        date2num(times), nave, fmt="k:", tz=None, xdate=True, ydate=False, markersize=2
    )

    # Format the xaxis.
    if xlim is not None:
        ax.set_xlim(xlim)
    if xticks is not None:
        ax.set_xticks(xticks)

    # Add labels to identify the nave axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]
    fig.text(
        pos[0] + pos[2] + 0.01, pos[1] - 0.004, "0", ha="left", va="bottom", size=8
    )
    fig.text(pos[0] + pos[2] + 0.01, pos[1] + pos[3], "80", ha="left", va="top", size=8)
    fig.text(
        pos[0] + pos[2] + 0.06,
        pos[1] + pos[3] / 2.0,
        "Nave",
        ha="center",
        va="center",
        size=8.5,
        rotation="vertical",
    )

    l = Line2D(
        [pos[0] + pos[2] + 0.07, pos[0] + pos[2] + 0.07],
        [pos[1] + 0.01, pos[1] + pos[3] - 0.01],
        transform=fig.transFigure,
        clip_on=False,
        ls=":",
        color="k",
        lw=1.5,
    )
    ax.add_line(l)
    ax.set_xticklabels([])
    # use only 2 major yticks
    ax.set_yticks([0, 80])
    ax.set_yticklabels([" ", " "])
    if ytickside == "right":
        ax.yaxis.tick_right()
    return


def plot_skynoise(ax, times, sky, xlim=None, xticks=None):

    import numpy as np
    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator

    # Format the yaxis.
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction="out")
    ax.set_ylim(bottom=0, top=6)
    ax.yaxis.set_minor_locator(MultipleLocator())
    ax.yaxis.set_tick_params(direction="out", which="minor")

    # Plot the sky noise data.
    ax.plot_date(
        date2num(times), np.log10(sky), fmt="k-", tz=None, xdate=True, ydate=False
    )

    # Format the xaxis.
    if xlim is not None:
        ax.set_xlim(xlim)
    if xticks is not None:
        ax.set_xticks(xticks)

    # Add labels to identify the noise axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]
    fig.text(pos[0] - 0.01, pos[1] + 0.004, r"$10^0$", ha="right", va="bottom", size=8)
    fig.text(pos[0] - 0.01, pos[1] + pos[3], r"$10^6$", ha="right", va="top", size=8)
    fig.text(
        pos[0] - 0.07,
        pos[1] + pos[3] / 2.0,
        "N.Sky",
        ha="center",
        va="center",
        size=8.5,
        rotation="vertical",
    )
    l = Line2D(
        [pos[0] - 0.06, pos[0] - 0.06],
        [pos[1] + 0.01, pos[1] + pos[3] - 0.01],
        transform=fig.transFigure,
        clip_on=False,
        ls="-",
        color="k",
        lw=1.5,
    )
    ax.add_line(l)
    ax.set_xticklabels([])
    # Only use 2 major yticks.
    ax.set_yticks([0, 6])
    ax.set_yticklabels([" ", " "])
    return


def plot_searchnoise(ax, times, search, xlim=None, xticks=None, ytickside="right"):

    import numpy as np
    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator

    # Format the yaxis.
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction="out")
    ax.set_ylim(bottom=0, top=6)
    ax.yaxis.set_minor_locator(MultipleLocator())
    ax.yaxis.set_tick_params(direction="out", which="minor")

    # Plot the search noise data.
    ax.plot_date(
        date2num(times),
        np.log10(search),
        fmt="k:",
        tz=None,
        xdate=True,
        ydate=False,
        lw=1.5,
    )

    # Format the xaxis.
    if xlim is not None:
        ax.set_xlim(xlim)
    if xticks is not None:
        ax.set_xticks(xticks)

    # Add labels to identify the noise axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]

    fig.text(
        pos[0] + pos[2] + 0.01,
        pos[1] + 0.004,
        r"$10^0$",
        ha="left",
        va="bottom",
        size=8,
    )
    fig.text(
        pos[0] + pos[2] + 0.01, pos[1] + pos[3], r"$10^6$", ha="left", va="top", size=8
    )
    fig.text(
        pos[0] + pos[2] + 0.06,
        pos[1] + pos[3] / 2.0,
        "N.Sch",
        ha="center",
        va="center",
        size=8.5,
        rotation="vertical",
    )

    l = Line2D(
        [pos[0] + pos[2] + 0.07, pos[0] + pos[2] + 0.07],
        [pos[1] + 0.01, pos[1] + pos[3] - 0.01],
        transform=fig.transFigure,
        clip_on=False,
        ls=":",
        color="k",
        lw=1.5,
    )
    ax.add_line(l)
    ax.set_xticklabels([])
    # use only 2 major yticks
    ax.set_yticks([0, 6])
    ax.set_yticklabels([" ", " "])
    if ytickside == "right":
        ax.yaxis.tick_right()

    return
