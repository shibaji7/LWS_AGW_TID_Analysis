# CONUS DIDBASE download, RSF scaling, and diagnostic plots

All defaults (years, stations, SZA regimes, station/year scaler overrides,
plot windows, output paths) live in [`config.toml`](config.toml). Every
script also still accepts the environment-variable overrides described
below. Use the `apep` conda environment, which has the DIDBASE connector.

Prefer `make download` etc. (see [Makefile](#makefile)) or invoke the env's
python binary directly, **not** `conda run -n apep python ...`: `conda run`
buffers the whole child process's output internally and only flushes it once
the subprocess exits, which silently swallows tqdm's live progress bar (it
all appears at once, at the very end, instead of updating in place).

```bash
/home/chakras4/anaconda3/envs/apep/bin/python -u didbase_pipeline.py --list-stations
/home/chakras4/anaconda3/envs/apep/bin/python -u didbase_pipeline.py
```

## 1. Download: `didbase_pipeline.py`

Downloads every day from April 1 through September 30 in 2024, 2025, and
2026 (capped at the current UTC time) for every CONUS Digisonde station.

```text
Data/didbase_giro_scaled/<year>/<station>/<YYYY-MM-DD>.csv          scaled params
Data/didbase_giro_scaled/<year>/<station>/<YYYY-MM-DD>.json         completion record
```

Each CSV contains native-cadence observations, all catalogued scaled
parameters, quality flags, scaler ID, equipment ID, and program ID.
DIDBASE's connector chooses one internally consistent scaler result per
measurement, preferring human-reviewed results when available. No median,
quality rejection, or local ionogram scaling is applied here.

Stations are discovered by matching DIDBASE URSI codes against the country
and coordinates in the [Digisonde instrument catalogue](https://www.digisonde.com/stationlist.php).
Only US locations within the contiguous-US coordinate bounds
(`config.toml`'s `[stations]`) are selected. Stations missing from that
catalogue (e.g. `WS833`) are supplemented from `[[stations.extra]]`.
Supply `--station-catalog locations.csv` (columns `URSI,Country,Lat,Long`)
to extend/correct geography.

Historical files are reused only when their checksum, parameter list, and
completion record match. Today's partial data are always queried again.
Use `--overwrite` to retrieve late arrivals or revised historical
observations. Failed requests are retried (`--retries`, default 3), reported
in `failures.csv`, and cause a nonzero exit status; rerunning resumes.

Optional flags: `--years 2024`, `--stations AU930 BC840`, `--retries 5`,
`--output-dir PATH`. Connections use `DIDBASE_USER`, `DIDBASE_PASSWORD`,
`DIDBASE_URL`, and `DIDBASE_ROLE`.

### Raw ionograms: `--raw`

```bash
/home/chakras4/anaconda3/envs/apep/bin/python -u didbase_pipeline.py --raw
```

Fetches the underlying ionogram files (RSF, falling back to MMM) that
`plot_viterbi.py` needs to draw traces -- DIDBASE's scaled-parameter query
alone cannot supply these. Same station/date selection, caching, and retry
behavior as the scaled download.

```text
Data/didbase_giro_scaled/<year>/<station>/raw/<YYYY-MM-DD>/*.RSF
Data/didbase_giro_scaled/<year>/<station>/raw/<YYYY-MM-DD>.json
```

## 2. Local RSF Viterbi scaling: `rsf_batch_scaler.py`

Batch-scales raw RSF files (from `--raw` above) with the generic Digisonde
Viterbi scaler, driven by SZA regime (day/twilight/night, `config.toml`
`[sza_regimes]`) and any per-station/year override
(`config.toml` `[rsf_scaler.station_overrides.<STATION>.<YEAR>]`, e.g.
amplitude floors, RFI filtering, masked frequency bands, trailing-gap
trimming/bridging). Outputs `Data/processing/scaled/rsf_scaled_<date>_<station>.csv`.
Loops over `Data/didbase_giro_scaled/<year>/<station>/raw/<date>/` by
default; set `RSF_SCALER_GLOB` to scale an arbitrary directory instead.

## 3. Diagnostic plots

- **`plot_viterbi.py`** -- raw/cleaned/Viterbi-trace overlay figures for any
  UT window (`RSF_VITERBI_UT_START`/`RSF_VITERBI_UT_END`, station, date
  range). Reuses `rsf_batch_scaler.py`'s scaling logic directly, so
  diagnostics always match the batch-scaled output.
- **`plot_param_vs_sza.py`** -- any scaled parameter (`RSF_PARAM_SZA_PARAM`, default
  foE) versus solar zenith angle, colored by UT, reading the DIDBASE CSVs
  from step 1 directly (no local scaling needed).
- **`plot_param_tsplot.py`** -- full native-cadence time series of one scaled
  parameter (`RSF_TSPLOT_PARAM`, default foF2) for one station over any date
  range, also reading the DIDBASE CSVs directly. `make plot-ts` (or
  `RSF_TSPLOT_DAILY=1`) instead batch-plots **foF2 and hmF2 for every
  downloaded station, one figure per day**, into
  `Data/figures/giro_scaled/<station>/`; `RSF_TSPLOT_ALL=1` batch-plots one
  figure per station spanning the whole season instead.

## Shared SZA/scaling logic: `utils.py`

Solar zenith angle calculation, SZA regime classification, and the
station/year override, amplitude-prefilter, gap-trimming, and reliability-gate
helpers used by both `rsf_batch_scaler.py` and `plot_viterbi.py`. All numeric
settings come from `config.toml`.

## Tests

```bash
/home/chakras4/anaconda3/envs/apep/bin/python -m unittest discover -s tests -v
```

## Makefile

```bash
make help          # list targets
make download       # scaled DIDBASE params
make download-raw   # raw RSF/MMM ionograms
make plot-ts        # foF2/hmF2 TS, one figure per station per DAY -> Data/figures/giro_scaled/<station>/
make test           # offline unit tests
make clean          # remove __pycache__ / .pyc (safe)
make clean-data     # DESTRUCTIVE: remove Data/
```
