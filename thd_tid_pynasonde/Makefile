# CONUS DIDBASE download / scaling / plotting pipeline.
# Targets are added incrementally alongside the scripts they wrap.

CONDA_ENV := apep
# `conda run -n ... python` buffers the whole child process's stdout/stderr
# internally and only flushes it once the subprocess exits, which silently
# swallows tqdm's live progress bar (it all appears at once, at the end).
# Invoking the env's python binary directly avoids that wrapper entirely.
PY        := /home/chakras4/anaconda3/envs/$(CONDA_ENV)/bin/python -u
DATA_DIR  := Data
# Defaults for the analysis / event-day targets; override on the command line, e.g.
#   make mstid STATIONS=AL945,MHJ45
#   make skymap-day STATION=AL945 DATE=2024-04-02
STATIONS  ?= AL945
# skymap-day takes one station: STATION=..., or the first entry of STATIONS=...
comma     := ,
STATION   ?= $(word 1,$(subst $(comma), ,$(STATIONS)))
DAYS_CSV  ?= Data/analysis/mstid_days.csv

.DEFAULT_GOAL := help

.PHONY: help download download-raw download-sky quality plot-ts quiet-curve mstid skymap-day analysis clean clean-pyc clean-data test

help:
	@echo "Targets:"
	@echo "  make download      - scaled DIDBASE params, all CONUS stations, config.toml years/months"
	@echo "  make download-raw  - raw RSF/MMM ionograms (feeds plot_viterbi.py)"
	@echo "  make download-sky  - DDM skymaps: one event day (STATIONS=MHJ45 DATE=2023-04-24) or all enhanced days in DAYS_CSV=$(DAYS_CSV)"
	@echo "  make quality       - QC scaled params -> Data/quality/day_quality.csv (pass/fail + reason)"
	@echo "  make plot-ts       - QC-passed foF2/hmF2/foE/hE + band-pass per station per DAY -> Data/figures/giro_scaled/<station>/ (all stations; STATIONS=AU930,AL945 to restrict)"
	@echo "  make quiet-curve   - Kp/SYM-H/AE quiet days + monthly median curves (STATIONS=$(STATIONS)) -> Data/analysis/, Data/figures/quiet_curve/"
	@echo "  make mstid         - per-day MSTID signatures vs quiet baseline (STATIONS=$(STATIONS)) -> Data/analysis/mstid_*.csv, Data/figures/mstid/"
	@echo "  make skymap-day    - skymap tilt/wavelet figure set for one event day (STATION=$(STATION) DATE=YYYY-MM-DD) -> Data/figures/skymap/<station>/<date>/"
	@echo "  make analysis      - quality -> quiet-curve -> mstid in order (STATIONS=$(STATIONS))"
	@echo "  make test          - offline unit tests (tests/)"
	@echo "  make clean         - remove __pycache__ / .pyc (safe, no downloaded data)"
	@echo "  make clean-data    - DESTRUCTIVE: remove $(DATA_DIR)/ (all downloads + plots)"

download:
	$(PY) didbase_pipeline.py

download-raw:
	$(PY) didbase_pipeline.py --raw

# DATE=YYYY-MM-DD (or several, space-separated and quoted) downloads just those event days;
# without DATE, every enhanced day in DAYS_CSV for STATIONS is fetched.
download-sky:
	$(PY) didbase_pipeline.py --sky --stations $(subst $(comma), ,$(STATIONS)) $(if $(DATE),--dates $(DATE),--days-csv $(DAYS_CSV))

quality:
	$(PY) quality.py

# STATIONS: comma-separated URSI list for the analysis scripts (e.g. STATIONS=AL945,MHJ45)
quiet-curve: quality
	ANALYSIS_STATIONS=$(STATIONS) $(PY) quiet_curve.py

mstid: quiet-curve
	ANALYSIS_STATIONS=$(STATIONS) $(PY) mstid.py

analysis: mstid

skymap-day:
	@test -n "$(DATE)" || (echo "usage: make skymap-day STATION=AL945 DATE=YYYY-MM-DD   (several consecutive days: DATE=\"2026-04-14 2026-04-15\")" && exit 1)
	$(PY) skymap.py $(STATION) $(DATE)

# STATIONS restricts plot-ts only when given explicitly (command line or environment);
# otherwise every station is plotted. Dates: RSF_TSPLOT_DAILY_START / RSF_TSPLOT_DAILY_END (YYYY-MM-DD).
PLOT_STATIONS := $(if $(filter-out file default,$(origin STATIONS)),$(STATIONS),)

plot-ts: quality
	RSF_TSPLOT_DAILY=1 $(if $(PLOT_STATIONS),RSF_TSPLOT_DAILY_STATIONS=$(PLOT_STATIONS)) $(PY) plot_param_tsplot.py

test:
	$(PY) -B -m unittest discover -s tests -v

clean: clean-pyc

clean-pyc:
	find . -name '__pycache__' -type d -exec rm -rf {} +
	find . -name '*.pyc' -delete

clean-data:
	rm -rf $(DATA_DIR)
