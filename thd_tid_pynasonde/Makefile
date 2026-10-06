# CONUS DIDBASE download / scaling / plotting pipeline.
# Targets are added incrementally alongside the scripts they wrap.

CONDA_ENV := apep
# `conda run -n ... python` buffers the whole child process's stdout/stderr
# internally and only flushes it once the subprocess exits, which silently
# swallows tqdm's live progress bar (it all appears at once, at the end).
# Invoking the env's python binary directly avoids that wrapper entirely.
PY        := /home/chakras4/anaconda3/envs/$(CONDA_ENV)/bin/python -u
DATA_DIR  := Data

.DEFAULT_GOAL := help

.PHONY: help download download-raw quality plot-ts clean clean-pyc clean-data test

help:
	@echo "Targets:"
	@echo "  make download      - scaled DIDBASE params, all CONUS stations, config.toml years/months"
	@echo "  make download-raw  - raw RSF/MMM ionograms (feeds plot_viterbi.py)"
	@echo "  make quality       - QC scaled params -> Data/quality/day_quality.csv (pass/fail + reason)"
	@echo "  make plot-ts       - QC-passed foF2/hmF2/foE/hE + 15-60 min band-pass, per station per DAY -> Data/figures/giro_scaled/<station>/"
	@echo "  make test          - offline unit tests (tests/)"
	@echo "  make clean         - remove __pycache__ / .pyc (safe, no downloaded data)"
	@echo "  make clean-data    - DESTRUCTIVE: remove $(DATA_DIR)/ (all downloads + plots)"

download:
	$(PY) didbase_pipeline.py

download-raw:
	$(PY) didbase_pipeline.py --raw

quality:
	$(PY) quality.py

plot-ts: quality
	RSF_TSPLOT_DAILY=1 $(PY) plot_param_tsplot.py

test:
	$(PY) -B -m unittest discover -s tests -v

clean: clean-pyc

clean-pyc:
	find . -name '__pycache__' -type d -exec rm -rf {} +
	find . -name '*.pyc' -delete

clean-data:
	rm -rf $(DATA_DIR)
