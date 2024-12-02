.PHONY: all

all: data/raw/event_start_stop_times.csv data/raw/covariates.csv data/raw/hemodynamics/individual_timeseries/*.csv \
data/processed/hemo_timeseries_all.rds data/processed/hemo_timeseries_cohort.rds \
data/processed/hemo_timeseries_interp_cohort.rds.rds \
data/analytic/exclusion_summary_cohort.csv data/analytic/covariates_post_exclusion.csv \
data/analytic/hemo_timeseries_interp_post_exclusion.rds

data/raw/hemodynamics/individual_timeseries/*.csv data/raw/event_start_stop_times.csv data/raw/covariates.csv: code/00_generate_data.R
	R --no-save < code/00_generate_data.R	

data/processed/hemo_timeseries_all.rds data/processed/hemo_timeseries_cohort.rds: data/raw/hemodynamics/individual_timeseries/*.csv data/raw/event_start_stop_times.csv
	R --no-save < code/01_process_ts.R

data/processed/hemo_timeseries_interp_cohort.rds.rds: data/processed/hemo_timeseries_cohort.rds
	R --no-save < code/02_interpolate.R

data/analytic/exclusion_summary_cohort.csv data/analytic/covariates_post_exclusion.csv data/analytic/hemo_timeseries_interp_post_exclusion.rds: data/processed/hemo_timeseries_interp_cohort.rds.rds data/raw/event_start_stop_times.csv data/raw/covariates.csv
	R --no-save < code/03_apply_exclusion_criteria.R


.PHONY: rebuild
rebuild:
	$(MAKE) clean
	$(MAKE) all	

.PHONY: clean
clean:
	rm -rf data/raw/*	
	rm -rf data/analytic/*	
	rm -rf data/processed/*		