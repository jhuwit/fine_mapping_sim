# CABG 

`00_generate_data.R` will simulate data in the same structure as the hemodynamics data we would have received from patients. Creates: 
- `data/raw/hemodynamics/individual_timeseries`
- `data/processed/covariates_cohort.csv`
- `data/procesed/event_start_stop_times_cohort.csv`

# data_processing

`code/data_processing` folder

## Step One

-   Run `code/data_processing/01_process_ts.R`

    -   Inputs

        -   `data/raw/hemodynamics/individual_timeseries`

    -   Outputs:

        -   `data/processed/hemo_timeseries_all.rds` : hemodynamic time series from raw files bound together

        -   `data/processed/hemo_timeseries_cohort.rds` : hemodynamic time series filtered to just individuals eligible for analysis (i.e. individuals with IDs in `full_cohort_covariates_2023-10-31`)

        -   `data/processed/hemo_timeseries_interp_cohort.rds`: hemodynamic time series for eligible cohort after interpolation



## Step Two

-   Run `code/data_processing/02_apply_exclusion_criteria.csv`

    -   Inputs

        -   `data/processed/covariates_cohort.csv`

        -   `data/processed/hemo_timeseries_cohort.rds`

        -   `data/processed/event_start_stop_times_cohort.csv`

    -   Outputs

        -   `data/processed/exclusion_summary_cohort.csv` summary of excluded and reasons for exclusion

        -   `data/processed/covariates_post_exclusion.csv` covariate file with only individuals not excluded

        -   `data/processed/hemo_timeseries_interp_post_exclusion.rds` interpoalted hemodynamic time series with only included individuals

    ## Interpolation Details

    -   Linearly interpolate periods of missingness $\leq 15$ minutes for both MAP and CVP

        -   Interpolation start: average of 3 minutes prior to missingness

        -   Interpolation end: average of 3 minutes after missingness

        -   Only interpolate if period of missingness is less than or equal to 15 minutes, otherwise leave values as `NA`

# analysis

-   `utilities.R` : helpful functions for generating time in range variables and running regression. Source before running any other analysis code

-   `CABG` all code for generation of analyses for CABG manuscript

    -   `generate_cabg_manuscript_results_as_df.R`

        -   Saves regression results to `results/data_cabg_manuscript` for univariate and adjusted regressions for MAP (`map_regressions.rds`), CVP (`cvp_regressions.rds`), MAP+CVP (`bricks_regressions.rds`) predictors and 48h and composite AKI outcomes. Column names:

            -   `term`: hemodynamic range

            -   `estimate`: associated odds ratio for spending additional 5 minutes in range

            -   `std.error`: std. error of odds ratio

            -   `statistic`: t-statistic for testing `estimate = 1`

            -   `p.value`: pvalue associated with test (unadjusted)

            -   `lb`: lower bound of 95% CI on odds ratio scale

            -   `ub`: upper bound of 95% CI on odds ratio scale

            -   `lb_bonf`: lower bound of 95% bonferroni-adjusted CI

            -   `ub_bonf`: upper bound of 95% bonferroni-adjusted CI

            -   `fdr_p`: fdr-adjusted p-value

            -   `type`: adjusted or univariate, i.e. whether covariates are present in model

            -   `outcome`: "AKI 48" or "AKI Composite", either 48hr or composite AKI outcome

            -   `lb_cma`: lower bound of 95% correlation and multiplicity adjusted (CMA) CI

            -   `ub_cma`: upper bound of 95% CMA CI
