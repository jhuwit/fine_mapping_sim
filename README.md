---
editor_options: 
  markdown: 
    wrap: 72
---

This repository accompanies the manuscript ''Fine-Mapping the
Association of Acute Kidney Injury with Mean Arterial and Central Venous
Pressures during Coronary Artery Bypass Surgery'' in Anesthesia &
Analgesia.

It will generate data to mimic hemodynamic time series, then apply the
data processing steps used in our analysis and generate some (fake)
results.

# code

`00_generate_data.R` will simulate data in the same structure as the
hemodynamics data we would have received from `n` patients. Outputs: -
`data/raw/hemodynamics/individual_timeseries/*.csv` -
`data/raw/covariates.csv` - `data/raw/event_start_stop_times.csv`

`01_process_ts.R` processes individual time series data Inputs: -
`data/raw/hemodynamics/individual_timeseries/*.csv` -
`data/raw/event_start_stop_times.csv` Outputs: -
`data/processed/hemo_timeseries_all.rds` -
`data/processed/hemo_timeseries_cohort.rds`

`02_interpolate.R` interpolate missing time series data Input: -
`data/processed/hemo_timeseries_cohort.rds` Output: -
`data/processed/hemo_timeseries_interp_cohort.rds`

## Interpolation Details

```         
-   Linearly interpolate periods of missingness $\leq 15$ minutes for both MAP and CVP

    -   Interpolation start: average of 3 minutes prior to missingness

    -   Interpolation end: average of 3 minutes after missingness

    -   Only interpolate if period of missingness is less than or equal to 15 minutes, otherwise leave values as `NA`
```

`03_apply_exclusion_criteria.R`: exclude subjects based on missingness
and exclusion criteria Inputs: -
`data/processed/hemo_timeseries_interp_cohort.rds` -
`data/raw/covariates.csv` - `data/raw/event_start_stop_times.csv`
Outputs: - `data/analytic/hemo_timeseries_interp_post_exclusion.rds` -
`data/analytic/covariates_post_exclusion.csv` -
`data/analytic/exclusion_summary_cohort.csv`

## analysis

-   `utilities.R` : helpful functions for generating time in range
    variables and running regression. Source before running any other
    analysis code

-   `01_generate_cabg_manuscript_results_as_df.R`

    ```         
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
    ```

    Note: some models may not converge given the simulated nature of the
    dataset.

The rest of the steps generate all of the results included in the
manuscript; again, some of the models may not converge or make sense due
to the nature of the dataset.

# data 

### raw 

Simulated data from `00_generate_data`

### processed

Processed data from `01_process_ts`

### analytic

Data used in analysis (post-interpolation and exclusion)

# results

Results from analysis
