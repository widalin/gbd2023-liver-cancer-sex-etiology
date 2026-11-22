# Supplementary Code — sex differences in liver cancer trends by etiology in China, 1990–2023 (GBD 2023)

Code for the manuscript "Sex differences in liver cancer incidence and mortality trends by etiology in China, 1990–2023: a GBD 2023 analysis" (submission version, November 2026).

GBD estimates are **not** redistributed here. Download them from the GBD Results Tool (https://vizhub.healthdata.org/gbd-results/) under IHME's terms, using the queries below.

## Requirements

Python 3.12 with the packages pinned in `requirements.txt` (tested 2026-11-21 in a clean virtual environment: `pip install -r requirements.txt`; `python-calamine` is the Excel engine used by `apc_nci_import.py` and `apc_nci_mapping_check.py`). Paths are relative to the package root; to run from elsewhere, set the environment variable `GBD_ROOT` to the package root (read by `consistency_2021_2023.py` and `robust_ratio_contrast.py`). Joinpoint models were fitted in NCI Joinpoint Regression Program 6.1.0, and the main age–period–cohort models in the NCI Age-Period-Cohort Web Tool (analysistools.cancer.gov/apc). Both are external programs; this code reads their exported results.

Run every script from the repository root (`python3 gbd/<script>.py`). Outputs are written to `out/`.

## Input data (place under `gbd/`)

| File | GBD Results Tool query |
|---|---|
| `gbd/GBD2023_liver_cancer_etiology_Global_China_1990-2023.csv` | GBD 2023; Global and China; liver cancer and its etiologic subcauses (HBV, HCV, alcohol, NASH, other, hepatoblastoma) and all causes; incidence, prevalence, deaths, DALYs; male, female, both; all ages (number, rate) and age-standardised rate; 1990–2023 |
| `gbd/raw/coauthor_A_china_age_etiology_2023/*.csv` | GBD 2023; China; the five etiologies; incidence and deaths; male and female; 5-year age groups 15–19 to 95+; number and rate; 1990–2023 |
| `gbd/raw/gbd2023_china_pop/pop_derived_allcause.csv` | Age- and sex-specific population derived as all-cause death counts ÷ rates (GBD 2023, China, same ages and years) |
| `gbd/raw/jkus_*` | GBD 2023; Japan, Republic of Korea, USA; liver cancer etiologies; incidence and deaths; male, female, both; all-age number and age-standardised rate; 1990–2023 |
| `gbd/raw/wpp2024/wpp2024_china_pop5_medium_1950_2100.csv` | UN World Population Prospects 2024, medium variant, China, population by sex and 5-year age group (population.un.org/wpp); used only by `pop_denominator_check.py` |
| `gbd/raw/gbd_std_pop_2019.csv` | GBD world standard population, percentage weights by age group (**included in this package**); used by `robust_female_hcv.py`, `equal_weight_asr.py`, `apc_sensitivity.py` and `consistency_2021_2023.py` |
| `gbd/raw/gbd2021_china_age/`, `gbd/raw/gbd2021_china_pop/` | GBD 2021 Results Tool; China; male and female; incidence and deaths by etiology and population; 5-year ages 15–19 to 95+; 1990–2021 (used by `consistency_2021_2023.py` and `robust_ratio_contrast.py`) |

File names inside `gbd/raw/` must match the paths in the scripts; edit the path constants at the top of each script if yours differ.

## External model outputs

- **NCI Joinpoint 6.1.0**: export the results of the 270 main-analysis series, plus the two sensitivity fits of 24 series (run A, 1990–2019; run B, 1990–2023 with at most three joinpoints; 48 sensitivity fits in total), and place them in `out/nci_joinpoint/rerun_output/02_raw_exports/`, `out/nci_joinpoint/rerun2_output/` and `out/nci_joinpoint/rerun3_output/` (read by `gbd/nci.py`). Settings: log-linear model, uncorrelated errors, heteroscedastic weighting by SE = (upper − lower UI)/3.92, 0–5 joinpoints, minimum two observations from either end and between joinpoints, weighted BIC, empirical-quantile CIs (5001 resamples).
- **NCI Age-Period-Cohort Web Tool** (default settings; ages 20–94 in 15 groups; periods 1994–1998 to 2019–2023; counts and person-years are 5-year sums of annual GBD estimates): save the 24 exported workbooks to `data/nci_apc_exports/` under their original export names. Copies are supplied in the Supplementary Data folder `NCI_APC_exports_S3/` (renamed by sex, etiology and measure); the original name of each is given in the Source column of Supplementary Table S3.

## Pipeline (in order)

1. `table1.py` — Table 1 (counts, ASRs, shares, 1990 vs 2023).
2. `main_figs.py` — integrated long dataset, Table 2 (NCI AAPCs via `nci.py`), Figures 1–3 and 5 (Figure 5 was Figure 4 in earlier drafts). `jp2.py` is an in-house joinpoint fallback used only if an NCI series is missing (none in the final version).
3. `supp_figs.py`, `fig2_v25.py` — supplementary figures and final Figure 2.
4. `apc_nci_import.py`, `apc_nci_mapping_check.py` — import NCI Age-Period-Cohort Web Tool exports and check workbook-to-series identity; `apc.py`, `apc_sensitivity.py` — in-house age–period–cohort sensitivity analyses (age and period ranges) and Supplementary Figure S2 (cohort rate ratios re-referenced to the 1953.5 birth cohort); `fig5_v25.py`, `apc_fig.py` — Figure 6 and Supplementary Figure S10.
5. `decomp.py`, `fig6_v25.py` — Das Gupta decomposition (ages 15+) and Figure 7. Defaults reproduce the main analysis (GBD 2023 counts; derived GBD 2023 population); pass a population file as the second argument for sensitivity runs.
6. `robust_female_hcv.py`, `age_specific_mash.py`, `equal_weight_asr.py` — robustness analyses (female HCV/HBV crossover, age-specific MASH trends, age weighting).
7. `sex_contrasts.py` — Supplementary Table S2 (AAPC differences and calendar-year paired log-linear contrasts with Newey–West SE, lag 3); `robust_ratio_contrast.py` — Supplementary Table S2 part B (Newey–West lags 0–5, periods ending 2019 and 2021, GBD 2021 vs GBD 2023 at ages ≥15); `fig_ratio_sex.py` — Figure 4 (HCV/HBV ratio by sex). `pop_denominator_table.py` — Supplementary Table S6 part B (run after `pop_denominator_check.py`, which writes the WPP 2024 population file). `pop_denominator_check.py` — checks of the derived GBD 2023 population denominator (internal identity, GBD 2021, UN WPP 2024); run `decomp.py <counts> out/audit/pop_wpp2024_15plus.csv sens_wpp2024` for the WPP decomposition sensitivity (Supplementary Table S6 part C).
8. `supp_workbook_v26.py` — supplementary tables workbook, working numbering (see Table numbering note below; working S3 and S4 sheets are carried over from the previous workbook version; point `A` to your copy).
9. `consistency_2021_2023.py` — cross-round comparison of 1990 female age-standardised death rates (GBD 2021 vs GBD 2023; Supplementary Methods).
10. `jk_female_hbv_hcv_contrast.py` — female HCV − HBV mortality AAPC differences for Japan and Korea (same method as Supplementary Table S2; quoted in the Results); `r16_workbook.py` — adds sheet "FigS2 net drift sensitivity" (net drifts behind Supplementary Figure S2, from `out/robust/apc_sensitivity.xlsx`) to the supplementary tables workbook.
11. `audit_manuscript.py`, `audit_triplets.py` — check numbers in the manuscript against the source data (they read the manuscript Markdown file; edit its path).

## Data sources by analysis

- Main analysis (Table 1, Figures 1–7, Supplementary Tables S1–S6 and the NCI Joinpoint and age–period–cohort models): GBD 2023 only.
- Cross-round comparison (ages ≥15 years, 1990–2021; Supplementary Table S2 part B and Supplementary Methods): GBD 2021 and GBD 2023.
- Population-denominator sensitivity analyses (Supplementary Methods; Supplementary Table S6 parts B and C): additionally UN World Population Prospects 2024, medium variant.
- Comparison with published estimates (Supplementary Table S7): values as reported in the cited publications.

## Notes

- All confidence intervals are conditional on GBD point estimates; GBD posterior draws were not used.
- "MASH" labels GBD's NASH cause.
- A smoke test from raw inputs reproduced Table 1, the AAPC table (now Supplementary Table S4) and the sex-contrast table (now Supplementary Table S2) identically.

## License

MIT License.

## Citation

The article citation will be added on publication.


## Table numbering note
Script comments, output file names and internal sheet names use the working numbering used during analysis. Mapping to the submitted Supplementary Material: working Table S2 → S3 (age–period–cohort), S3 → S5 (EAPC), S4 → S6 (decomposition, parts B and C), S5 → S2 (sex contrasts, part B), S6 → S7, S7/"Table 2" → S4 (AAPC); working Supplementary Figures S9 and S10 → S10 and S9. Workbook `GBD2023_Supplementary_Tables_S2-S5.xlsx` written by `supp_workbook_v26.py` is supplied as `GBD2023_Supplementary_Tables_workbook.xlsx` with renumbered sheets.
