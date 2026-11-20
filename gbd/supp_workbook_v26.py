"""v2.6 supplementary workbook (2nd review): single file out/main/GBD2023_Supplementary_Tables_S2-S5.xlsx.
S2 from NCI APC Web Tool exports (apc_summary_nci2023.csv), df15 only; S5 from calendar-year paired contrasts."""
import pandas as pd, numpy as np, shutil
B=''; F=B+'out/main/GBD2023_Supplementary_Tables_S2-S5.xlsx'
A=B+'out/main/archive_superseded/GBD2023_Supplementary_Tables_S2-S5_v2.5.xlsx'
old=pd.read_excel(A,sheet_name=None,engine='openpyxl')
S=pd.read_csv(B+'out/apc/apc_summary_nci2023.csv')
pf=lambda p:'<0.001' if p<0.001 else f'{p:.3f}'
S2=pd.DataFrame({'Measure':S.measure,'Sex':S.sex,'Etiology':S.etiology,
 'Net drift, %/year (95% CI)':[f'{a:.2f} ({b:.2f} to {c:.2f})' for a,b,c in zip(S.net_drift,S.nd_lo,S.nd_hi)],
 'P, net drift = 0':S.p_netdrift.map(pf),'Wald χ², local drifts = net drift':S.X2_localdrift.round(1),'df':S.ld_df,
 'P, local drifts = net drift':S.p_localdrift_equal.map(pf),'P, period deviations = 0 (df 4)':S.p_period_dev.map(pf),'P, cohort deviations = 0 (df 18)':S.p_cohort_dev.map(pf),
 'Reference groups (age; period; cohort)':'55-59 y; 2004-2008; 1949','Source':'NCI APC Web Tool export: '+S.source})
I=pd.read_csv(B+'out/apc/apc_summary_gbd2023.csv')
S2b=pd.DataFrame({'Measure':I.measure,'Sex':I.sex,'Etiology':I.etiology,
 'Net drift, %/year (95% CI)':[f'{a:.2f} ({b:.2f} to {c:.2f})' for a,b,c in zip(I.net_drift,I.nd_lo,I.nd_hi)],
 'P, net drift = 0':I.p_netdrift.map(pf),'P, local drifts = net drift (df 15)':I.p_localdrift_equal.map(pf),
 'P, period deviations = 0':I.p_period_dev.map(pf),'P, cohort deviations = 0':I.p_cohort_dev.map(pf),
 'Source':'In-house Python implementation (gbd/apc.py); SENSITIVITY ONLY, not the main analysis'})
c=pd.read_csv(B+'out/contrast/sex_contrasts.csv')
S5=pd.DataFrame({'Measure':c.measure,'Contrast':c.contrast,
 'AAPC difference, percentage points/year (95% CI)':[f'{a:.2f} ({b:.2f} to {d:.2f})' for a,b,d in zip(c.AAPC_diff,c.AAPC_lo,c.AAPC_hi)],'P (AAPC)':c.AAPC_p.map(pf),
 'Slope difference ×100, calendar-year paired (95% CI)':[f'{a:.2f} ({b:.2f} to {d:.2f})' for a,b,d in zip(c.slope_diff_pct,c.slope_lo,c.slope_hi)],'P (slope)':c.slope_p.map(pf)})
notes=pd.DataFrame({'Notes':[
 'Source: GBD 2023, China. All CIs are conditional on GBD point estimates and do not propagate GBD uncertainty.',
 'S2: NCI Age-Period-Cohort Web Tool (analysistools.cancer.gov/apc; Rosenberg et al. 2014), default settings; ages 20-94 (15 groups), periods 1994-1998 to 2019-2023; reference age 55-59, period 2004-2008, cohort 1949. Input: 5-year sums of annual GBD 2023 counts and person-years. Local-drift Wald test with 15 df as reported by the tool. Values exported from the tool (export file named per row); workbook-to-series mapping verified (gbd/apc_nci_mapping_check.py). These 24 NCI estimates are the main APC results reported in the manuscript.',
 'S2 in-house sensitivity (separate sheet): same inputs fitted with an in-house Python implementation of the APC model (gbd/apc.py), used only for the age- and period-range sensitivity analyses (Supplementary Figure S2). Its estimates differ from the NCI tool in estimation and CI method (e.g. male MASH incidence net drift -0.22% [-0.56 to 0.12] in-house vs -0.20% [-0.38 to -0.01] NCI) and must not be read as the main results.',
 'S3: EAPC from ln(ASR)~year OLS, 1990-2023, 95% CI from t distribution (32 df); GBD 2023 age-standardised rates, China.',
 'S4: Das Gupta three-factor decomposition, ages 15+; counts summed over 5-year age groups and may differ slightly from GBD all-age totals. "All five etiologies" = HBV + HCV + alcohol + MASH + other causes (hepatoblastoma excluded). Each column is rounded independently from full-precision results, so the net change may differ by 1 from the difference of the rounded counts.',
 'S5 AAPC difference: difference of NCI Joinpoint AAPCs (percentage points per year); SE of each AAPC approximated as (upper - lower)/3.92 from its empirical-quantile CI, AAPCs assumed independent (approximate).',
 'S5 slope difference: for each calendar year t (1990-2023), D_t = sum of w_i x ln(ASR_i,t) for the contrast weights; OLS of D_t on year (34 observations); slope x100 (approximately percentage points per year); Newey-West (Bartlett kernel, lag 3) standard error, normal approximation, no small-sample correction. Pairing by calendar year accounts for year-to-year variation shared by the series in a contrast.',
 'S5 contrasts: overall and four etiologies (HBV, HCV, alcohol, MASH); other causes and hepatoblastoma were not contrasted.',
 'Code: gbd/sex_contrasts.py (S5), gbd/apc_nci_import.py (S2); to be released with the article (repository URL to be added).','S2 Wald tests as exported by the NCI tool: NetDrift = 0 (df 1); All Local Drifts = Net Drift (df 15); All Period Deviations = 0 (df 4); All Cohort Deviations = 0 (df 18).']})
with pd.ExcelWriter(F) as w:
    S2.to_excel(w,sheet_name='Table S2 APC',index=False); S2b.to_excel(w,sheet_name='S2 in-house sensitivity',index=False); old['Table S3 EAPC'].to_excel(w,sheet_name='Table S3 EAPC',index=False)
    old['Table S4 Decomposition'].rename(columns={'Epidemiological change':'Rate change (change in age-specific rates)'}).to_excel(w,sheet_name='Table S4 Decomposition',index=False); S5.to_excel(w,sheet_name='Table S5 Sex contrasts',index=False); notes.to_excel(w,sheet_name='Notes',index=False)
print(S2[(S2.Sex=='Male')&(S2.Etiology=='MASH')].to_string()); print(S5.to_string())
