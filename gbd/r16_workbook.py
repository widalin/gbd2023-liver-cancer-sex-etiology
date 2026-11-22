"""2026-11-22 review16 (M2, M5): add the in-house age-period-cohort net drifts behind Supplementary Figure S2
(out/robust/apc_sensitivity.xlsx, sheet Net_drift) to the Supplementary Tables workbook as sheet
'FigS2 net drift sensitivity'; write age-period-cohort in full in the Notes sheet (APC = annual percentage change only)."""
import pandas as pd, openpyxl
from openpyxl.styles import Font
B=''
WB=B+'out/submission_supp_v3.0/GBD2023_Supplementary_Tables_workbook.xlsx'
nd=pd.read_excel(B+'out/robust/apc_sensitivity.xlsx','Net_drift')
cols=['measure','sex','etiology','Main 20-94, 1994-2023','Ages 30-84','Ages 20-84','Periods 1994-2018','Periods 1991-2020']
nd=nd[cols].sort_values(['measure','sex','etiology'])
wb=openpyxl.load_workbook(WB)
SH='FigS2 net drift sensitivity'
if SH in wb.sheetnames: del wb[SH]
ws=wb.create_sheet(SH, index=wb.sheetnames.index('S3 in-house sensitivity')+1)
ws.append(['Supplementary Figure S2 (sensitivity analyses): net drift, %/year (95% CI), in-house Python age-period-cohort implementation (gbd/apc.py; gbd/apc_sensitivity.py), GBD 2023 China. Main = ages 20-94, periods 1994-2023; ages 30-84 and 20-84 keep periods 1994-2023; periods 1994-2018 drops 2019-2023; periods 1991-2020 shifts the 5-year grouping by 3 years. ns = 95% CI includes 0. Sensitivity estimates only; the main age-period-cohort results are the NCI Web Tool estimates in Table S3.'])
ws['A1'].font=Font(italic=True)
ws.append(['Measure','Sex','Etiology','Main: ages 20-94, 1994-2023','Ages 30-84','Ages 20-84','Periods 1994-2018','Periods 1991-2020'])
for c in ws[2]: c.font=Font(bold=True)
for r in nd.itertuples(index=False): ws.append(list(r))
for col,w in zip('ABCDEFGH',[11,8,9,26,22,22,22,22]): ws.column_dimensions[col].width=w
# Notes: age-period-cohort in full (M2) + pointer to the new sheet
n=wb['Notes']; done=0
for row in n.iter_rows():
    for c in row:
        if isinstance(c.value,str):
            v=c.value.replace('These 24 NCI estimates are the main APC results reported in the manuscript.','These 24 NCI estimates are the main age-period-cohort results reported in the manuscript.')
            v=v.replace('in-house Python implementation of the APC model (gbd/apc.py)','in-house Python implementation of the age-period-cohort model (gbd/apc.py)')
            if v!=c.value: c.value=v; done+=1
n.append(["FigS2 net drift sensitivity (separate sheet): in-house net drifts for every age- and period-range specification shown in Supplementary Figure S2 (male MASH-related incidence: -0.43% to -0.20% per year; main in-house specification -0.22%). In Supplementary Figure S2, each cohort rate-ratio curve is re-referenced to the 1953.5 birth cohort (log-linear interpolation between 1950.5 and 1955.5 for periods 1991-2020). In this workbook APC denotes annual percentage change only."])
if 'Table S3 APC' in wb.sheetnames: wb['Table S3 APC'].title='Table S3 age-period-cohort'
wb.save(WB); print('notes edited',done,'sheets',wb.sheetnames)
