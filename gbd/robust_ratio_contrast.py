"""2026-11-19 third-party review: robustness of the female-minus-male HCV/HBV log-ratio slope contrast (abstract 0.42).
(a) official GBD 2023 all-age ASRs: Newey-West lags 0-5, period 1990-2019, 1990-2021.
(b) cross-round: GBD 2021 vs GBD 2023 ASRs for ages >=15 (GBD 2019 std pop renormalised), same procedure (consistency_2021_2023.py), 1990-2021."""
import pandas as pd,numpy as np,statsmodels.api as sm,glob,runpy,io,contextlib
B=''
d=pd.read_csv(B+'out/main/master_dataset_long.csv')
d=d[(d.location=='China')&(d.age=='Age-standardized')&(d.metric=='Rate')&(d.sex!='Both')]
def ser(src,m,s,e): return src[(m,s,e)]
off={}
for (m,s,e),g in d.groupby(['measure','sex','etiology']): off[(m,s,e)]=g.set_index('year').val.sort_index()
def contrast(src,m,y0,y1,lag):
    D=np.log(src[(m,'Female','HCV')])-np.log(src[(m,'Female','HBV')])-np.log(src[(m,'Male','HCV')])+np.log(src[(m,'Male','HBV')])
    D=D.loc[y0:y1]; t=D.index.values-1990
    f=sm.OLS(D.values,sm.add_constant(t)).fit(cov_type='HAC',cov_kwds={'maxlags':lag}) if lag>0 else sm.OLS(D.values,sm.add_constant(t)).fit()
    b,se=f.params[1],f.bse[1]; return 100*b,100*(b-1.96*se),100*(b+1.96*se)
def within(src,m,s,y0,y1,lag=3):
    D=(np.log(src[(m,s,'HCV')])-np.log(src[(m,s,'HBV')])).loc[y0:y1]; t=D.index.values-1990
    f=sm.OLS(D.values,sm.add_constant(t)).fit(cov_type='HAC',cov_kwds={'maxlags':lag}); b,se=f.params[1],f.bse[1]; return 100*b,100*(b-1.96*se),100*(b+1.96*se)
rows=[]
for m in ['Incidence','Deaths']:
    for lag in [0,1,2,3,5]: rows.append(('GBD 2023, all-age ASR',m,'1990-2023',f'lag {lag}' if lag else 'OLS',*contrast(off,m,1990,2023,lag)))
    rows.append(('GBD 2023, all-age ASR',m,'1990-2019','lag 3',*contrast(off,m,1990,2019,3)))
    rows.append(('GBD 2023, all-age ASR',m,'1990-2021','lag 3',*contrast(off,m,1990,2021,3)))
# cross-round, ages 15+
buf=io.StringIO()
with contextlib.redirect_stdout(buf): ns=runpy.run_path(B+'gbd/consistency_2021_2023.py')
r21,r23=ns['r21'],ns['r23']
def todict(r):
    out={}
    for (m,s,e),g in r.groupby(level=[0,1,2]): out[(m,s,e)]=g.droplevel([0,1,2]).sort_index()
    return out
s21,s23=todict(r21),todict(r23)
for lab,src in [('GBD 2021, ages ≥15 ASR',s21),('GBD 2023, ages ≥15 ASR',s23)]:
    for m in ['Incidence','Deaths']:
        rows.append((lab,m,'1990-2021','lag 3',*contrast(src,m,1990,2021,3)))
R=pd.DataFrame(rows,columns=['source','measure','period','SE','diff_x100','lo','hi']).round(2)
W=[]
for lab,src,y1 in [('GBD 2023, all-age ASR',off,2023),('GBD 2021, ages ≥15 ASR',s21,2021),('GBD 2023, ages ≥15 ASR',s23,2021)]:
    for m in ['Incidence','Deaths']:
        for s in ['Female','Male']: W.append((lab,m,s,f'1990-{y1}',*within(src,m,s,1990,y1)))
W=pd.DataFrame(W,columns=['source','measure','sex','period','HCV_minus_HBV_slope_x100','lo','hi']).round(2)
# 1990 and last-year female ratio by round
rat=[]
for lab,src in [('GBD 2021 ≥15',s21),('GBD 2023 ≥15',s23)]:
    for m in ['Incidence','Deaths']:
        for s in ['Female','Male']:
            r=src[(m,s,'HCV')]/src[(m,s,'HBV')]; rat.append((lab,m,s,round(r.loc[1990],3),round(r.loc[2021],3)))
rat=pd.DataFrame(rat,columns=['source','measure','sex','HCV/HBV 1990','HCV/HBV 2021'])
R.to_csv(B+'out/contrast/ratio_contrast_robustness.csv',index=False); W.to_csv(B+'out/contrast/ratio_within_sex_by_round.csv',index=False); rat.to_csv(B+'out/contrast/ratio_levels_by_round.csv',index=False)
pd.set_option('display.width',250); print(R.to_string()); print(W.to_string()); print(rat.to_string())
