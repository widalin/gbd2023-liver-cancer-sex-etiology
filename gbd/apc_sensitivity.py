"""APC sensitivity (2026-11-10): main (20-94, 1994-2023) vs ages 30-84, 20-84, periods 1994-2018 (pre-COVID), 1991-2020.
Also reconciles APC net drift with EAPC of truncated (20-94) ASR over the same years. Outputs out/robust/apc_sensitivity.xlsx + FigS8."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B=''; S=B+'out/robust/apc_sens/'
V={'Main 20-94, 1994-2023':B+'out/apc/apc_%s_gbd2023.csv','Ages 30-84':S+'apc_%s_sens_age30_84.csv','Ages 20-84':S+'apc_%s_sens_age20_84.csv',
   'Periods 1994-2018':S+'apc_%s_sens_1994_2018.csv','Periods 1991-2020':S+'apc_%s_sens_1991_2020.csv'}
rows=[]
for k,f in V.items():
    d=pd.read_csv(f%'summary'); d['variant']=k; rows.append(d)
R=pd.concat(rows); R=R[R.etiology!='Other']
R['nd']=R.apply(lambda r:f"{r.net_drift:.2f} ({r.nd_lo:.2f}, {r.nd_hi:.2f})"+('' if (r.nd_lo>0 or r.nd_hi<0) else ' ns'),axis=1)
T=R.pivot_table(index=['measure','sex','etiology'],columns='variant',values='nd',aggfunc='first')[list(V)]
# cohort RR at selected cohorts
F=[];
for k,f in V.items():
    d=pd.read_csv(f%'functions'); d=d[d.function=='cohort_RR']; d['variant']=k; F.append(d)
F=pd.concat(F)
# truncated ASR EAPC reconciliation (20-94, GBD std, 1994-2023)
a=pd.read_csv(B+'gbd/raw/coauthor_A_china_age_etiology_2023/IHME-GBD_2023_DATA-fc723b63-1.csv')
ET={'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV','Liver cancer due to alcohol use':'Alcohol','Liver cancer due to NASH':'MASH'}
std=pd.read_csv(B+'gbd/raw/gbd_std_pop_2019.csv').set_index('age_name').percent
AG=[f'{x}-{x+4} years' for x in range(20,95,5)]; w=std.loc[AG]/std.loc[AG].sum()
a=a[a.cause_name.isin(ET)&a.age_name.isin(AG)&(a.metric_name=='Rate')]; a['et']=a.cause_name.map(ET); a['wv']=a.val*a.age_name.map(w)
s=a.groupby(['measure_name','sex_name','et','year']).wv.sum().reset_index()
rec=[]
for (m,sx,e),g in s.groupby(['measure_name','sex_name','et']):
    for y0,y1 in [(1994,2023),(1994,2018),(1990,2023)]:
        h=g[(g.year>=y0)&(g.year<=y1)]; b,c=np.polyfit(h.year,np.log(h.wv),1)
        X=np.c_[np.ones(len(h)),h.year]; r=np.log(h.wv)-X@[c,b]; se=np.sqrt(r@r/(len(h)-2)*np.linalg.inv(X.T@X)[1,1])
        rec.append(dict(measure=m,sex=sx,etiology=e,window=f'{y0}-{y1}',EAPC_trunc20_94=f"{100*(np.exp(b)-1):.2f} ({100*(np.exp(b-1.96*se)-1):.2f}, {100*(np.exp(b+1.96*se)-1):.2f})"))
REC=pd.DataFrame(rec).pivot_table(index=['measure','sex','etiology'],columns='window',values='EAPC_trunc20_94',aggfunc='first')
nd_main=R[R.variant=='Main 20-94, 1994-2023'].set_index(['measure','sex','etiology']).nd
REC['APC_net_drift_main']=nd_main
with pd.ExcelWriter(B+'out/robust/apc_sensitivity.xlsx') as x:
    pd.DataFrame({'Note':['Net drift %/yr (95% CI); ns = CI includes 0. GBD 2023 China, in-house Python APC (Rosenberg 2014 parameterisation, Poisson quasi-likelihood).',
      'Periods 1994-2018 = drops 2019-2023 (COVID-era) period. Periods 1991-2020 = shifts 5-yr grouping by 3 years.',
      'Reconciliation: EAPC of ASR restricted to ages 20-94 (GBD standard, 2019 weights), log-linear OLS, to compare with APC net drift on same age range.']}).to_excel(x,sheet_name='Notes',index=False)
    T.reset_index().to_excel(x,sheet_name='Net_drift',index=False)
    F.pivot_table(index=['measure','sex','etiology','x'],columns='variant',values='est').reset_index().round(3).to_excel(x,sheet_name='Cohort_RR',index=False)
    REC.reset_index().to_excel(x,sheet_name='Reconcile_EAPC',index=False)
pd.set_option('display.width',300); pd.set_option('display.max_columns',20); pd.set_option('display.max_colwidth',24)
print(T.to_string()); print(REC.to_string())
# figure (Supplementary Figure S2): cohort RR, incidence; every curve re-referenced to the 1953.5 birth cohort
# (divided by its own RR at 1953.5; log-linear interpolation between 1950.5 and 1955.5 for periods 1991-2020, whose grid is offset by 3 years)
def _reref(d,ref=1953.5):
    d=d.sort_values('x'); x=d.x.values.astype(float); y=np.log(d.est.values.astype(float))
    r=y[x==ref][0] if (x==ref).any() else np.interp(ref,x,y)
    return x,np.exp(y-r)
fig,ax=plt.subplots(2,2,figsize=(10.5,7.6),sharey=False); ax=ax.ravel()
col=dict(zip(V,['#222','#C2703D','#E0A458','#4BA3D9','#7FA67A'])); RR=[]
for i,(sx,e) in enumerate([('Male','MASH'),('Male','Alcohol'),('Female','HCV'),('Female','HBV')]):
    for k in V:
        d=F[(F.variant==k)&(F.measure=='Incidence')&(F.sex==sx)&(F.etiology==e)]
        x_,y_=_reref(d); RR+=[dict(sex=sx,etiology=e,variant=k,cohort=a,RR_ref1953_5=round(b,4)) for a,b in zip(x_,y_)]
        ax[i].plot(x_,y_,color=col[k],lw=2 if k.startswith('Main') else 1.2,label=k.replace('-','\u2013'))
    ax[i].axhline(1,color='#1d1d1f55',lw=.8,ls=':'); ax[i].set_yscale('log'); ax[i].set_xlabel('Birth cohort')
    ax[i].set_title(f'{"ABCD"[i]}  {"Men" if sx=="Male" else "Women"}, {e}: cohort RR (incidence)',loc='left')
    ax[i].set_yticks([0.1,0.2,0.5,1,2,5]); ax[i].set_yticklabels(['0.1','0.2','0.5','1','2','5']); ax[i].minorticks_off()
ax[0].legend(frameon=False,fontsize=7.5,loc='lower left')
plt.tight_layout(); plt.savefig(B+'out/robust/FigS2_APC_sensitivity_cohort.png',dpi=600)
pd.DataFrame(RR).to_csv(B+'out/robust/apc_sensitivity_cohortRR_ref1953.5.csv',index=False)
