import sys,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B=''
TAG=sys.argv[1] if len(sys.argv)>1 else 'gbd2021'; MEAS=sys.argv[2] if len(sys.argv)>2 else 'Incidence'
COL={"HBV":"#7f93b8","HCV":"#607868","Alcohol":"#c9935e","MASH":"#8a7fa6","Other":"#c2a94a","All":"#333333"}
K=pd.read_csv(B+f'out/apc/apc_functions_{TAG}.csv'); S=pd.read_csv(B+f'out/apc/apc_summary_{TAG}.csv')
YR={'gbd2021':'1992-2021','gbd2023':'1994-2023','nci2023':'1994-2023'}[TAG]
REF={'gbd2021':('2007-2011','1951.5'),'gbd2023':('2009-2013','1953.5'),'nci2023':('2004-2008','1949')}[TAG]
K=K[K.measure==MEAS]; S=S[S.measure==MEAS]
fn=[('local_drift','Local drift (% per year)','Age (years)',False),('long_age_rate','Longitudinal age curve (per 100,000)','Age (years)',True),('period_RR','Period rate ratio','Period (midpoint)',False),('cohort_RR','Cohort rate ratio','Birth cohort (midpoint)',True)]
fig,ax=plt.subplots(2,4,figsize=(17,8.2))
for r,sex in enumerate(['Male','Female']):
    for cI,(f,t,xl,lg) in enumerate(fn):
        a=ax[r,cI]
        for e in ['HBV','HCV','Alcohol','MASH','Other','All']:
            g=K[(K.sex==sex)&(K.etiology==e)&(K.function==f)]
            nd=S[(S.sex==sex)&(S.etiology==e)].iloc[0]
            lab=f"{e} (net drift {nd.net_drift:+.2f}%)" if cI==0 else e
            a.plot(g.x,g.est,color=COL[e],lw=2.2 if e in('HCV','MASH') else 1.4,ls='--' if e=='All' else '-',label=lab)
            if e in ('HCV','MASH'): a.fill_between(g.x,g.lo,g.hi,color=COL[e],alpha=.15,lw=0)
        if lg:
            a.set_yscale('log'); from matplotlib.ticker import FuncFormatter, NullFormatter; a.yaxis.set_major_formatter(FuncFormatter(lambda v,_: f'{v:g}')); a.yaxis.set_minor_formatter(NullFormatter())
        a.axhline(0 if f=='local_drift' else 1,color='grey',lw=.8,ls=':') if f!='long_age_rate' else None
        a.set_title(f"{sex}: {t}",fontsize=10.5,loc='left'); a.set_xlabel(xl,fontsize=9); a.tick_params(labelsize=8.5)
        a.spines[['top','right']].set_visible(False)
        if cI==0: a.legend(fontsize=7.5,frameon=False)
fig.suptitle(f"Age-period-cohort effects on liver cancer {MEAS.lower()} by etiology and sex in China ({'GBD 2023' if TAG=='nci2023' else TAG.upper().replace('GBD','GBD ')}, ages 20-94, {YR})",fontsize=13,x=0.01,ha='left')
fig.text(0.01,0.005,(f"Shaded bands: 95% CI for HCV and MASH. Reference period {REF[0]}; reference birth cohort midpoint {REF[1]}. Six 5-year periods; "+('estimated with the NCI Age-Period-Cohort Web Tool (default settings).' if TAG=='nci2023' else 'Poisson APC with overdispersion.')),fontsize=8.5,color='#555')
plt.tight_layout(rect=(0,0.02,1,0.95)); fn_=B+f'out/apc/Fig5_APC_{MEAS}_{TAG}.png'; plt.savefig(fn_,dpi=200); print(fn_)
