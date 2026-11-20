"""A3: age-specific trends of MASH-related liver cancer in China (GBD 2023), focus on older men.
Uses GBD 2023 age-specific Rate (per 100,000) with 95% UI from the co-author download.
Outputs: out/robust/age_specific_mash.xlsx, out/robust/FigS9_age_specific_MASH.png"""
import numpy as np, pandas as pd, statsmodels.api as sm, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
B=''
ET={'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV','Liver cancer due to alcohol use':'Alcohol','Liver cancer due to NASH':'MASH','Liver cancer due to other causes':'Other'}
a=pd.read_csv(B+'gbd/raw/coauthor_A_china_age_etiology_2023/IHME-GBD_2023_DATA-fc723b63-1.csv')
a=a[a.cause_name.isin(ET)].copy(); a['et']=a.cause_name.map(ET)
a['age_lo']=a.age_name.str.extract(r'^(\d+)').astype(int)
R=a[a.metric_name=='Rate']; N=a[a.metric_name=='Number']

def eapc(g,y0,y1):
    g=g[(g.year>=y0)&(g.year<=y1)].sort_values('year'); y=np.log(g.val.values); X=sm.add_constant(g.year.values)
    f=sm.OLS(y,X).fit(); b=f.params[1]; se=f.bse[1]; t=stats.t.ppf(.975,len(g)-2)
    return 100*(np.exp(b)-1),100*(np.exp(b-t*se)-1),100*(np.exp(b+t*se)-1)

rows=[]
for (m,s,e,ag,lo),g in R.groupby(['measure_name','sex_name','et','age_name','age_lo']):
    r={'measure':m,'sex':s,'etiology':e,'age':ag,'age_lo':lo}
    for tag,(y0,y1) in {'1990_2023':(1990,2023),'1994_2023':(1994,2023),'1990_2019':(1990,2019)}.items():
        r[f'EAPC_{tag}'],r[f'lo_{tag}'],r[f'hi_{tag}']=eapc(g,y0,y1)
    g0=g[g.year==1990].iloc[0]; g1=g[g.year==2023].iloc[0]
    r.update(rate1990=g0.val,rate1990_lo=g0.lower,rate1990_hi=g0.upper,rate2023=g1.val,rate2023_lo=g1.lower,rate2023_hi=g1.upper,
             pct_change=100*(g1.val/g0.val-1),UI_overlap=bool(not(g1.lower>g0.upper or g1.upper<g0.lower)))
    rows.append(r)
E=pd.DataFrame(rows).sort_values(['measure','sex','etiology','age_lo'])

# APC local drift (main model) for comparison
f=pd.read_csv(B+'out/apc/apc_functions_nci2023.csv'); ld=f[f.function=='local_drift'].copy()
ld['age_lo']=(ld.x-2.5).round().astype(int)
E=E.merge(ld[['measure','sex','etiology','age_lo','est','lo','hi']].rename(columns={'est':'APC_local_drift','lo':'LD_lo','hi':'LD_hi'}),on=['measure','sex','etiology','age_lo'],how='left')

# share of male MASH cases by age band
n=N[(N.et=='MASH')].copy(); n['band']=pd.cut(n.age_lo,[0,45,60,70,80,200],right=False,labels=['15-44','45-59','60-69','70-79','80+'])
S=n.groupby(['measure_name','sex_name','year','band'],observed=True).val.sum().reset_index()
S['share']=S.val/S.groupby(['measure_name','sex_name','year']).val.transform('sum')*100

# pop consistency check: Number/Rate implied pop vs derived pop
chk=N.merge(R,on=['measure_name','sex_name','et','age_name','year'],suffixes=('_n','_r'))
chk['pop_implied']=chk.val_n/chk.val_r*1e5
pc=chk.groupby(['sex_name','age_name','year']).pop_implied.agg(['min','max']); pc['rel_spread_%']=100*(pc['max']/pc['min']-1)

with pd.ExcelWriter(B+'out/robust/age_specific_mash.xlsx') as w:
    E[E.etiology=='MASH'].round(3).to_excel(w,sheet_name='MASH_age_EAPC',index=False)
    E.round(3).to_excel(w,sheet_name='All_etiologies_age_EAPC',index=False)
    S.round(2).to_excel(w,sheet_name='MASH_case_share_by_age',index=False)
    pc.reset_index().round(4).to_excel(w,sheet_name='Check_pop_consistency',index=False)
    pd.DataFrame({'Notes':[
     'Source: GBD 2023 China, age-specific Rate (per 100,000) and Number with 95% UI, co-author download.',
     'EAPC = 100*(exp(b)-1), ln(rate) ~ year OLS, t-based 95% CI. Three windows: 1990-2023, 1994-2023 (APC window), 1990-2019 (excludes 2020-2023 rebound).',
     'APC_local_drift = main APC model (ages 20-94, 1994-2023), for comparison; 15-19 and 95+ not in APC.',
     'UI_overlap: whether 1990 and 2023 95% UIs overlap (descriptive; not a formal test).',
     'Check_pop_consistency: population implied by Number/Rate across 5 etiologies x 2 measures for the same sex/age/year; spread should be ~0.']}).to_excel(w,sheet_name='Notes',index=False)

# ---------- figure ----------
from matplotlib.ticker import FuncFormatter,NullFormatter,FixedLocator
ECOL={'HBV':'#7f93b8','HCV':'#607868','Alcohol':'#c9935e','MASH':'#6d5a96','Other':'#c8b04a'}
ages_show=['50-54 years','60-64 years','70-74 years','80-84 years','90-94 years','95+ years']
cmap=plt.get_cmap('viridis'); cc={a_:cmap(i/(len(ages_show)-1)*0.9) for i,a_ in enumerate(ages_show)}
fig,ax=plt.subplots(2,2,figsize=(13,9.5))
A=ax[0,0]
for ag in ages_show:
    g=R[(R.measure_name=='Incidence')&(R.sex_name=='Male')&(R.et=='MASH')&(R.age_name==ag)].sort_values('year')
    A.plot(g.year,g.val,color=cc[ag],lw=1.9,label=ag.replace(' years',''))
    A.fill_between(g.year,g.lower,g.upper,color=cc[ag],alpha=.10,lw=0)
A.set_yscale('log'); A.yaxis.set_major_locator(FixedLocator([0.5,1,2,3,5,8]))
A.yaxis.set_major_formatter(FuncFormatter(lambda v,_: f'{v:g}')); A.yaxis.set_minor_formatter(NullFormatter())
A.set_title('A  Male MASH-related incidence rate by age (per 100,000)',loc='left',fontsize=10.5)
A.set_xlabel('Year'); A.legend(fontsize=7.5,ncol=3,title='Age (years)',title_fontsize=8,frameon=False,loc='lower right')
Bx=ax[0,1]
sty={('Male','Incidence'):('#6d5a96','-','o'),('Male','Deaths'):('#6d5a96','--','s'),('Female','Incidence'):('#c9935e','-','o'),('Female','Deaths'):('#c9935e','--','s')}
for (s_,m),(c,ls,mk) in sty.items():
    g=E[(E.etiology=='MASH')&(E.sex==s_)&(E.measure==m)&(E.age_lo>=20)].sort_values('age_lo')
    x=g.age_lo+2.5+(0.6 if m=='Deaths' else -0.6)
    Bx.errorbar(x,g.EAPC_1990_2023,yerr=[g.EAPC_1990_2023-g.lo_1990_2023,g.hi_1990_2023-g.EAPC_1990_2023],color=c,ls=ls,marker=mk,ms=3.5,lw=1.3,capsize=1.5,label=f'{s_}, {m.lower()}')
Bx.axhline(0,color='grey',lw=.8,ls=':'); Bx.set_xlabel('Age group (years, midpoint; 95+ plotted at 97.5)'); Bx.set_ylabel('EAPC 1990-2023 (%/year, 95% CI)')
Bx.set_title('B  MASH-related rates: EAPC by age, men vs women',loc='left',fontsize=10.5); Bx.legend(fontsize=8,frameon=False)
C=ax[1,0]
for e,c in ECOL.items():
    g=E[(E.etiology==e)&(E.sex=='Male')&(E.measure=='Incidence')&(E.age_lo>=20)].sort_values('age_lo')
    C.plot(g.age_lo+2.5,g.EAPC_1990_2023,color=c,marker='o',ms=3,lw=2.2 if e=='MASH' else 1.4,label=e)
    if e=='MASH': C.fill_between(g.age_lo+2.5,g.lo_1990_2023,g.hi_1990_2023,color=c,alpha=.15,lw=0)
C.axhline(0,color='grey',lw=.8,ls=':'); C.axvspan(80,100,color='#eeeeee',zorder=0)
C.set_xlabel('Age group (years, midpoint)'); C.set_ylabel('EAPC 1990-2023 (%/year)')
C.set_title('C  Male incidence: EAPC by age, all etiologies',loc='left',fontsize=10.5); C.legend(fontsize=8,frameon=False,ncol=2,loc='lower right')
D=ax[1,1]; bc={'15-44':'#d9d9d9','45-59':'#a6a6a6','60-69':'#7f93b8','70-79':'#c9935e','80+':'#6d5a96'}
for b,c in bc.items():
    g=S[(S.measure_name=='Incidence')&(S.sex_name=='Male')&(S.band==b)].sort_values('year')
    D.plot(g.year,g.share,color=c,lw=2.4 if b=='80+' else 1.5,label=b)
D.set_ylabel('% of male MASH-related incident cases'); D.set_xlabel('Year')
D.set_title('D  Age distribution of male MASH-related incident cases',loc='left',fontsize=10.5); D.legend(fontsize=8,frameon=False,title='Age (years)',title_fontsize=8,ncol=5,loc='upper center',bbox_to_anchor=(.5,-.13))
for A_ in ax.flat: A_.spines[['top','right']].set_visible(False); A_.tick_params(labelsize=8.5)
fig.suptitle('Age-specific trends in MASH-related liver cancer in China, 1990-2023 (GBD 2023)',x=.01,ha='left',fontsize=12.5)
fig.text(.01,.005,'Shaded bands: 95% uncertainty intervals (A) and 95% CI of EAPC for MASH (C). EAPC from log-linear regression of age-specific rates on calendar year.',fontsize=8,color='#555')
plt.tight_layout(rect=(0,.02,1,.96)); plt.savefig(B+'out/robust/FigS9_age_specific_MASH.png',dpi=200)
pd.set_option('display.width',250)
print(E[(E.etiology=='MASH')&(E.sex=='Male')][['measure','age','rate1990','rate2023','pct_change','UI_overlap','EAPC_1990_2023','lo_1990_2023','hi_1990_2023','EAPC_1994_2023','EAPC_1990_2019','APC_local_drift']].round(2).to_string(index=False))
print(S[(S.sex_name=='Male')&S.year.isin([1990,2023])].round(1).to_string(index=False))
print('pop consistency max spread %',pc['rel_spread_%'].max().round(3))
