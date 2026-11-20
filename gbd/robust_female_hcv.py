"""Robustness of the female HCV>HBV crossover, China, GBD 2023 (2026-11-10).
R1  ASR & counts with 95% UI, 4 measures: crossover year, years with non-overlapping UI.
R2  HCV/HBV ratio with approximate 95% CI (log-normal from UI, independence assumed -> conservative-ish caveat).
R3  Age-specific HCV/HBV ratio (5-yr ages 15-95+), 1990/2000/2010/2023, and crossover year by age.
R4  Alternative age-standardisation (15+): GBD standard, China female 2023 population, truncated 30-74 (Segi-like truncated), and age 15-64 vs 65+.
R5  Men as contrast.
Outputs: out/robust/female_hcv_robustness.xlsx, out/robust/FigS7_female_HCV_HBV_robustness.png"""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, os
B=''; OUT=B+'out/robust/'; os.makedirs(OUT,exist_ok=True)
m=pd.read_csv(B+'gbd/GBD2023_liver_cancer_etiology_Global_China_1990-2023.csv')
m=m[(m.location_name=='China')&m.cause_name.isin(['Liver cancer due to hepatitis B','Liver cancer due to hepatitis C'])]
m['et']=m.cause_name.map({'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV'})
m.loc[m.measure_name.str.startswith('DALY'),'measure_name']='DALYs'
z=1.959964
rows=[];ratio=[]
for sex in ['Female','Male']:
  for meas in ['Incidence','Deaths','DALYs','Prevalence']:
    for metric,age in [('Rate','Age-standardized'),('Number','All ages')]:
      g=m[(m.sex_name==sex)&(m.measure_name==meas)&(m.metric_name==metric)&(m.age_name==age)].pivot(index='year',columns='et',values=['val','lower','upper'])
      hv,cv=g['val']['HBV'],g['val']['HCV']
      above=cv>hv
      cross=next((y for y in g.index if above.loc[y:].all()),None)
      nonov=list(g.index[g['lower']['HCV']>g['upper']['HBV']])
      se=lambda e:(np.log(g['upper'][e])-np.log(g['lower'][e]))/(2*z)
      lr=np.log(cv/hv); s=np.sqrt(se('HCV')**2+se('HBV')**2)
      lo,hi=np.exp(lr-z*s),np.exp(lr+z*s)
      for y in g.index: ratio.append(dict(sex=sex,measure=meas,metric='ASR' if metric=='Rate' else 'Count',year=y,HBV=hv[y],HCV=cv[y],ratio=cv[y]/hv[y],ratio_lo=lo[y],ratio_hi=hi[y]))
      rows.append(dict(sex=sex,measure=meas,metric='ASR' if metric=='Rate' else 'Count',
        HBV_2023=f"{hv[2023]:.2f} ({g['lower']['HBV'][2023]:.2f}-{g['upper']['HBV'][2023]:.2f})",
        HCV_2023=f"{cv[2023]:.2f} ({g['lower']['HCV'][2023]:.2f}-{g['upper']['HCV'][2023]:.2f})",
        ratio_2023=f"{cv[2023]/hv[2023]:.2f} ({lo[2023]:.2f}-{hi[2023]:.2f})",
        sustained_crossover_year=cross if cross else 'none',
        years_UI_nonoverlap=(f"{nonov[0]}-{nonov[-1]} (n={len(nonov)})" if nonov else 'none'),
        years_ratio_CI_excl_1=int(((lo>1)).sum())))
R1=pd.DataFrame(rows); R2=pd.DataFrame(ratio)
# ---------- age-specific
a=pd.read_csv(B+'gbd/raw/coauthor_A_china_age_etiology_2023/IHME-GBD_2023_DATA-fc723b63-1.csv')
a=a[a.cause_name.isin(['Liver cancer due to hepatitis B','Liver cancer due to hepatitis C'])]
a['et']=a.cause_name.map({'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV'})
AG=[f'{x}-{x+4} years' for x in range(15,95,5)]+['95+ years']
r3=[]
for sex in ['Female','Male']:
  for meas in ['Incidence','Deaths']:
    g=a[(a.sex_name==sex)&(a.measure_name==meas)&(a.metric_name=='Rate')].pivot_table(index=['age_name','year'],columns='et',values='val')
    for ag in AG:
      s=g.loc[ag]; r=s.HCV/s.HBV; ab=r>1
      cross=next((y for y in s.index if ab.loc[y:].all()),None)
      r3.append(dict(sex=sex,measure=meas,age=ag,**{f'ratio_{y}':round(r[y],2) for y in [1990,2000,2010,2023]},sustained_crossover_year=cross if cross else ('always' if ab.all() else 'none')))
R3=pd.DataFrame(r3)
# ---------- alternative standardisation (15+)
std=pd.read_csv(B+'gbd/raw/gbd_std_pop_2019.csv').set_index('age_name').percent
pop=pd.read_csv(B+'gbd/raw/gbd2023_china_pop/pop_derived_allcause.csv')
def asr(sex,meas,ages,w):
  g=a[(a.sex_name==sex)&(a.measure_name==meas)&(a.metric_name=='Rate')&a.age_name.isin(ages)]
  ww=pd.Series(w).loc[ages]; ww=ww/ww.sum()
  g=g.assign(w=g.age_name.map(ww)); g['wv']=g.val*g.w
  return g.groupby(['year','et']).wv.sum().unstack()
r4=[];curves={}
fp=pop[(pop.sex_name=='Female')&(pop.year==2023)].set_index('age_name').val
SCH={'GBD std, 15+':(AG,std),'China female 2023 pop, 15+':(AG,fp),'GBD std, 30-74':([f'{x}-{x+4} years' for x in range(30,75,5)],std),
     'GBD std, 15-64':([f'{x}-{x+4} years' for x in range(15,65,5)],std),'GBD std, 65+':([f'{x}-{x+4} years' for x in range(65,95,5)]+['95+ years'],std)}
for meas in ['Incidence','Deaths']:
  for k,(ages,w) in SCH.items():
    s=asr('Female',meas,ages,w); ab=s.HCV>s.HBV
    cross=next((y for y in s.index if ab.loc[y:].all()),None)
    curves[(meas,k)]=s
    r4.append(dict(measure=meas,standard=k,HBV_2023=round(s.HBV[2023],2),HCV_2023=round(s.HCV[2023],2),ratio_2023=round(s.HCV[2023]/s.HBV[2023],2),sustained_crossover_year=cross if cross else 'none'))
R4=pd.DataFrame(r4)
with pd.ExcelWriter(OUT+'female_hcv_robustness.xlsx') as w:
  pd.DataFrame({'Note':['GBD 2023, China. Crossover year = first year from which HCV>HBV every year through 2023.','Ratio CI: log-normal from 95% UI, HCV and HBV assumed independent (GBD draws not available); approximation only.','Age-specific rates from co-author download (5-yr ages). ASR in R4 recomputed from age-specific rates, 15+ only (no <15 data).','Population for China-female standard: derived from all-cause deaths/rate (pop_derived_allcause.csv).']}).to_excel(w,sheet_name='Notes',index=False)
  R1.to_excel(w,sheet_name='R1_crossover_UI',index=False); R2.round(4).to_excel(w,sheet_name='R2_ratio_by_year',index=False)
  R3.to_excel(w,sheet_name='R3_age_specific',index=False); R4.to_excel(w,sheet_name='R4_alt_standard',index=False)
# ---------- figure
C={'HBV':'#4BA3D9','HCV':'#7FA67A'}
fig,ax=plt.subplots(1,3,figsize=(15,4.4))
for meas,ls in [('Incidence','-'),('Deaths','--')]:
  d=R2[(R2.sex=='Female')&(R2.measure==meas)&(R2.metric=='ASR')]
  ax[0].plot(d.year,d.ratio,ls,color='#333',label=meas); ax[0].fill_between(d.year,d.ratio_lo,d.ratio_hi,color='#999',alpha=.25)
ax[0].axhline(1,color='r',lw=.8); ax[0].set_title('A  Women: HCV/HBV ratio of ASR (approx. 95% CI)'); ax[0].legend(frameon=False); ax[0].set_xlabel('Year')
f=R3[(R3.sex=='Female')&(R3.measure=='Incidence')]
x=np.arange(len(f)); 
for y,c in [(1990,'#bbb'),(2010,'#888'),(2023,'#222')]: ax[1].plot(x,f[f'ratio_{y}'],'o-',ms=3,color=c,label=str(y))
ax[1].axhline(1,color='r',lw=.8); ax[1].set_yscale('log'); ax[1].set_yticks([0.02,0.05,0.1,0.2,0.5,1,2,5]); ax[1].set_yticklabels(['0.02','0.05','0.1','0.2','0.5','1','2','5']); ax[1].minorticks_off(); ax[1].set_xticks(x[::2]); ax[1].set_xticklabels([s.split(' ')[0] for s in f.age][::2],rotation=45)
ax[1].set_title('B  Women: age-specific HCV/HBV incidence ratio'); ax[1].legend(frameon=False); ax[1].set_xlabel('Age group')
for k,c in [('GBD std, 15+','#222'),('China female 2023 pop, 15+','#C2703D'),('GBD std, 15-64','#4BA3D9'),('GBD std, 65+','#7FA67A')]:
  s=curves[('Incidence',k)]; ax[2].plot(s.index,s.HCV/s.HBV,color=c,label=k)
ax[2].axhline(1,color='r',lw=.8); ax[2].set_ylim(0,2.4); ax[2].set_title('C  Women: HCV/HBV incidence ratio by standard'); ax[2].legend(frameon=False,fontsize=8,loc='upper left',ncol=2); ax[2].set_xlabel('Year')
plt.tight_layout(); plt.savefig(OUT+'FigS7_female_HCV_HBV_robustness.png',dpi=300)
pd.set_option('display.width',250); pd.set_option('display.max_columns',20)
print(R1.to_string()); print(R3[R3.sex=='Female'].to_string()); print(R4.to_string())
print(R3[(R3.sex=='Male')].sustained_crossover_year.value_counts())
