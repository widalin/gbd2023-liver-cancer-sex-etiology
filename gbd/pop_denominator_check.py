"""2026-11-19 review 9: check of GBD 2023 age-sex population derived as all-cause deaths / all-cause death rate.
(1) internal: liver cancer count / (rate/1e5) from coauthor_A (Number & Rate, same cell) vs derived pop.
(2) external: vs GBD 2021 official population (1990-2021) and UN WPP 2024 (1990-2023), ages 15+ by 5-y group and sex."""
import pandas as pd,numpy as np,glob
P=''
pop=pd.read_csv(P+'gbd/raw/gbd2023_china_pop/pop_derived_allcause.csv').rename(columns={'val':'pop23'})
AG=[f'{a}-{a+4} years' for a in range(15,95,5)]+['95+ years']
pop=pop[pop.age_name.isin(AG)]
a=pd.read_csv(glob.glob(P+'gbd/raw/coauthor_A*/*.csv')[0])
a=a[a.age_name.isin(AG)]
n=a[a.metric_name=='Number'].set_index(['measure_name','sex_name','age_name','year','cause_name']).val
r=a[a.metric_name=='Rate'].set_index(['measure_name','sex_name','age_name','year','cause_name']).val
imp=(n/(r/1e5)).dropna(); imp=imp[n.reindex(imp.index)>=50].rename('imp').reset_index()
m=imp.merge(pop,on=['sex_name','age_name','year']); m['rel']=m.imp/m.pop23-1
print('internal (count/rate, cells with >=50 cases) n=%d; |rel diff| median %.4f%%, 95th pct %.4f%%, max %.4f%%'%(len(m),100*m.rel.abs().median(),100*m.rel.abs().quantile(.95),100*m.rel.abs().max()))
g=pd.read_csv(P+'gbd/raw/gbd2021_china_pop/pop.csv'); g=g[g.age_name.isin(AG)][['sex_name','age_name','year','val']].rename(columns={'val':'pop21'})
c=pop.merge(g,on=['sex_name','age_name','year']); c['rel']=c.pop23/c.pop21-1
w=pd.read_csv(P+'gbd/raw/wpp2024/wpp2024_china_pop5_medium_1950_2100.csv')
w=w[(w.Time>=1990)&(w.Time<=2023)&(w.AgeGrpStart>=15)].copy()
w['age_name']=np.where(w.AgeGrpStart>=95,'95+ years',w.AgeGrpStart.astype(str)+'-'+(w.AgeGrpStart+4).astype(str)+' years')
w=w.groupby(['age_name','Time'],as_index=False)[['PopMale','PopFemale']].sum()
w=w.melt(['age_name','Time'],var_name='sex_name',value_name='wpp').replace({'PopMale':'Male','PopFemale':'Female'}).rename(columns={'Time':'year'}); w['wpp']*=1000
w.rename(columns={'wpp':'val'})[['sex_name','age_name','year','val']].to_csv(P+'out/audit/pop_wpp2024_15plus.csv',index=False)  # denominator for decomp.py WPP sensitivity
c2=pop.merge(w,on=['sex_name','age_name','year']); c2['rel']=c2.pop23/c2.wpp-1
def summ(df,lab):
    tot=df.groupby(['sex_name','year'])[[df.columns[3],df.columns[4]]].sum(); t=(tot.iloc[:,0]/tot.iloc[:,1]-1)
    old=df[df.age_name.isin(AG[-4:])]
    print(f'{lab}: 15+ totals rel diff range {100*t.min():.2f}% to {100*t.max():.2f}%; age-group |rel| median {100*df.rel.abs().median():.2f}%, 15-79 max {100*df[~df.age_name.isin(AG[-4:])].rel.abs().max():.2f}%, 80+ |rel| median {100*old.rel.abs().median():.2f}% max {100*old.rel.abs().max():.2f}%')
summ(c[['sex_name','age_name','year','pop23','pop21','rel']],'vs GBD 2021 (1990-2021)')
summ(c2[['sex_name','age_name','year','pop23','wpp','rel']],'vs UN WPP 2024 (1990-2023)')
# share of 15+ change: ageing proxy - age structure 1990 vs 2023 under each source
for lab,df,col in [('GBD2023 derived',pop,'pop23'),('WPP 2024',w,'wpp')]:
    s=df[df.year.isin([1990,2023])].groupby(['sex_name','year'])[col].sum().unstack(); print(lab,'15+ growth 1990->2023 %',(100*(s[2023]/s[1990]-1)).round(1).to_dict())
pd.concat([c.assign(ref='GBD2021'),c2.rename(columns={'wpp':'pop21'}).assign(ref='WPP2024')]).to_csv(P+'out/audit/pop_denominator_check.csv',index=False)
