"""GBD 2021 vs GBD 2023 consistency: ages 15+ ASR (GBD 2019 std pop, renormalised to 15+) from age-specific counts/pop, EAPC 1990-2021."""
import pandas as pd,numpy as np,glob
from scipy import stats
B=''
ET={'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV','Liver cancer due to alcohol use':'Alcohol','Liver cancer due to NASH':'MASH','Liver cancer due to other causes':'Other'}
AG=[f'{a}-{a+4} years' for a in range(15,95,5)]+['95+ years']
sp=pd.read_csv(B+'gbd/raw/gbd_std_pop_2019.csv')
w=sp.set_index('age_name').percent/100
def asr(cnt,pop):
    c=cnt[cnt.cause_name.isin(ET)&cnt.age_name.isin(AG)&(cnt.metric_name=='Number')].copy(); c['et']=c.cause_name.map(ET)
    a=c.groupby(['measure_name','sex_name','age_name','year'],as_index=False).val.sum(); a['et']='All'
    c=pd.concat([c[['measure_name','sex_name','age_name','year','et','val']],a])
    c=c.merge(pop[['sex_name','age_name','year','val']].rename(columns={'val':'pop'}),on=['sex_name','age_name','year'])
    c['x']=c.val/c['pop']*1e5*c.age_name.map(w)
    return c.groupby(['measure_name','sex_name','et','year']).x.sum()
r21=asr(pd.read_csv(B+'gbd/raw/gbd2021_china_age/IHME-GBD_2021_DATA-2cc0e067-1.csv'),pd.read_csv(B+'gbd/raw/gbd2021_china_pop/pop.csv'))
r23=asr(pd.read_csv(glob.glob(B+'gbd/raw/coauthor_A*/*.csv')[0]),pd.read_csv(B+'gbd/raw/gbd2023_china_pop/pop_derived_allcause.csv'))
def eapc(s):
    s=s[(s.index>=1990)&(s.index<=2021)]; r=stats.linregress(s.index,np.log(s.values)); t=stats.t.ppf(.975,len(s)-2); f=lambda v:100*(np.exp(v)-1)
    return f(r.slope),f(r.slope-t*r.stderr),f(r.slope+t*r.stderr)
rows=[]
for k in r21.index.droplevel(3).unique():
    a=r21.loc[k]; b=r23.loc[k]; e1=eapc(a); e2=eapc(b)
    rows.append(dict(measure=k[0],sex=k[1],etiology=k[2],ASR15_1990_g21=a.loc[1990],ASR15_1990_g23=b.loc[1990],ASR15_2021_g21=a.loc[2021],ASR15_2021_g23=b.loc[2021],
      EAPC_g21=e1[0],EAPC_g21_lo=e1[1],EAPC_g21_hi=e1[2],EAPC_g23=e2[0],EAPC_g23_lo=e2[1],EAPC_g23_hi=e2[2],
      same_direction=np.sign(e1[0])==np.sign(e2[0]),corr_log=np.corrcoef(np.log(a.loc[1990:2021]),np.log(b.loc[1990:2021]))[0,1]))
R=pd.DataFrame(rows); R.to_csv(B+'out/apc/consistency_gbd2021_vs_gbd2023.csv',index=False)
pd.set_option('display.width',250); print(R.round(2).to_string())
# validation vs official GBD 2023 ASR (all ages)
o=pd.read_csv(B+'gbd/IHME-GBD_2023_DATA-4f95895a-1.csv'); o=o[(o.location_name=='China')&(o.age_name=='Age-standardized')&(o.metric_name=='Rate')&o.cause_name.isin(ET)]
o['et']=o.cause_name.map(ET); o=o.set_index(['measure_name','sex_name','et','year']).val
v=(r23.reindex(o.index)/o-1).dropna(); print('15+ derived vs official GBD2023 ASR rel diff:',v.describe().round(4).to_dict())
