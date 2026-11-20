"""A5: test whether the gap between APC net drift and ASR trend is due to age weighting.
For each sex x etiology x measure (GBD 2023, China, ages 20-94, 1994-2023, same as APC main model):
 (a) EAPC of ASR with GBD 2019 standard weights; (b) EAPC of equal-weight arithmetic mean rate;
 (c) EAPC of equal-weight geometric mean rate (= mean of age-specific log trends, the APC net drift analogue); vs APC net drift.
Output: out/robust/equal_weight_asr.xlsx"""
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
B=''
ET={'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV','Liver cancer due to alcohol use':'Alcohol','Liver cancer due to NASH':'MASH','Liver cancer due to other causes':'Other'}
a=pd.read_csv(B+'gbd/raw/coauthor_A_china_age_etiology_2023/IHME-GBD_2023_DATA-fc723b63-1.csv')
a=a[a.cause_name.isin(ET)&(a.metric_name=='Rate')].copy(); a['et']=a.cause_name.map(ET)
AG=[f'{x}-{x+4} years' for x in range(20,95,5)]
a=a[a.age_name.isin(AG)&a.year.between(1994,2023)]
std=pd.read_csv(B+'gbd/raw/gbd_std_pop_2019.csv').set_index('age_name').percent; w=std.loc[AG]/std.loc[AG].sum()
ap=pd.read_csv(B+'out/apc/apc_summary_gbd2023.csv')
def ea(y,yr):
    f=sm.OLS(np.log(y),sm.add_constant(yr)).fit(); b,se=f.params[1],f.bse[1]; t=stats.t.ppf(.975,len(y)-2)
    return [100*(np.exp(v)-1) for v in (b,b-t*se,b+t*se)]
rows=[]
for (m,s,e),g in a.groupby(['measure_name','sex_name','et']):
    P=g.pivot(index='year',columns='age_name',values='val')[AG]; yr=P.index.values
    std_asr=(P*w.values).sum(1); eq=P.mean(1); geo=np.exp(np.log(P).mean(1))
    # share of std-weighted ASR contributed by ages >=60 in 2023
    c=(P.loc[2023]*w.values); sh60=100*c[[x for x in AG if int(x[:2])>=60]].sum()/c.sum()
    nd=ap[(ap.measure==m)&(ap.sex==s)&(ap.etiology==e)].iloc[0]
    r=dict(measure=m,sex=s,etiology=e)
    for k,v in (('std',std_asr),('equal_arith',eq),('equal_geo',geo)):
        r[f'EAPC_{k}'],r[f'{k}_lo'],r[f'{k}_hi']=ea(v.values,yr)
    r.update(APC_net_drift=nd.net_drift,nd_lo=nd.nd_lo,nd_hi=nd.nd_hi,share_ASR_from_60plus_2023=sh60)
    rows.append(r)
R=pd.DataFrame(rows)
R['gap_std_minus_netdrift']=R.EAPC_std-R.APC_net_drift; R['gap_geo_minus_netdrift']=R.EAPC_equal_geo-R.APC_net_drift
R.to_excel(B+'out/robust/equal_weight_asr.xlsx',index=False)
pd.set_option('display.width',250); print(R[['measure','sex','etiology','EAPC_std','EAPC_equal_arith','EAPC_equal_geo','APC_net_drift','gap_std_minus_netdrift','gap_geo_minus_netdrift','share_ASR_from_60plus_2023']].round(2).to_string())
print('corr(gap_geo vs nd) max abs', R.gap_geo_minus_netdrift.abs().max().round(3), ' mean abs gap std', R.gap_std_minus_netdrift.abs().mean().round(3))
