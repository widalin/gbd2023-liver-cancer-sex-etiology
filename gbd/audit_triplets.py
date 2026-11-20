"""Triplet audit: every 'est (lo to hi)' / 'est [lo to hi]' / 'est (lo--hi)' in Abstract..Conclusions and figure/table legends
must match (to the printed precision, all three numbers) a triplet from the source pool. Unmatched -> out/audit/triplets_unmatched.csv"""
import re, sys, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
B=''; sys.path.insert(0,B+'gbd'); import nci
md=open(B+'manuscript/manuscript_v2.9.md').read().replace('\u2212','-')
pool=[]
def add(t,src):
    if all(pd.notna(v) for v in t): pool.append((tuple(float(v) for v in t),src))
# NCI AAPC + APC segments (all batches)
for _,r in nci.AAPC.iterrows(): add((r['AAPC'],r['AAPC C.I. Low'],r['AAPC C.I. High']),'NCI AAPC '+str(r['name']))
for _,r in nci.APC.iterrows(): add((r['APC'],r['APC 95% LCL'],r['APC 95% UCL']),'NCI APC '+str(r['name']))
S3=pd.read_csv(B+'out/nci_joinpoint/rerun3_output/China_24_sensitivity_NCI_wide.csv')
for _,r in S3.iterrows():
    for c in [c for c in S3.columns if c.startswith('AAPC_')]:
        s=c[5:]; 
        if 'lo_'+s in S3: add((r[c],r['lo_'+s],r['hi_'+s]),'S1 '+r['name']+' '+s)
# EAPC from raw, all loc/meas/sex/cause, 1990-2023 and 1990-2019
M=pd.read_csv(B+'gbd/GBD2023_liver_cancer_etiology_Global_China_1990-2023.csv')
A=M[(M.metric_name=='Rate')&(M.age_name=='Age-standardized')]
for k,x in A.groupby(['location_name','measure_name','sex_name','cause_name']):
    x=x.sort_values('year'); add((x.val.iloc[0],x.lower.iloc[0],x.upper.iloc[0]),f'ASR1990 {k}'); add((x.val.iloc[-1],x.lower.iloc[-1],x.upper.iloc[-1]),f'ASR2023 {k}')
    for y0,y1 in [(1990,2023),(1990,2019)]:
        z=x[x.year.between(y0,y1)]; f=sm.OLS(np.log(z.val.values),sm.add_constant(z.year.values)).fit(); b,se=f.params[1],f.bse[1]; t=stats.t.ppf(.975,len(z)-2)
        add([100*(np.exp(v)-1) for v in (b,b-t*se,b+t*se)],f'EAPC {y0}-{y1} {k}')
N=M[(M.metric_name=='Number')&(M.age_name=='All ages')]
for _,r in N[N.year.isin([1990,2023])].iterrows(): add((r.val,r.lower,r.upper),f'N {r.year} {r.location_name} {r.measure_name} {r.sex_name} {r.cause_name}')
# APC
ap=pd.read_csv(B+'out/apc/apc_summary_nci2023.csv')
for _,r in ap.iterrows(): add((r.net_drift,r.nd_lo,r.nd_hi),f'net drift {r.measure} {r.sex} {r.etiology}')
fn=pd.read_csv(B+'out/apc/apc_functions_nci2023.csv')
for _,r in fn.iterrows(): add((r.est,r.lo,r.hi),f'{r.function} {r.measure} {r.sex} {r.etiology} {r.x}')
for sh in ['MASH_age_EAPC','All_etiologies_age_EAPC']:
    x=pd.read_excel(B+'out/robust/age_specific_mash.xlsx',sheet_name=sh)
    for _,r in x.iterrows():
        for p in ['1990_2023','1994_2023','1990_2019']: add((r['EAPC_'+p],r['lo_'+p],r['hi_'+p]),f'ageEAPC {p} {r.measure} {r.sex} {r.etiology} {r.age}')
        add((r.APC_local_drift,r.LD_lo,r.LD_hi),f'local drift {r.measure} {r.sex} {r.etiology} {r.age}')
        add((r.rate1990,r.rate1990_lo,r.rate1990_hi),f'agerate1990 {r.sex} {r.etiology} {r.age}'); add((r.rate2023,r.rate2023_lo,r.rate2023_hi),f'agerate2023 {r.sex} {r.etiology} {r.age}')
ew=pd.read_excel(B+'out/robust/equal_weight_asr.xlsx')
for _,r in ew.iterrows():
    for a,b_,c in [('EAPC_std','std_lo','std_hi'),('EAPC_equal_arith','equal_arith_lo','equal_arith_hi'),('EAPC_equal_geo','equal_geo_lo','equal_geo_hi'),('APC_net_drift','nd_lo','nd_hi')]: add((r[a],r[b_],r[c]),f'A5 {a} {r.measure} {r.sex} {r.etiology}')
rr=pd.read_excel(B+'out/robust/female_hcv_robustness.xlsx',sheet_name='R2_ratio_by_year')
for _,r in rr.iterrows(): add((r.ratio,r.ratio_lo,r.ratio_hi),f'HCV/HBV ratio {r.sex} {r.measure} {r.metric} {r.year}')
P=np.array([p[0] for p in pool])
sec=md[md.index('# Abstract'):md.index('# Declarations')]+md[md.index('# References'):]  # Abstract..Conclusions + tables/legends
pat=re.compile(r'(-?\d[\d,]*\.?\d*)%?(?: per year| per 100,000)?\s*\\?[\(\[]\s*(?:95% (?:UI|CI) |approximate 95% CI )?(-?\d[\d,]*\.?\d*)%?\s*(?:to|--|,)\s*(-?\d[\d,]*\.?\d*)%?\s*\\?[\)\]]')
out=[]
for m in pat.finditer(sec):
    s=[g for g in m.groups()]; v=np.array([float(g.replace(',','')) for g in s])
    d=np.array([len(g.split('.')[1]) if '.' in g else 0 for g in s]); f=10.0**d
    Rd=np.sign(P)*np.floor(np.abs(P)*f+0.5+1e-7)/f   # strict: round half away from zero to printed precision
    ok=np.all(np.abs(Rd-v)<1e-9,axis=1)
    if v[1]>v[2]: continue  # not a CI (e.g. years)
    if 1900<=v[0]<=2100: continue
    out.append(dict(text=m.group(0),matched=bool(ok.any()),source=pool[int(np.argmax(ok))][1] if ok.any() else '',ctx=sec[max(0,m.start()-70):m.start()].replace('\n',' ')))
R=pd.DataFrame(out); R.to_csv(B+'out/audit/triplets_all.csv',index=False); R[~R.matched].to_csv(B+'out/audit/triplets_unmatched.csv',index=False)
print('pool',len(pool),'triplets',len(R),'unmatched',(~R.matched).sum()); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',90)
print(R[~R.matched][['text','ctx']].to_string())
