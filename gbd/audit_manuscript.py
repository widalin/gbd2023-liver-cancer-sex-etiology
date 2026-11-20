"""Full numeric audit of manuscript_v2.5.md (Abstract, Results, Table 1, Table 2, Table S1) against source data.
Writes out/audit/audit_v2.4.csv (one row per checked number) and prints mismatches."""
import re, sys, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
B=''
sys.path.insert(0,B+'gbd'); import nci
md=open(B+'manuscript/manuscript_v2.9.md').read().replace('\u2212','-')
M=pd.read_csv(B+'gbd/GBD2023_liver_cancer_etiology_Global_China_1990-2023.csv')
CN={'All liver cancer':'Liver cancer','HBV':'Liver cancer due to hepatitis B','HCV':'Liver cancer due to hepatitis C','Alcohol':'Liver cancer due to alcohol use',
    'MASH':'Liver cancer due to NASH','Other causes':'Liver cancer due to other causes','Other':'Liver cancer due to other causes','Hepatoblastoma':'Hepatoblastoma'}
def g(loc,meas,sex,cause,metric,year,age='All ages'):
    if metric=='ASR': metric,age='Rate','Age-standardized'
    x=M[(M.location_name==loc)&(M.measure_name==meas)&(M.sex_name==sex)&(M.cause_name==CN.get(cause,cause))&(M.metric_name==metric)&(M.year==year)&(M.age_name==age)]
    assert len(x)==1,(loc,meas,sex,cause,metric,year,len(x)); return x.iloc[0]
def eapc(loc,meas,sex,cause,y0=1990,y1=2023):
    x=M[(M.location_name==loc)&(M.measure_name==meas)&(M.sex_name==sex)&(M.cause_name==CN[cause])&(M.metric_name=='Rate')&(M.age_name=='Age-standardized')&(M.year.between(y0,y1))].sort_values('year')
    f=sm.OLS(np.log(x.val.values),sm.add_constant(x.year.values)).fit(); b,se=f.params[1],f.bse[1]; t=stats.t.ppf(.975,len(x)-2)
    return [100*(np.exp(v)-1) for v in (b,b-t*se,b+t*se)]
NC={'All liver cancer':'All','Other causes':'Other'}
def aapc(loc,meas,sex,cause): k,A,s,f,p=nci.result(loc,meas,sex,NC.get(cause,cause)); return list(A)
rows=[]
def chk(where,claim,ms,src,tol=None):
    ms=float(str(ms).replace(',','')); d=abs(ms-src)
    if tol is None:
        s=str(ms); dec=len(str(ms).split('.')[1]) if '.' in str(ms) else 0; tol=0.5*10**(-dec)+1e-9
    rows.append(dict(where=where,claim=claim,manuscript=ms,source=round(src,4),ok=d<=tol))
def num(s): return float(s.replace(',',''))
# ---------- Table 1 ----------
t1=re.findall(r'^\s+(Incidence|Deaths)\s+(Both|Male|Female)\s+(All liver cancer|HBV|HCV|Alcohol|MASH|Other causes|Hepatoblastoma)\s+([\d,]+) \(([\d,]+)--([\d,]+)\)\s+([\d.]+) \(([\d.]+)--([\d.]+)\)\s+([\d,]+) \(([\d,]+)--([\d,]+)\)\s+([\d.]+) \(([\d.]+)--([\d.]+)\)\s+([\d.]+)\s+(-?[\d.]+)',md,re.M)
for r in t1:
    m,s,e=r[:3]; w=f'Table1 {m} {s} {e}'
    for y,(i0) in ((1990,3),(2023,9)):
        n=g('China',m,s,e,'Number',y); a=g('China',m,s,e,'ASR',y)
        for lab,v,src in (('N',r[i0],n.val),('N lo',r[i0+1],n.lower),('N hi',r[i0+2],n.upper),('ASR',r[i0+3],a.val),('ASR lo',r[i0+4],a.lower),('ASR hi',r[i0+5],a.upper)):
            chk(w,f'{y} {lab}',num(v),src,tol=0.5+1e-9 if lab.startswith('N') else None)
    tot=g('China',m,s,'All liver cancer','Number',2023).val
    chk(w,'share 2023 %',num(r[15]),100*g('China',m,s,e,'Number',2023).val/tot)
    chk(w,'ASR change %',num(r[16]),100*(g('China',m,s,e,'ASR',2023).val/g('China',m,s,e,'ASR',1990).val-1))
print('Table1 rows parsed',len(t1))
# ---------- Table 2 ----------
t2=re.findall(r'^\s+(All liver cancer|HBV|HCV|Alcohol|MASH|Other|Hepatoblastoma)\s+(-?[\d.]+) \((-?[\d.]+) to (-?[\d.]+)\)\s+(-?[\d.]+) \((-?[\d.]+) to (-?[\d.]+)\)\s+(-?[\d.]+) \((-?[\d.]+) to (-?[\d.]+)\)\s+(-?[\d.]+) \((-?[\d.]+) to (-?[\d.]+)\)',md[md.index('**Table 2.**'):md.index('# Supplementary tables')],re.M)
for r in t2:
    e=r[0]
    for j,(m,s) in enumerate([('Incidence','Male'),('Incidence','Female'),('Deaths','Male'),('Deaths','Female')]):
        A=aapc('China',m,s,e)
        for k,lab in enumerate(['AAPC','lo','hi']): chk(f'Table2 {e} {m} {s}',lab,num(r[1+3*j+k]),A[k])
print('Table2 rows parsed',len(t2))
# ---------- Table S1: all columns vs NCI (main + batch 3 sensitivity) and EAPC ----------
S3=pd.read_csv(B+'out/nci_joinpoint/rerun3_output/China_24_sensitivity_NCI_wide.csv').set_index('name')
P='(\\S+) \\((\\S+), (\\S+)\\)[\\s|]+'
s1=re.findall(r'^[\s|]+(Incidence|Deaths) (Male|Female) (MASH|Alcohol|HBV|HCV)[\s|]+'+P*4+r'(\d{4})--(\d{4}) APC (\S+) \((\S+), (\S+)\)',md,re.M)
for r in s1:
    m,s,e=r[:3]; n=f'China_{m}_{s}_{e}'; k_,A,segs,_,_=nci.result('China',m,s,e); x=S3.loc[n]; E=eapc('China',m,s,e); sg=segs[-1]
    src=[A,[x.AAPC_B,x.lo_B,x.hi_B],[x.AAPC_A,x.lo_A,x.hi_A],E]
    for j,lab in enumerate(['Main NCI','Max3 NCI','1990-2019 NCI','EAPC']):
        for k,kk in enumerate(['','lo','hi']): chk(f'S1 {m} {s} {e}',f'{lab} {kk}',num(r[3+3*j+k]),src[j][k])
    chk(f'S1 {m} {s} {e}','last seg start',num(r[15]),sg[0],0); chk(f'S1 {m} {s} {e}','last seg end',num(r[16]),sg[1],0)
    for k,kk in enumerate(['APC','lo','hi']): chk(f'S1 {m} {s} {e}',f'last {kk}',num(r[17+k]),sg[2+k])
print('S1 rows parsed',len(s1))
# ---------- Abstract + Results text claims ----------
ap=pd.read_csv(B+'out/apc/apc_summary_nci2023.csv'); dc=pd.read_csv(B+'out/apc/decomposition_gbd2023.csv')
def nd(m,s,e): x=ap[(ap.measure==m)&(ap.sex==s)&(ap.etiology==e)].iloc[0]; return [x.net_drift,x.nd_lo,x.nd_hi]
def share(m,s,e,y=2023): return 100*g('China',m,s,e,'Number',y).val/g('China',m,s,'All liver cancer','Number',y).val
def mfr(y): return g('China','Incidence','Male','All liver cancer','ASR',y).val/g('China','Incidence','Female','All liver cancer','ASR',y).val
T=[('Abstract','cases 2023',165113,g('China','Incidence','Both','All liver cancer','Number',2023).val,0.5),
 ('Abstract','cases lo',143144,g('China','Incidence','Both','All liver cancer','Number',2023).lower,0.5),
 ('Abstract','cases hi',190095,g('China','Incidence','Both','All liver cancer','Number',2023).upper,0.5),
 ('Abstract','deaths 2023',143638,g('China','Deaths','Both','All liver cancer','Number',2023).val,0.5),
 ('Abstract','deaths lo',126474,g('China','Deaths','Both','All liver cancer','Number',2023).lower,0.5),
 ('Abstract','deaths hi',163971,g('China','Deaths','Both','All liver cancer','Number',2023).upper,0.5),
 ('Abstract','M:F incidence ratio 1990',2.21,mfr(1990),None),('Abstract','M:F 2023',3.20,mfr(2023),None),
 ('Abstract','F HBV inc AAPC',-2.41,aapc('China','Incidence','Female','HBV')[0],None),('Abstract','F HCV inc AAPC',-1.56,aapc('China','Incidence','Female','HCV')[0],None),
 ('Abstract','F HCV share cases',39.0,share('Incidence','Female','HCV'),None),('Abstract','M HBV share cases',72.7,share('Incidence','Male','HBV'),None),
 ('Abstract','M MASH ASR 1990',0.48,g('China','Incidence','Male','MASH','ASR',1990).val,None),('Abstract','M MASH ASR 2023',0.53,g('China','Incidence','Male','MASH','ASR',2023).val,None),
 ('Abstract','M MASH net drift',-0.20,nd('Incidence','Male','MASH')[0],None),('Abstract','nd lo',-0.38,nd('Incidence','Male','MASH')[1],None),('Abstract','nd hi',-0.01,nd('Incidence','Male','MASH')[2],None),
 ('Abstract','M MASH cases 1990',2145,g('China','Incidence','Male','MASH','Number',1990).val,0.5),('Abstract','M MASH cases 2023',5795,g('China','Incidence','Male','MASH','Number',2023).val,0.5),
 ('Results','Both HBV share inc',62.6,share('Incidence','Both','HBV'),None),('Results','Both HBV share deaths',61.1,share('Deaths','Both','HBV'),None),
 ('Results','M HBV cases',89910,g('China','Incidence','Male','HBV','Number',2023).val,0.5),('Results','M HBV deaths share',71.8,share('Deaths','Male','HBV'),None),
 ('Results','M HBV deaths',75800,g('China','Deaths','Male','HBV','Number',2023).val,0.5),
 ('Results','M HCV share',12.0,share('Incidence','Male','HCV'),None),('Results','M alcohol share',6.7,share('Incidence','Male','Alcohol'),None),('Results','M MASH share',4.7,share('Incidence','Male','MASH'),None),
 ('Results','F HCV cases',16177,g('China','Incidence','Female','HCV','Number',2023).val,0.5),('Results','F all cases',41449,g('China','Incidence','Female','All liver cancer','Number',2023).val,0.5),
 ('Results','F HCV death share',41.1,share('Deaths','Female','HCV'),None),('Results','F HCV deaths',15625,g('China','Deaths','Female','HCV','Number',2023).val,0.5),
 ('Results','F all deaths',38009,g('China','Deaths','Female','All liver cancer','Number',2023).val,0.5),
 ('Results','F HBV inc share',32.6,share('Incidence','Female','HBV'),None),('Results','F HBV death share',31.5,share('Deaths','Female','HBV'),None),
 ('Results','F HBV death share 1990',39.6,share('Deaths','Female','HBV',1990),None),('Results','F HCV death share 1990',32.5,share('Deaths','Female','HCV',1990),None),
 ('Results','F MASH inc share',10.6,share('Incidence','Female','MASH'),None),
 ('Results','F HBV ASIR 1990',2.61,g('China','Incidence','Female','HBV','ASR',1990).val,None),('Results','F HBV ASIR 2023',1.18,g('China','Incidence','Female','HBV','ASR',2023).val,None),
 ('Results','F HCV ASIR 1990',2.19,g('China','Incidence','Female','HCV','ASR',1990).val,None),('Results','F HCV ASIR 2023',1.33,g('China','Incidence','Female','HCV','ASR',2023).val,None),
 ('Results','F HCV ASDR 2023',1.28,g('China','Deaths','Female','HCV','ASR',2023).val,None),('Results','F HBV ASDR 2023',1.02,g('China','Deaths','Female','HBV','ASR',2023).val,None),
 ('Results','F EAPC inc HBV',-2.88,eapc('China','Incidence','Female','HBV')[0],None),('Results','F EAPC inc HCV',-1.93,eapc('China','Incidence','Female','HCV')[0],None),
 ('Results','M HBV ASIR 1990',11.24,g('China','Incidence','Male','HBV','ASR',1990).val,None),('Results','M HBV ASIR 2023',8.34,g('China','Incidence','Male','HBV','ASR',2023).val,None),
 ('Results','M alcohol EAPC',-0.27,eapc('China','Incidence','Male','Alcohol')[0],None),
 ('Results','M MASH max3',0.04,S3.loc['China_Incidence_Male_MASH','AAPC_B'],None),('Results','M MASH max3 lo',-0.14,S3.loc['China_Incidence_Male_MASH','lo_B'],None),('Results','M MASH max3 hi',0.23,S3.loc['China_Incidence_Male_MASH','hi_B'],None),('Results','M MASH to2019',-0.07,S3.loc['China_Incidence_Male_MASH','AAPC_A'],None),('Results','M MASH to2019 lo',-0.10,S3.loc['China_Incidence_Male_MASH','lo_A'],None),('Results','M MASH to2019 hi',-0.02,S3.loc['China_Incidence_Male_MASH','hi_A'],None),('Results','M alc to2019',-0.54,S3.loc['China_Incidence_Male_Alcohol','AAPC_A'],None),
 ('Results','M MASH EAPC',0.07,eapc('China','Incidence','Male','MASH')[0],None),('Results','M MASH EAPC lo',-0.10,eapc('China','Incidence','Male','MASH')[1],None),('Results','M MASH EAPC hi',0.24,eapc('China','Incidence','Male','MASH')[2],None),
 ('Results','M MASH ASIR 2015',0.54,g('China','Incidence','Male','MASH','ASR',2015).val,None),('Results','M MASH ASIR 2019',0.48,g('China','Incidence','Male','MASH','ASR',2019).val,None),
 ('Results','M MASH ASDR 1990',0.47,g('China','Deaths','Male','MASH','ASR',1990).val,None),('Results','M MASH ASDR 2023',0.47,g('China','Deaths','Male','MASH','ASR',2023).val,None),
 ('Results','Global F HCV ASIR 1990',1.54,g('Global','Incidence','Female','HCV','ASR',1990).val,None),('Results','Global F HBV ASIR 1990',1.07,g('Global','Incidence','Female','HBV','ASR',1990).val,None),
 ('Results','Global F HCV ASIR 2023',1.55,g('Global','Incidence','Female','HCV','ASR',2023).val,None),('Results','Global F HBV ASIR 2023',0.83,g('Global','Incidence','Female','HBV','ASR',2023).val,None),
]
for m,s,e,v in [('Incidence','Female','HBV',[-3.15,-3.29,-3.01]),('Incidence','Female','HCV',[-2.77,-3.06,-2.48]),('Incidence','Female','Alcohol',[-2.22,-2.43,-2.00]),('Incidence','Female','MASH',[-2.02,-2.17,-1.87]),('Incidence','Female','Other',[-2.85,-3.00,-2.70]),
                ('Incidence','Male','HBV',[-1.45,-1.65,-1.25]),('Incidence','Male','HCV',[-1.07,-1.37,-0.78]),('Incidence','Male','Alcohol',[-0.42,-0.71,-0.12]),('Incidence','Male','MASH',[-0.20,-0.38,-0.01]),('Deaths','Male','MASH',[-0.55,-0.74,-0.36])]:
    for k,lab in enumerate(['nd','lo','hi']): T.append(('Results APC',f'{m} {s} {e} {lab}',v[k],nd(m,s,e)[k],None))
x=ap[ap.etiology!='All']; rows_extra=[('Results APC','max p period dev <0.05',x.p_period_dev.max()),('Results APC','max p cohort dev <0.05',x.p_cohort_dev.max()),('Results APC','max p local drift df15 <0.05',x.p_localdrift_equal.max()),('Results APC','max p local drift df13 <0.05',x.p_ld_alt.max())]
def d(m,s,e,c): return dc[(dc.measure==m)&(dc.sex==s)&(dc.etiology==e)].iloc[0][c]
for lab,v,src in [('M all change',54151,d('Incidence','Male','All','change')),('F all change',11419,d('Incidence','Female','All','change')),
    ('M ageing',51142,d('Incidence','Male','All','ageing')),('F ageing',23514,d('Incidence','Female','All','ageing')),('M pop',28881,d('Incidence','Male','All','population')),('F pop',12061,d('Incidence','Female','All','population')),
    ('M rate',-25872,d('Incidence','Male','All','epidemiological')),('F rate',-24156,d('Incidence','Female','All','epidemiological')),('M MASH rate',231,d('Incidence','Male','MASH','epidemiological')),
    ('M MASH change',3651,d('Incidence','Male','MASH','change')),('M MASH pop',1124,d('Incidence','Male','MASH','population')),('M MASH ageing',2296,d('Incidence','Male','MASH','ageing')),
    ('F HCV change',6900,d('Incidence','Female','HCV','change')),('F HBV change',1192,d('Incidence','Female','HBV','change')),('F HCV rate',-7055,d('Incidence','Female','HCV','epidemiological')),('F HBV rate',-11141,d('Incidence','Female','HBV','epidemiological'))]:
    T.append(('Results decomp',lab,v,src,0.5))
for w,c,ms,src,tol in T: chk(w,c,ms,src,tol)
for w,c,v in rows_extra: rows.append(dict(where=w,claim=c,manuscript=float(c.split('<')[1]),source=round(v,4),ok=v<float(c.split('<')[1])))
for w,c,ms,tol in [('Results text','F MASH inc AAPC hi',-1.13,None),('Results text','F alcohol inc AAPC hi',-1.51,None)]: chk(w,c,ms,aapc('China','Incidence','Female','MASH' if 'MASH' in c else 'Alcohol')[2],tol)
R=pd.DataFrame(rows); import os; os.makedirs(B+'out/audit',exist_ok=True); R.to_csv(B+'out/audit/audit_v2.4.csv',index=False)
print('checked',len(R),'mismatches',(~R.ok).sum()); pd.set_option('display.width',200); print(R[~R.ok].to_string())

# ---- 2026-11-10 21:30 additions: numbers added in Results rewrite / Discussion pass
ew=pd.read_excel(B+'out/robust/equal_weight_asr.xlsx'); ewr=lambda m,s,e: ew[(ew.measure==m)&(ew.sex==s)&(ew.etiology==e)].iloc[0]
a5=ewr('Incidence','Male','MASH'); a5d=ewr('Deaths','Male','MASH')
rr=pd.read_excel(B+'out/robust/female_hcv_robustness.xlsx',sheet_name='R2_ratio_by_year'); rq=lambda m: rr[(rr.sex=='Female')&(rr.measure==m)&(rr.metric=='ASR')&(rr.year==2023)].iloc[0] if 'metric' in rr else None
ags=pd.read_excel(B+'out/robust/age_specific_mash.xlsx',sheet_name='All_etiologies_age_EAPC')
def ld(sex,e,age): return ags[(ags.measure=='Incidence')&(ags.sex==sex)&(ags.etiology==e)&(ags.age==age)].iloc[0].APC_local_drift
X=[('A5','std EAPC',0.04,a5.EAPC_std),('A5','geo',-0.26,a5.EAPC_equal_geo),('A5','geo lo',-0.43,a5.equal_geo_lo),('A5','geo hi',-0.09,a5.equal_geo_hi),
   ('A5','arith',0.41,a5.EAPC_equal_arith),('A5','arith lo',0.19,a5.equal_arith_lo),('A5','arith hi',0.63,a5.equal_arith_hi),
   ('A5','share 60+ %',64,a5.share_ASR_from_60plus_2023),('A5 deaths','std',-0.21,a5d.EAPC_std),('A5 deaths','geo',-0.63,a5d.EAPC_equal_geo),
   ('Results','M HCV death AAPC',-0.91,aapc('China','Deaths','Male','HCV')[0]),('Results','M other death AAPC',-1.02,aapc('China','Deaths','Male','Other')[0]),
   ('Results','M MASH death AAPC',0.12,aapc('China','Deaths','Male','MASH')[0]),
   ('Results','F MASH LD min 80+',0.38,min(ld('Female','MASH',a) for a in ags.age.unique() if str(a)[:2] in ('80','85','90'))),
   ('Results','F MASH LD max 80+',1.68,max(ld('Female','MASH',a) for a in ags.age.unique() if str(a)[:2] in ('80','85','90'))),
]
for w,c,ms,src in X: chk(w,c,ms,src,0.5 if c=='share 60+ %' else None)
R=pd.DataFrame(rows); R.to_csv(B+'out/audit/audit_v2.4.csv',index=False)
print('checked',len(R),'mismatches',(~R.ok).sum()); print(R[~R.ok].to_string())
