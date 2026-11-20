# Direct sex contrasts (review point 2). Conditional/model-based CIs.
import pandas as pd, numpy as np, statsmodels.formula.api as smf
from scipy import stats
d=pd.read_csv('out/main/master_dataset_long.csv')
d=d[(d.location=='China')&(d.age=='Age-standardized')&(d.metric=='Rate')&(d.sex!='Both')]
t2=pd.read_csv('out/main/table2_china.csv')
def aapc(m,s,e):
    r=t2[(t2.measure==m)&(t2.sex==s)&(t2.etiology==e)].iloc[0]
    return r.AAPC,(r.AAPC_hi-r.AAPC_lo)/3.92
def eapc_contrast(m,terms,weights):
    # 2026-11-11 (2nd review): calendar-year paired contrast. For each year t, D_t = sum_i w_i*ln(ASR_i,t);
    # OLS D_t ~ t (34 years); Newey-West HAC (Bartlett kernel, maxlags 3) on the single annual series.
    # The slope equals the same linear contrast of group-specific log-linear slopes, and shared year-to-year
    # fluctuations are handled by pairing on calendar year. (Earlier version stacked the four series and
    # applied HAC over row order - wrong; archived gbd/sex_contrasts_pre_paired.py.)
    import statsmodels.api as sm
    D=None
    for (s,e),w in zip(terms,weights):
        x=d[(d.measure==m)&(d.sex==s)&(d.etiology==e)].set_index('year').val.sort_index()
        D=w*np.log(x) if D is None else D+w*np.log(x)
    t=D.index.values-1990
    f=sm.OLS(D.values,sm.add_constant(t)).fit(cov_type='HAC',cov_kwds={'maxlags':3})
    return f.params[1],f.bse[1]
out=[]
CONTRASTS=[('Overall: male − female',[('Male','All'),('Female','All')],[1,-1]),
   ('MASH: male − female',[('Male','MASH'),('Female','MASH')],[1,-1]),
   ('HBV: male − female',[('Male','HBV'),('Female','HBV')],[1,-1]),
   ('HCV: male − female',[('Male','HCV'),('Female','HCV')],[1,-1]),
   ('Alcohol: male − female',[('Male','Alcohol'),('Female','Alcohol')],[1,-1]),
   ('Women: HCV − HBV',[('Female','HCV'),('Female','HBV')],[1,-1]),
   ('Men: HCV − HBV',[('Male','HCV'),('Male','HBV')],[1,-1]),
   ('(HCV − HBV) women − men',[('Female','HCV'),('Female','HBV'),('Male','HCV'),('Male','HBV')],[1,-1,-1,1])]
for m in ['Incidence','Deaths']:
    for lab,terms,w in CONTRASTS:
        a=[aapc(m,s,e) for s,e in terms]
        ae=sum(wi*v for wi,(v,_) in zip(w,a)); ase=np.sqrt(sum(s**2 for _,s in a))
        eb,es=eapc_contrast(m,terms,w)
        # slope difference on log scale -> percentage-point approx via 100*b
        pz=2*stats.norm.sf(abs(ae/ase)); pe=2*stats.norm.sf(abs(eb/es))
        out.append(dict(measure=m,contrast=lab,AAPC_diff=round(ae,2),AAPC_lo=round(ae-1.96*ase,2),AAPC_hi=round(ae+1.96*ase,2),AAPC_p=pz,
            slope_diff_pct=round(100*eb,2),slope_lo=round(100*(eb-1.96*es),2),slope_hi=round(100*(eb+1.96*es),2),slope_p=pe))
o=pd.DataFrame(out); o.to_csv('out/contrast/sex_contrasts.csv',index=False)
pd.set_option('display.width',250); print(o.round(4).to_string())
