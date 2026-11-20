"""Verify NCI APC workbook -> series mapping (2026-11-11). NCI workbooks hold no input data, so each workbook's fitted
cross-sectional age curve (reference period 2004-08) and fitted temporal trend (reference age 55-59) are matched against
Poisson APC fitted rates (identifiable) of all 24 input files (out/apc_benchmark/all24). A workbook is confirmed if its mapped series is the
best match by log-RMSE and clearly separated from the runner-up."""
import pandas as pd, numpy as np, glob, os, re
B=''; D='data/nci_apc_exports/'  # the 24 NCI APC Web Tool exports under their original names (see Table S2, Source column)
exec(re.search(r"MAP=\{.*?\}\n", open(B+'gbd/apc_nci_import.py').read(), re.S).group(0))
import statsmodels.api as sm
def fitted(Y,N):
    """Poisson APC fitted rates (identifiable), per 100,000; ages x periods."""
    A,P=Y.shape; ai,pi=np.meshgrid(range(A),range(P),indexing='ij'); ci=pi-ai+A-1
    X=np.column_stack([np.ones(A*P)]+[(ai.ravel()==k)*1. for k in range(1,A)]+[(pi.ravel()==k)*1. for k in range(1,P)]+[(ci.ravel()==k)*1. for k in range(1,A+P-1)])
    m=sm.GLM(Y.ravel(),X,family=sm.families.Poisson(),offset=np.log(N.ravel())).fit()
    return 1e5*(m.fittedvalues/N.ravel()).reshape(A,P)
obs={}
for f in sorted(glob.glob(B+'out/apc_benchmark/all24/*.csv')):
    rows=[l for l in open(f).read().splitlines() if re.match(r'^[\d.]+,',l)]
    a=np.array([[float(v) for v in r.split(',')] for r in rows]); obs[os.path.basename(f)[:-4]]=fitted(a[:,0::2],a[:,1::2])
out=[]
for wb,mapped in MAP.items():
    x=pd.read_excel(D+wb+'.xlsx',sheet_name=None,engine='calamine',header=None)
    ca=x['CrossAge'].iloc[1:,1].astype(float).values; ft=x['FittedTemporalTrends'].iloc[1:,1].astype(float).values
    sc={k:np.sqrt(np.mean(np.r_[(np.log(ca)-np.log(o[:,2]))**2,(np.log(ft)-np.log(o[7,:]))**2])) for k,o in obs.items()}
    r=sorted(sc.items(),key=lambda t:t[1])
    out.append(dict(workbook=wb,mapped=mapped,best=r[0][0],rmse_best=r[0][1],second=r[1][0],rmse_second=r[1][1],rmse_mapped=sc[mapped],ok=r[0][0]==mapped))
o=pd.DataFrame(out); o.to_csv(B+'out/apc_nci/mapping_check.csv',index=False)
pd.set_option('display.width',250); print(o.round(3).to_string()); print('confirmed',o.ok.sum(),'/',len(o),' unique mapped',o.mapped.nunique())
