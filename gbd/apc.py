"""Age-period-cohort analysis (Rosenberg et al. 2014 NCI APC Web Tool parameterisation, Poisson GLM, overdispersion via Pearson X2).
Inputs: counts by sex x etiology x 5-yr age x year; population by sex x age x year.
Periods: 5-yr groups; ages 20-24..90-94 (15 groups). Outputs out/apc/."""
import sys, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
B=''
CNT=sys.argv[1] if len(sys.argv)>1 else B+'gbd/raw/gbd2021_china_age/IHME-GBD_2021_DATA-2cc0e067-1.csv'
POP=sys.argv[2] if len(sys.argv)>2 else B+'gbd/raw/gbd2021_china_pop/pop.csv'
TAG=sys.argv[3] if len(sys.argv)>3 else 'gbd2021'
Y0=int(sys.argv[4]) if len(sys.argv)>4 else 1992   # first year of first 5-yr period
OUT=B+'out/apc/'; import os; os.makedirs(OUT,exist_ok=True)
ET={'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV','Liver cancer due to alcohol use':'Alcohol','Liver cancer due to NASH':'MASH','Liver cancer due to other causes':'Other'}
import os as _os
_A0,_A1=[int(x) for x in _os.environ.get('APC_AGES','20,94').split(',')]
AGES=[f'{a}-{a+4} years' for a in range(_A0,_A1+1,5)]
c=pd.read_csv(CNT); c=c[c.cause_name.isin(ET)&c.age_name.isin(AGES)&(c.metric_name=='Number')]
c['et']=c.cause_name.map(ET)
allc=c.groupby(['measure_name','sex_name','age_name','year'],as_index=False).val.sum(); allc['et']='All'
c=pd.concat([c[['measure_name','sex_name','age_name','year','et','val']],allc])
p=pd.read_csv(POP); p=p[p.age_name.isin(AGES)&p.sex_name.isin(['Male','Female'])][['sex_name','age_name','year','val']].rename(columns={'val':'pop'})
Y1=int(_os.environ.get('APC_Y1',c.year.max())); c=c[c.year<=Y1]; p=p[p.year<=Y1]; nP=(Y1-Y0+1)//5; Y0=Y1-5*nP+1
def per(y): return (y-Y0)//5
c=c[c.year>=Y0]; p=p[p.year>=Y0]
c['pi']=per(c.year); p['pi']=per(p.year)
c['ai']=c.age_name.map({a:i for i,a in enumerate(AGES)}); p['ai']=p.age_name.map({a:i for i,a in enumerate(AGES)})
A=len(AGES); P=nP; C=A+P-1
amid=np.array([_A0+2.5+5*i for i in range(A)]); pmid=np.array([Y0+2+5*j for j in range(P)])
cmid=pmid[0]-amid[-1]+5*np.arange(C)   # cohort midpoint birth years
def curv(K):
    x=np.arange(K)-(K-1)/2; V=np.vander(x,K,increasing=True); Q,_=np.linalg.qr(V); return Q[:,2:]   # orthonormal, orthogonal to const+linear
QA,QP,QC=curv(A),curv(P),curv(C)
ab=(A-1)/2; pb=(P-1)/2
def row(i,j,pdev=True):
    k=j-i+A-1
    r=[1,i-ab,j-pb]+list(QA[i])+(list(QP[j]) if pdev else [0]*(P-2))
    kk=k if 0<=k<C else None
    r+=list(QC[kk]) if kk is not None else [0]*(C-2)
    return np.array(r,float)
nA,nPd,nC=QA.shape[1],QP.shape[1],QC.shape[1]
iA=slice(3,3+nA); iPd=slice(3+nA,3+nA+nPd); iC=slice(3+nA+nPd,3+nA+nPd+nC)
res=[];curves=[]
z=stats.norm.ppf(.975)
for (m,s,e),g in c.groupby(['measure_name','sex_name','et']):
    d=g.groupby(['ai','pi'],as_index=False).val.sum().merge(p[p.sex_name==s].groupby(['ai','pi'],as_index=False)['pop'].sum(),on=['ai','pi'])
    X=np.vstack([row(int(i),int(j)) for i,j in zip(d.ai,d.pi)])
    fit=sm.GLM(d.val.values,X,family=sm.families.Poisson(),offset=np.log(d['pop'].values)).fit(scale='X2')
    b=fit.params; V=fit.cov_params()
    def est(L,scale=1):
        v=L@b; se=np.sqrt(L@V@L); return v*scale,se*scale
    def wald(idx):
        bb=b[idx]; VV=V[idx][:,idx]; W=bb@np.linalg.solve(VV,bb); df=len(bb); return stats.chi2.sf(W,df)
    pct=lambda v: 100*(np.exp(v)-1)
    L=np.zeros(len(b)); L[2]=1/5
    v,se=est(L); nd=(pct(v),pct(v-z*se),pct(v+z*se))
    # local drifts: slope of fitted log rate across periods within age i (per year)
    w=(np.arange(P)-pb)/np.sum((np.arange(P)-pb)**2)/5
    ld=[]
    for i in range(A):
        L=sum(w[j]*row(i,j) for j in range(P)); v,se=est(L); ld.append((pct(v),pct(v-z*se),pct(v+z*se)))
    # local drifts all equal to net drift <=> cohort curvature = 0
    # FIX 2026-11-10: local-drift test was wald(iC) = cohort test. Now a proper joint Wald test of
    # H0: local drift_i - net drift = 0 for all ages (contrast matrix D, df = rank, Rosenberg 2014 df = A-2)
    Lnd=np.zeros(len(b)); Lnd[2]=1/5
    D=np.vstack([sum(w[j]*row(i,j) for j in range(P))-Lnd for i in range(A)])
    dv=D@b; DV=D@V@D.T; rk=np.linalg.matrix_rank(DV,tol=1e-10*np.abs(DV).max())
    W=dv@np.linalg.pinv(DV,rcond=1e-10)@dv; p_ld=stats.chi2.sf(W,rk); ld_df=rk
    # Rosenberg-style: test deviations of local drifts from their own mean and linear age trend (df=A-2)
    xa=np.arange(A)-ab; H=np.column_stack([np.ones(A),xa]); M=np.eye(A)-H@np.linalg.pinv(H)
    dv2=M@dv; DV2=M@DV@M.T; rk2=np.linalg.matrix_rank(DV2,tol=1e-10*np.abs(DV2).max()); W2=dv2@np.linalg.pinv(DV2,rcond=1e-10)@dv2
    p_ld13=stats.chi2.sf(W2,rk2); ld_df13=rk2; p_ld_df13_sameW=stats.chi2.sf(W,A-2)
    _ldm=np.mean(dv)
    p_per=wald(np.r_[iPd]); p_coh=wald(np.r_[iC]); p_nd=2*stats.norm.sf(abs(b[2]/np.sqrt(V[2,2])))
    res.append(dict(measure=m,sex=s,etiology=e,net_drift=nd[0],nd_lo=nd[1],nd_hi=nd[2],p_netdrift=p_nd,p_localdrift_equal=p_ld,ld_df=ld_df,p_ld_alt=p_ld13,ld_df_alt=ld_df13,p_ld_sameW_dfA2=p_ld_df13_sameW,ld_mean_minus_nd=_ldm,p_period_dev=p_per,p_cohort_dev=p_coh,deviance=fit.deviance,df=fit.df_resid,scale=fit.scale))
    # period RR (ref = middle period), cohort RR (ref = middle cohort): drift*(x-x0)+deviations
    j0=P//2; k0=C//2
    for j in range(P):
        L=np.zeros(len(b)); L[2]=j-j0; L[iPd]=QP[j]-QP[j0]; v,se=est(L)
        curves.append((m,s,e,'period_RR',pmid[j],np.exp(v),np.exp(v-z*se),np.exp(v+z*se)))
    for k in range(C):
        L=np.zeros(len(b)); L[2]=k-k0; L[iC]=QC[k]-QC[k0]; v,se=est(L)
        curves.append((m,s,e,'cohort_RR',cmid[k],np.exp(v),np.exp(v-z*se),np.exp(v+z*se)))
    for i in range(A):
        curves.append((m,s,e,'local_drift',amid[i],*ld[i]))
        # longitudinal age curve: reference cohort k0, period deviations set to 0; per 100k
        j=k0+i-(A-1); L=np.zeros(len(b)); L[0]=1; L[1]=i-ab; L[2]=j-pb; L[iA]=QA[i]; L[iC]=QC[k0]
        v,se=est(L); curves.append((m,s,e,'long_age_rate',amid[i],1e5*np.exp(v),1e5*np.exp(v-z*se),1e5*np.exp(v+z*se)))
R=pd.DataFrame(res); K=pd.DataFrame(curves,columns=['measure','sex','etiology','function','x','est','lo','hi'])
R.to_csv(OUT+f'apc_summary_{TAG}.csv',index=False); K.to_csv(OUT+f'apc_functions_{TAG}.csv',index=False)
print(f'periods {pmid}, cohorts {cmid[0]}-{cmid[-1]}')
print(R[['measure','sex','etiology','net_drift','nd_lo','nd_hi','p_localdrift_equal','p_period_dev','p_cohort_dev','scale']].round(3).to_string())
