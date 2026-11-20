"""Das Gupta three-factor decomposition of change in counts (population size, age structure, age-specific rate). Ages 15-95+."""
import sys,numpy as np,pandas as pd,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B=''
CNT=sys.argv[1] if len(sys.argv)>1 else B+'gbd/raw/coauthor_A_china_age_etiology_2023/IHME-GBD_2023_DATA-fc723b63-1.csv'  # main analysis: GBD 2023 counts
POP=sys.argv[2] if len(sys.argv)>2 else B+'gbd/raw/gbd2023_china_pop/pop_derived_allcause.csv'  # main analysis: derived GBD 2023 population (all-cause deaths / rate)
TAG=sys.argv[3] if len(sys.argv)>3 else 'gbd2023'
ET={'Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV','Liver cancer due to alcohol use':'Alcohol','Liver cancer due to NASH':'MASH','Liver cancer due to other causes':'Other'}
c=pd.read_csv(CNT); c=c[c.cause_name.isin(ET)]
if 'metric_name' in c: c=c[c.metric_name=='Number']  # FIX 2026-11-10: GBD 2023 file has Number+Rate rows
c['et']=c.cause_name.map(ET)
a=c.groupby(['measure_name','sex_name','age_name','year'],as_index=False).val.sum(); a['et']='All'
c=pd.concat([c[['measure_name','sex_name','age_name','year','et','val']],a])
p=pd.read_csv(POP)[['sex_name','age_name','year','val']].rename(columns={'val':'pop'})
y0,y1=1990,c.year.max(); rows=[]
for (m,s,e),g in c.groupby(['measure_name','sex_name','et']):
    d={}
    for y in (y0,y1):
        x=g[g.year==y].merge(p[(p.sex_name==s)&(p.year==y)],on=['age_name']).sort_values('age_name')
        P=x['pop'].sum(); d[y]=(P,(x['pop']/P).values,(x.val/x['pop']).values)
    (A0,b0,c0),(A1,b1,c1)=d[y0],d[y1]
    eA=np.sum((A1-A0)*((b0*c0+b1*c1)/3+(b0*c1+b1*c0)/6))
    eB=np.sum((b1-b0)*((A0*c0+A1*c1)/3+(A0*c1+A1*c0)/6))
    eC=np.sum((c1-c0)*((A0*b0+A1*b1)/3+(A0*b1+A1*b0)/6))
    N0,N1=np.sum(A0*b0*c0),np.sum(A1*b1*c1)
    rows.append(dict(measure=m,sex=s,etiology=e,N_1990=N0,N_end=N1,change=N1-N0,population=eA,ageing=eB,epidemiological=eC,check=eA+eB+eC-(N1-N0),
        pct_population=100*eA/(N1-N0),pct_ageing=100*eB/(N1-N0),pct_epi=100*eC/(N1-N0)))
R=pd.DataFrame(rows); R.to_csv(B+f'out/apc/decomposition_{TAG}.csv',index=False)
print(R.drop(columns=['check']).round(0).to_string()); print('max check',R.check.abs().max())
fig,ax=plt.subplots(1,2,figsize=(13,5.2),sharey=True)
order=['All','HBV','HCV','Alcohol','MASH','Other']; cols={'population':'#7f93b8','ageing':'#c9935e','epidemiological':'#607868'}
for k,(m) in enumerate(['Incidence','Deaths']):
    a=ax[k]; ylab=[]; yi=0
    for e in order:
        for s in ['Male','Female']:
            r=R[(R.measure==m)&(R.sex==s)&(R.etiology==e)].iloc[0]; pos=0;neg=0
            for f in ['population','ageing','epidemiological']:
                v=r[f]/1000
                if v>=0: a.barh(yi,v,left=pos,color=cols[f],height=.75,label=f.capitalize() if yi==0 else None); pos+=v
                else: a.barh(yi,v,left=neg,color=cols[f],height=.75,label=f.capitalize() if yi==0 else None); neg+=v
            a.plot(r.change/1000,yi,'D',color='black',ms=5,label='Net change' if yi==0 else None)
            ylab.append(f"{e}, {s[0]}"); yi+=1
        yi+=0.6
    a.axvline(0,color='grey',lw=.8); a.set_title(f"{m}: change in counts 1990-{y1} (thousands)",loc='left',fontsize=11)
    a.set_yticks([i+0.6*(i//2) for i in range(len(ylab))]); a.set_yticklabels(ylab,fontsize=8.5); a.invert_yaxis()
    a.spines[['top','right']].set_visible(False)
ax[0].legend(fontsize=8.5,frameon=False,loc='upper right')
fig.suptitle(f"Decomposition of change in liver cancer incident cases and deaths by etiology and sex in China ({TAG.upper().replace('GBD','GBD ')}, ages 15+)",x=.01,ha='left',fontsize=12.5)
fig.text(.01,.005,"Das Gupta three-factor decomposition: population growth, population ageing, and change in age-specific rates.",fontsize=8,color='#555')
plt.tight_layout(rect=(0,.03,1,.94)); plt.savefig(B+f'out/apc/Fig6_decomposition_{TAG}.png',dpi=200)
