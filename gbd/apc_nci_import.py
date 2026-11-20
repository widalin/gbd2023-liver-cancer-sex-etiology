# Import NCI APC Web Tool (analysistools.cancer.gov/apc, Default reference) workbooks -> in-house schema (TAG nci2023)
import pandas as pd,shutil,os
B=''; D='data/nci_apc_exports/'  # the 24 NCI APC Web Tool exports under their original names (see Table S2, Source column)
MAP={'APC_Analysis_13_9_4-24a067e6':'Female_HBV_incidence','APC_Analysis_13_11_7-1e106dce':'Female_HCV_incidence','APC_Analysis_13_11_41-e6abd267':'Male_MASH_incidence','APC_Analysis_13_12_7-333aaf39':'Male_HBV_incidence',
'APC_Analysis_13_28_0-bccd9ecc':'Female_Other_deaths','APC_Analysis_13_28_14-d8710f4d':'Female_Other_incidence','APC_Analysis_13_28_29-e0717946':'Male_Alcohol_deaths','APC_Analysis_13_28_45-0236aa9d':'Male_Alcohol_incidence',
'APC_Analysis_13_28_1-8210dda2':'Female_HBV_deaths','APC_Analysis_13_28_19-69d6e29e':'Female_HCV_deaths','APC_Analysis_13_28_34-a6b01a08':'Female_MASH_deaths','APC_Analysis_13_28_50-18e9bd6d':'Female_MASH_incidence',
'APC_Analysis_13_28_3-a1af91f3':'Female_Alcohol_incidence','APC_Analysis_13_28_22-d70c4335':'Female_Alcohol_deaths','APC_Analysis_13_28_38-a858fa90':'Female_All_incidence','APC_Analysis_13_28_57-212773fc':'Female_All_deaths',
'APC_Analysis_13_28_5-2a6dd791':'Male_HCV_incidence','APC_Analysis_13_28_29-cc4c3227':'Male_MASH_deaths','APC_Analysis_13_28_48-59dac936':'Male_Other_deaths','APC_Analysis_13_29_7-bb34d5d2':'Male_Other_incidence',
'APC_Analysis_13_28_19-3f87e4ae':'Male_All_deaths','APC_Analysis_13_28_38-0aac7928':'Male_All_incidence','APC_Analysis_13_28_54-ef9f92d7':'Male_HBV_deaths','APC_Analysis_13_29_11-32ea8eb5':'Male_HCV_deaths'}
os.makedirs(B+'out/apc_nci/workbooks',exist_ok=True); S=[];K=[]
FN={'LocalDrifts':'local_drift','LongAge':'long_age_rate','PeriodRR':'period_RR','CohortRR':'cohort_RR'}
for f,n in MAP.items():
    src=D+f+'.xlsx'; shutil.copy(src,B+f'out/apc_nci/workbooks/{n}.xlsx')
    x=pd.read_excel(src,sheet_name=None,engine='calamine',header=None)
    sex,e,m=n.split('_'); m={'incidence':'Incidence','deaths':'Deaths'}[m]
    nd=x['NetDrift'].iloc[1]; w=x['Waldtests'].set_index(0)
    g=lambda k:(float(w.loc[k,1]),int(w.loc[k,2]),float(w.loc[k,3]))
    S.append(dict(measure=m,sex=sex,etiology=e,net_drift=nd[0],nd_lo=nd[1],nd_hi=nd[2],
      p_netdrift=g('NetDrift = 0')[2],X2_localdrift=g('All Local Drifts = Net Drift')[0],ld_df=g('All Local Drifts = Net Drift')[1],p_localdrift_equal=g('All Local Drifts = Net Drift')[2],
      p_period_dev=g('All Period Deviations = 0')[2],p_cohort_dev=g('All Cohort Deviations = 0')[2],source=f+'.xlsx'))
    for sh,fn in FN.items():
        t=x[sh].iloc[1:,:4].astype(float); t.columns=['x','est','lo','hi']; t.insert(0,'function',fn); t.insert(0,'etiology',e); t.insert(0,'sex',sex); t.insert(0,'measure',m); K.append(t)
S=pd.DataFrame(S).sort_values(['measure','sex','etiology']); K=pd.concat(K)
S.to_csv(B+'out/apc/apc_summary_nci2023.csv',index=False); K.to_csv(B+'out/apc/apc_functions_nci2023.csv',index=False)
h=pd.read_csv(B+'out/apc/apc_summary_gbd2023.csv')[['measure','sex','etiology','net_drift','nd_lo','nd_hi','p_localdrift_equal']]
c=S.merge(h,on=['measure','sex','etiology'],suffixes=('_nci','_inhouse')); c['diff']=c.net_drift_nci-c.net_drift_inhouse
c['sig_nci']=(c.nd_hi_nci<0)|(c.nd_lo_nci>0); c['sig_in']=(c.nd_hi_inhouse<0)|(c.nd_lo_inhouse>0)
c.to_csv(B+'out/apc_nci/nci_vs_inhouse_netdrift.csv',index=False)
pd.set_option('display.width',250); print(c[['measure','sex','etiology','net_drift_nci','nd_lo_nci','nd_hi_nci','net_drift_inhouse','diff','sig_nci','sig_in','ld_df','p_localdrift_equal_nci']].round(3).to_string())
