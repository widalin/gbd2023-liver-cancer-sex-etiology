# 2026-11-22 review 16 (M1): female HCV − HBV mortality AAPC difference in Japan and Korea,
# same method as China in Supplementary Table S2 (gbd/sex_contrasts.py): AAPCs from NCI Joinpoint 6.1.0,
# SE = (upper − lower 95% CI)/3.92, AAPCs treated as independent, normal approximation.
# Positive value = HBV-related ASDR declined faster than HCV-related ASDR.
import pandas as pd, numpy as np
from scipy import stats
import os
F=[f for f in ['out/figs_final/_mainfigs_side/fig4_JKUS_deaths_AAPC_by_sex.csv','out/main/fig4_JKUS_deaths_AAPC_by_sex.csv'] if os.path.exists(f)][0]  # written by main_figs.py
d=pd.read_csv(F)
rows=[]
for loc in ['China','Japan','Korea','USA']:
    g=d[(d.location==loc)&(d.sex=='Female')].set_index('cause')
    h,c=g.loc['HBV'],g.loc['HCV']
    se=np.hypot((h.AAPC_hi-h.AAPC_lo)/3.92,(c.AAPC_hi-c.AAPC_lo)/3.92)
    diff=c.AAPC-h.AAPC; z=diff/se; p=2*stats.norm.sf(abs(z))
    rows.append(dict(location=loc,sex='Female',measure='Deaths',HBV_AAPC=h.AAPC,HCV_AAPC=c.AAPC,diff_HCV_minus_HBV=diff,lo=diff-1.96*se,hi=diff+1.96*se,z=z,P=p))
os.makedirs('out/contrast',exist_ok=True)
o=pd.DataFrame(rows); o.to_csv('out/contrast/jk_female_hbv_hcv_deaths_aapc_diff.csv',index=False)
print(o.round(4).to_string())
