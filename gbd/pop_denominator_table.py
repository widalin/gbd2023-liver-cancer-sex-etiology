# Supplementary Table S4 part B: deviations of derived GBD 2023 population from GBD 2021 official and UN WPP 2024; outputs out/audit/pop_denominator_table.csv
import pandas as pd, numpy as np
B=''  # run from the package root
d=pd.read_csv(B+'gbd/raw/gbd2023_china_pop/pop_derived_allcause.csv')[['sex_name','age_name','year','val']]
g=pd.read_csv(B+'gbd/raw/gbd2021_china_pop/pop.csv')[['sex_name','age_name','year','val']]
w=pd.read_csv(B+'out/audit/pop_wpp2024_15plus.csv')
def band(a):
    if a=='95+ years': return '95+'
    lo=int(a.split('-')[0]); return '15–49' if lo<50 else ('50–79' if lo<80 else '80–94')
rows=[]
for name,ref,yr in [('GBD 2021 official',g,(1990,2021)),('UN WPP 2024',w,(1990,2023))]:
    m=d.merge(ref,on=['sex_name','age_name','year'],suffixes=('_d','_r')); m=m[m.year.between(*yr)&m.age_name.str.match(r'^(1[5-9]|[2-9]\d)')]
    m['band']=m.age_name.map(band); m['rel']=100*(m.val_d/m.val_r-1)
    for s in ['Female','Male']:
        x=m[m.sex_name==s]
        t=x.groupby('year')[['val_d','val_r']].sum(); r=100*(t.val_d/t.val_r-1)
        rows.append(dict(reference=name,years=f'{yr[0]}–{yr[1]}',sex=s,ages='≥15 total',median=r.median(),min=r.min(),max=r.max()))
        for b in ['15–49','50–79','80–94','95+']:
            y=x[x.band==b].rel; rows.append(dict(reference=name,years=f'{yr[0]}–{yr[1]}',sex=s,ages=b,median=y.median(),min=y.min(),max=y.max()))
T=pd.DataFrame(rows).round(1); T.to_csv(B+'out/audit/pop_denominator_table.csv',index=False); print(T.to_string())
