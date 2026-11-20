import pandas as pd
d=pd.read_csv('gbd/GBD2023_liver_cancer_etiology_Global_China_1990-2023.csv')
d=d[d.location_name=='China']
causes={'Liver cancer':'All liver cancer','Liver cancer due to hepatitis B':'HBV','Liver cancer due to hepatitis C':'HCV',
 'Liver cancer due to alcohol use':'Alcohol','Liver cancer due to NASH':'MASH','Liver cancer due to other causes':'Other causes','Hepatoblastoma':'Hepatoblastoma'}
meas={'Incidence':'Incidence','Deaths':'Deaths','DALYs (Disability-Adjusted Life Years)':'DALYs'}
def g(m,c,s,y,age,met):
    r=d[(d.measure_name==m)&(d.cause_name==c)&(d.sex_name==s)&(d.year==y)&(d.age_name==age)&(d.metric_name==met)]
    assert len(r)==1,(m,c,s,y,age,met,len(r)); r=r.iloc[0]; return r.val,r.lower,r.upper
def fn(v,dec=0):
    f=f'{{:,.{dec}f}}'; return f'{f.format(v[0])} ({f.format(v[1])}–{f.format(v[2])})'
rows=[];raw=[]
for m,ml in meas.items():
  for s in ['Both','Male','Female']:
    tot23=g(m,'Liver cancer',s,2023,'All ages','Number')[0]
    for c,cl in causes.items():
      n90=g(m,c,s,1990,'All ages','Number'); n23=g(m,c,s,2023,'All ages','Number')
      a90=g(m,c,s,1990,'Age-standardized','Rate'); a23=g(m,c,s,2023,'Age-standardized','Rate')
      dec=2 if ml!='DALYs' else 1
      rows.append({'Measure':ml,'Sex':s,'Etiology':cl,
        '1990 number (95% UI)':fn(n90),'1990 ASR per 100,000 (95% UI)':fn(a90,dec),
        '2023 number (95% UI)':fn(n23),'2023 ASR per 100,000 (95% UI)':fn(a23,dec),
        'Share of 2023 total (same sex and measure), %':f'{100*n23[0]/tot23:.1f}',
        'Change in number 1990–2023, %':f'{100*(n23[0]/n90[0]-1):+.1f}',
        'Change in ASR 1990–2023, %':f'{100*(a23[0]/a90[0]-1):+.1f}'})
      raw.append(dict(Measure=ml,Sex=s,Etiology=cl,n1990=n90[0],n1990_lo=n90[1],n1990_hi=n90[2],n2023=n23[0],n2023_lo=n23[1],n2023_hi=n23[2],
        asr1990=a90[0],asr1990_lo=a90[1],asr1990_hi=a90[2],asr2023=a23[0],asr2023_lo=a23[1],asr2023_hi=a23[2]))
t=pd.DataFrame(rows); r=pd.DataFrame(raw)
out='out/main/GBD2023_Table1_China_burden_by_sex_etiology.xlsx'
with pd.ExcelWriter(out) as w:
  for ml in meas.values(): t[t.Measure==ml].drop(columns='Measure').to_excel(w,sheet_name=f'Table1_{ml}',index=False)
  r.to_excel(w,sheet_name='raw_values',index=False)
  pd.DataFrame({'Note':['Source: GBD 2023 Results Tool, China, 1990 and 2023.','ASR = age-standardised rate per 100,000 (GBD world standard population).','95% UI = 95% uncertainty interval from GBD (25th and 975th of 1,000 draws).','Percentage change computed from point estimates; UIs not propagated.','"Share" = etiology number / all liver cancer number (incident cases, deaths or DALYs) in 2023 for the same sex and measure; etiologies sum to 100%.','MASH = GBD cause "Liver cancer due to NASH".']}).to_excel(w,sheet_name='Notes',index=False)
from openpyxl import load_workbook
wb=load_workbook(out)
for ws in wb.worksheets:
  for col in ws.columns: ws.column_dimensions[col[0].column_letter].width=min(34,max(len(str(c.value or '')) for c in col)+2)
wb.save(out)
t.to_csv('out/main/table1_china.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',20)
print(t[t.Sex!='Both'][['Measure','Sex','Etiology','2023 number (95% UI)','2023 ASR per 100,000 (95% UI)','Share of 2023 total (same sex and measure), %','Change in number 1990–2023, %','Change in ASR 1990–2023, %']].to_string())
