"""2026-11-19 review 8: new main Figure 4 - age-standardised HCV/HBV rate ratio by sex, China, 1990-2023 (GBD 2023).
Bands: approximate 95% intervals from marginal GBD UIs on log scale, assuming independence between etiologies."""
import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator
B=''
INK='#1d1d1f'; F='#b3707a'; M='#1d1d1f'; BG='#e8e8ea'
plt.rcParams.update({'figure.facecolor':BG,'axes.facecolor':BG,'savefig.facecolor':BG,'axes.edgecolor':'#1d1d1f80','xtick.color':'#1d1d1f99','ytick.color':'#1d1d1f99','xtick.major.size':0,'ytick.major.size':0,'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans'})
d=pd.read_csv(B+'out/main/master_dataset_long.csv')
d=d[(d.location=='China')&(d.age=='Age-standardized')&(d.metric=='Rate')]
R=pd.read_csv(B+'out/contrast/sex_contrasts.csv')
fig,axes=plt.subplots(1,2,figsize=(9.2,4.4),sharey=True)
for ax,m,lab in zip(axes,['Incidence','Deaths'],['A  Incidence','B  Deaths']):
    for s,c in [('Female',F),('Male',M)]:
        g={e:d[(d.measure==m)&(d.sex==s)&(d.etiology==e)].set_index('year').sort_index() for e in ['HCV','HBV']}
        r=g['HCV'].val/g['HBV'].val
        se=np.sqrt(sum(((np.log(g[e].upper)-np.log(g[e].lower))/3.92)**2 for e in g))
        ax.fill_between(r.index,r*np.exp(-1.96*se),r*np.exp(1.96*se),color=c,alpha=.15,lw=0)
        ax.plot(r.index,r,color=c,lw=2.2,label='Women' if s=='Female' else 'Men')
        ax.text(2023.6,r.iloc[-1],f'{r.iloc[-1]:.2f}',color=c,fontsize=9,va='center',weight='bold')
        ax.text(1989.6,r.iloc[0],f'{r.iloc[0]:.2f}',color=c,fontsize=9,va='center',ha='right')
    ax.axhline(1,color=INK,lw=.9,ls='--'); ax.set_yscale('log')
    ax.yaxis.set_major_locator(FixedLocator([0.1,0.2,0.5,1,2])); ax.yaxis.set_major_formatter(FixedFormatter(['0.1','0.2','0.5','1','2'])); ax.yaxis.set_minor_locator(NullLocator())
    ax.set_ylim(0.08,2.6); ax.set_xlim(1986,2028); ax.set_xticks([1990,2000,2010,2023])
    x=R[(R.measure==m)&(R.contrast=='(HCV − HBV) women − men')].iloc[0]
    ax.set_title(lab,loc='left',weight='bold',fontsize=11,color=INK)
    ax.text(0.03,0.04,f'Female − male slope difference ×100:\n{x.slope_diff_pct:.2f} per year (95% CI {x.slope_lo:.2f} to {x.slope_hi:.2f})',transform=ax.transAxes,fontsize=8.6,color=INK)
axes[0].set_ylabel('HCV/HBV ratio of age-standardised rates (log scale)')
axes[0].legend(frameon=False,loc='center left',bbox_to_anchor=(0.02,0.47))
fig.tight_layout()
fig.savefig(B+'manuscript/media_v2.5/media/media/fig_ratio.png',dpi=300)
fig.savefig(B+'out/main/Fig4_HCV_HBV_ratio_by_sex.png',dpi=300)
print('ok')
