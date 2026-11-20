"""Fig 6 v2.5 (2026-11-11, reviewer Fig 6 note): top = all etiologies (thousands); bottom = MASH and alcohol
rows at case scale with each component labelled and its share of the net change, so the small positive rate
contribution for male MASH (+231, 6.3%) is readable without relying on colour. Data: out/apc/decomposition_gbd2023.csv"""
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = ""
R = pd.read_csv(B + "out/apc/decomposition_gbd2023.csv")
F = ["population", "ageing", "epidemiological"]
LAB = {"population": "Population growth", "ageing": "Population ageing", "epidemiological": "Rate change"}
COL = {"population": "#7f93b8", "ageing": "#c9935e", "epidemiological": "#607868"}
fig = plt.figure(figsize=(14, 11)); gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1], hspace=0.38, wspace=0.32)
def bars(a, rows, scale, labels=False):
    yl = []
    for yi, (e, s, m) in enumerate(rows):
        r = R[(R.measure == m) & (R.sex == s) & (R.etiology == e)].iloc[0]; pos = neg = 0
        for f in F:
            v = r[f] / scale; left = pos if v >= 0 else neg
            a.barh(yi, v, left=left, color=COL[f], height=0.7, label=LAB[f] if yi == 0 else None)
            if labels:
                pct = f"\n({100*r[f]/r.change:.1f}%)" if (e, s) == ("MASH", "Male") and m == "Incidence" else ""
                if abs(v) >= 550:
                    a.text(left + v / 2, yi, f"{r[f]:+,.0f}{pct}", ha="center", va="center", fontsize=7.8, color="white")
                else:
                    a.annotate(f"{r[f]:+,.0f}{pct}", (left + v / 2, yi - 0.36), xytext=(0, 9 if not pct else 4), textcoords="offset points",
                               ha="center", va="bottom", fontsize=7.8, color="#1d1d1f", arrowprops=dict(arrowstyle="-", lw=0.6, color="#555"))
            if v >= 0: pos += v
            else: neg += v
        a.plot(r.change / scale, yi, "D", color="black", ms=5, label="Net change" if yi == 0 else None)
        if labels: a.text(max(pos, r.change / scale) + 120, yi, f"net {r.change:+,.0f}", va="center", fontsize=8.5, fontweight="bold")
        yl.append(f"{e}, {s.lower()}s")
    a.set_yticks(range(len(rows))); a.set_yticklabels(yl, fontsize=8.8); a.invert_yaxis(); a.axvline(0, color="grey", lw=0.8)
    a.spines[["top", "right"]].set_visible(False)
order = ["All", "HBV", "HCV", "Alcohol", "MASH", "Other"]
for k, m in enumerate(["Incidence", "Deaths"]):
    a = fig.add_subplot(gs[0, k]); bars(a, [(e, s, m) for e in order for s in ["Male", "Female"]], 1000)
    a.set_title(f"{'AB'[k]}  {m}, all etiologies (thousands)", loc="left", fontsize=10.5, fontweight="bold")
    if k == 0: fig.legend(*a.get_legend_handles_labels(), ncol=4, fontsize=9.5, frameon=False, loc="upper left", bbox_to_anchor=(0.01, 0.93))
    b = fig.add_subplot(gs[1, k]); bars(b, [(e, s, m) for e in ["MASH", "Alcohol"] for s in ["Male", "Female"]], 1, labels=True)
    b.set_xlim(-2000, 6400); b.set_ylim(3.6, -0.75); b.set_xlabel("Cases" if m == "Incidence" else "Deaths", fontsize=9)
    b.set_title(f"{'CD'[k]}  {m}, MASH and alcohol (number)", loc="left", fontsize=10.5, fontweight="bold")
fig.suptitle("Decomposition of the change in liver cancer cases and deaths by etiology and sex, China, 1990-2023 (GBD 2023, ages 15+)",
             x=0.01, ha="left", fontsize=12.5, fontweight="bold", y=0.965)
fig.text(0.01, 0.06, "Das Gupta three-factor decomposition. Components sum to the net change. For male MASH-related incidence, rate change added 231 cases (6.3% of the net increase of 3,651).",
         fontsize=8.5, color="#555")
fig.savefig(B + "out/apc/Fig6_decomposition_v2.5.png", dpi=200, bbox_inches="tight")
