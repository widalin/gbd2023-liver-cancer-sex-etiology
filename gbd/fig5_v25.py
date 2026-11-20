"""Fig 5 v2.5 (2026-11-11, reviewer Fig 5 note): main-text APC figure limited to the key comparisons
(men: HBV, alcohol, MASH; women: HBV, HCV), local drift / period RR / cohort RR with 95% CIs for every line.
Full 8-panel figure (all etiologies + longitudinal age curves) kept as Supplementary Figure S10.
Data: out/apc/apc_functions_nci2023.csv + apc_summary_nci2023.csv (NCI APC Web Tool, default settings)."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, shutil
B = ""
K = pd.read_csv(B + "out/apc/apc_functions_nci2023.csv"); S = pd.read_csv(B + "out/apc/apc_summary_nci2023.csv")
K = K[K.measure == "Incidence"]; S = S[S.measure == "Incidence"]
COL = {"HBV": "#7f93b8", "HCV": "#4f6a57", "Alcohol": "#c9935e", "MASH": "#8a7fa8"}
SETS = {"Male": ["HBV", "Alcohol", "MASH"], "Female": ["HBV", "HCV"]}
FN = [("local_drift", "Local drift (% per year)", "Age (years)"), ("period_RR", "Period rate ratio", "Period (midpoint)"),
      ("cohort_RR", "Cohort rate ratio", "Birth cohort (midpoint)")]
fig, ax = plt.subplots(2, 3, figsize=(14, 8.4))
for i, sex in enumerate(["Male", "Female"]):
    for j, (fn, title, xl) in enumerate(FN):
        a = ax[i, j]
        for e in SETS[sex]:
            g = K[(K.sex == sex) & (K.etiology == e) & (K.function == fn)].sort_values("x")
            if fn == "cohort_RR": g = g[g.x >= 1910]
            s = S[(S.sex == sex) & (S.etiology == e)].iloc[0]
            a.fill_between(g.x, g.lo, g.hi, color=COL[e], alpha=0.18, lw=0)
            a.plot(g.x, g.est, color=COL[e], lw=2.2, label=f"{e} (net drift {s.net_drift:.2f}% [{s.nd_lo:.2f} to {s.nd_hi:.2f}])")
        a.axhline(0 if fn == "local_drift" else 1, color="#555", lw=0.8, ls=":")
        if fn == "cohort_RR": a.set_yscale("log"); a.set_yticks([0.125, 0.25, 0.5, 1, 2]); a.set_yticklabels(["0.125", "0.25", "0.5", "1", "2"]); a.minorticks_off()
        a.set_title(f"{'ABCDEF'[3*i+j]}  {'Men' if sex == 'Male' else 'Women'}: {title}", loc="left", fontsize=11, fontweight="bold")
        a.set_xlabel(xl, fontsize=9); a.spines[["top", "right"]].set_visible(False)
        if j == 0: a.legend(fontsize=8.3, frameon=False, loc="lower right")
fig.suptitle("Age-period-cohort analysis of key etiologies: liver cancer incidence in China (GBD 2023, ages 20-94, 1994-2023)",
             x=0.01, ha="left", fontsize=12.5, fontweight="bold")
fig.text(0.01, 0.005, "Shaded bands: 95% CIs. Reference period 2004-2008; reference birth cohort 1949. NCI Age-Period-Cohort Web Tool (default settings). "
         "Cohorts before 1910 omitted (wide CIs); all etiologies in Supplementary Figure S10.", fontsize=8.3, color="#555")
fig.tight_layout(rect=(0, 0.025, 1, 0.95))
fig.savefig(B + "out/apc/Fig5_APC_key_etiologies_v2.5.png", dpi=200)
