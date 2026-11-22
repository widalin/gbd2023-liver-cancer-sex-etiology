"""Main figures/tables for the sex-divergent etiologic transition paper (plan v2).
Outputs to out/main/: Fig1-Fig4, FigS1 (2019-cutoff sensitivity), Table2 + crossover + sensitivity xlsx.
"""
import itertools, numpy as np, pandas as pd, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

BASE = ""
OUT = BASE + "out/main/"; os.makedirs(OUT, exist_ok=True)
a = pd.read_csv(BASE + "gbd/GBD2023_liver_cancer_etiology_Global_China_1990-2023.csv")
b = pd.read_csv(BASE + "gbd/raw/jkus_deaths_dalys_F_both/IHME-GBD_2023_DATA-6c99c8d9-1.csv")
import glob as _gl
_co = [pd.read_csv(_gl.glob(BASE + f"gbd/raw/{_k}/*.csv")[0]) for _k in
       ["coauthor_B1_jkus_male_inc_deaths_2023", "coauthor_B2_jkus_female_inc_2023", "coauthor_G5_jkus_both_inc_2023"]]
d = pd.concat([a, b] + _co)  # 2026-11-10: + co-author GBD 2023 JKUS male inc/deaths, female & both incidence
CAUSE = {"Liver cancer": "All", "Liver cancer due to hepatitis B": "HBV", "Liver cancer due to hepatitis C": "HCV",
         "Liver cancer due to alcohol use": "Alcohol", "Liver cancer due to NASH": "MASH",
         "Liver cancer due to other causes": "Other", "Hepatoblastoma": "Hepatoblastoma"}
MEAS = {"Deaths": "Deaths", "Incidence": "Incidence", "DALYs (Disability-Adjusted Life Years)": "DALYs"}
LOC = {"China": "China", "Japan": "Japan", "Republic of Korea": "Korea", "United States of America": "USA", "Global": "Global"}
d = d[d.cause_name.isin(CAUSE) & d.measure_name.isin(MEAS) & d.location_name.isin(LOC)].copy()
d["cause"] = d.cause_name.map(CAUSE); d["measure"] = d.measure_name.map(MEAS); d["loc"] = d.location_name.map(LOC)
d = d.drop_duplicates(["loc", "cause", "measure", "sex_name", "age_name", "metric_name", "year"])
ETI = ["HBV", "HCV", "Alcohol", "MASH", "Other", "Hepatoblastoma"]
COL = {"HBV": "#b1bcd2", "HCV": "#607868", "Alcohol": "#dbb691", "MASH": "#a79fb8", "Other": "#d9c489",
       "Hepatoblastoma": "#7e9c9b", "All": "#1d1d1f"}
SEXC = {"Male": "#1d1d1f", "Female": "#c9a3a8", "Both": "#7e9c9b"}
INK = "#1d1d1f"; BG = "#e8e8ea"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#1d1d1f80",
                     "axes.labelcolor": INK, "xtick.color": "#1d1d1f99", "ytick.color": "#1d1d1f99",
                     "xtick.major.size": 0, "ytick.major.size": 0, "axes.titleweight": "bold", "axes.titlesize": 11})

def get(loc, cause, meas, sex, age="Age-standardized", metric="Rate"):
    g = d[(d["loc"] == loc) & (d.cause == cause) & (d.measure == meas) & (d.sex_name == sex)
          & (d.age_name == age) & (d.metric_name == metric)].sort_values("year")
    return g

# ---------- joinpoint (same as joinpoint_jkus.py) ----------
def design(x, ks):
    return np.column_stack([np.ones_like(x), x] + [np.clip(x - k, 0, None) for k in ks])

def jp(x, y, maxk=3, gap=3):
    n = len(x); best = None
    cand = list(range(int(x[0]) + gap, int(x[-1]) - gap + 1))
    for nk in range(maxk + 1):
        for ks in itertools.combinations(cand, nk):
            if any(ks[i + 1] - ks[i] < gap for i in range(len(ks) - 1)): continue
            X = design(x, ks); beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            rss = float(((y - X @ beta) ** 2).sum()); p = X.shape[1] + nk
            bic = n * np.log(rss / n) + p * np.log(n)
            if best is None or bic < best[0]: best = (bic, ks, beta, rss, X)
    bic, ks, beta, rss, X = best
    p = X.shape[1]; df = n - p - len(ks); s2 = rss / df
    cov = s2 * np.linalg.inv(X.T @ X); tq = stats.t.ppf(0.975, df)
    bounds = [x[0], *ks, x[-1]]; segs = []
    f = lambda v: 100 * (np.exp(v) - 1)
    for i in range(len(bounds) - 1):
        c = np.zeros(p); c[1:2 + i] = 1
        sl = c @ beta; se = np.sqrt(c @ cov @ c)
        segs.append((int(bounds[i]), int(bounds[i + 1]), f(sl), f(sl - tq * se), f(sl + tq * se), c))
    w = np.array([s[1] - s[0] for s in segs], float); w /= w.sum()
    cA = sum(wi * s[5] for wi, s in zip(w, segs)); sl = cA @ beta; se = np.sqrt(cA @ cov @ cA)
    return ks, (f(sl), f(sl - tq * se), f(sl + tq * se)), [s[:5] for s in segs]

# ---- v2 (2026-11-06): NCI-like joinpoint with SE weights, cached; overrides jp() above ----
import sys, pickle; sys.path.insert(0, BASE + "gbd"); import jp2
_CF = BASE + "gbd/jp2_cache.pkl"
_cache = pickle.load(open(_CF, "rb")) if os.path.exists(_CF) else {}
def SE(g): return ((np.log(g.upper) - np.log(g.lower)) / 3.92).values
def jp(x, y, se=None, full=False):
    key = (tuple(np.round(x, 3)), tuple(np.round(y, 9)), se is None)
    if key not in _cache:
        _cache[key] = jp2.jp(x, y, se); pickle.dump(_cache, open(_CF, "wb"))
    r = _cache[key]
    return r if full else r[:3]

# ---- v3 (2026-11-08): official NCI Joinpoint 6.1.0 results for all full-period (1990-2023) series ----
import nci
_FP = {}
for (_L, _c, _m, _s), _g in d[(d.age_name == "Age-standardized") & (d.metric_name == "Rate")].groupby(["loc", "cause", "measure", "sex_name"]):
    _FP[tuple(np.round(np.log(_g.sort_values("year").val.values), 9))] = (_L, _m, _s, _c)
_jp_py = jp
def jp(x, y, se=None, full=False):
    k = _FP.get(tuple(np.round(y, 9)))
    if k is None or len(x) != 34: return _jp_py(x, y, se, full)
    try: ks, A, segs, fit, p = nci.result(*k, years=x)
    except IndexError: return _jp_py(x, y, se, full)   # series not yet run in NCI -> in-house jp2 (flag in outputs)
    return (ks, A, segs, fit) if full else (ks, A, segs)

def eapc(x, y):
    r = stats.linregress(x, np.log(y)); tq = stats.t.ppf(0.975, len(x) - 2)
    f = lambda v: 100 * (np.exp(v) - 1)
    return f(r.slope), f(r.slope - tq * r.stderr), f(r.slope + tq * r.stderr)

from decimal import Decimal, ROUND_HALF_UP
_hr = lambda v: float(Decimal(repr(float(v))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
fmt = lambda t: f"{_hr(t[0]):.2f} ({_hr(t[1]):.2f} to {_hr(t[2]):.2f})"  # half-up rounding (e.g. -0.415 -> -0.42)

# ---------- Table 2: China EAPC + AAPC by sex x etiology x measure ----------
t2, segrows = [], []
for meas in ["Incidence", "Deaths", "DALYs"]:
    for sex in ["Male", "Female", "Both"]:
        for c in ["All"] + ETI:
            g = get("China", c, meas, sex)
            x = g.year.values.astype(float); y = g.val.values
            e = eapc(x, y); ks, A, segs = jp(x, np.log(y), SE(g))
            t2.append(dict(measure=meas, sex=sex, etiology=c, ASR_1990=round(y[0], 2), ASR_2023=round(y[-1], 2),
                           pct_change=round(100 * (y[-1] / y[0] - 1), 1), EAPC=round(e[0], 2), EAPC_lo=round(e[1], 2),
                           EAPC_hi=round(e[2], 2), AAPC=round(A[0], 2), AAPC_lo=round(A[1], 2), AAPC_hi=round(A[2], 2),
                           EAPC_95CI=fmt(e), AAPC_95CI=fmt(A), joinpoints=", ".join(map(str, ks)) or "none"))
            for s in segs:
                segrows.append(dict(measure=meas, sex=sex, etiology=c, start=s[0], end=s[1], APC=round(s[2], 2),
                                    APC_lo=round(s[3], 2), APC_hi=round(s[4], 2)))
t2 = pd.DataFrame(t2); segdf = pd.DataFrame(segrows)

# ---------- Sensitivity: 1990-2019 vs 1990-2023 ----------
sens = []
for meas in ["Incidence", "Deaths"]:
    for sex in ["Male", "Female", "Both"]:
        for c in ["All", "HBV", "HCV", "Alcohol", "MASH", "Other"]:
            g = get("China", c, meas, sex); x = g.year.values.astype(float); y = g.val.values
            m = x <= 2019
            e1, e2 = eapc(x[m], y[m]), eapc(x, y)
            ks, A, _ = jp(x[m], np.log(y[m]), SE(g)[m])
            v19 = y[x == 2019][0]
            sens.append(dict(measure=meas, sex=sex, etiology=c, EAPC_1990_2019=fmt(e1), EAPC_1990_2023=fmt(e2),
                             AAPC_1990_2019_joinpoint=fmt(A), e19=e1[0], e19lo=e1[1], e19hi=e1[2], e23=e2[0],
                             e23lo=e2[1], e23hi=e2[2], ASR_2019=round(v19, 2), ASR_2023=round(y[-1], 2),
                             change_2019_2023_pct=round(100 * (y[-1] / v19 - 1), 1)))
sens = pd.DataFrame(sens)

# ---------- Crossover: HCV vs HBV (China) ----------
cross = []
for meas in ["Incidence", "Deaths", "DALYs"]:
    for sex in ["Female", "Male", "Both"]:
        for age, metric, lab in [("Age-standardized", "Rate", "ASR"), ("All ages", "Number", "Count")]:
            h = get("China", "HCV", meas, sex, age, metric).set_index("year").val
            bb = get("China", "HBV", meas, sex, age, metric).set_index("year").val
            diff = (h - bb)
            yrs = diff[diff > 0].index
            # first year after which HCV stays above HBV
            first = None
            for yv in diff.index:
                if (diff.loc[yv:] > 0).all(): first = int(yv); break
            cross.append(dict(measure=meas, sex=sex, metric=lab, HCV_exceeds_HBV_from=first if first else "not crossed",
                              HCV_2023=round(h.loc[2023], 2), HBV_2023=round(bb.loc[2023], 2),
                              HCV_share_2023_pct=round(100 * h.loc[2023] / get("China", "All", meas, sex, age, metric)
                                                       .set_index("year").val.loc[2023], 1) if lab == "Count" else ""))
cross = pd.DataFrame(cross)

with pd.ExcelWriter(OUT + "GBD2023_China_sex_etiology_tables.xlsx") as w:
    t2.drop(columns=["EAPC_lo", "EAPC_hi", "AAPC_lo", "AAPC_hi"]).to_excel(w, sheet_name="Table2_EAPC_AAPC", index=False)
    segdf.to_excel(w, sheet_name="Joinpoint_segments", index=False)
    cross.to_excel(w, sheet_name="HCV_HBV_crossover", index=False)
    sens[[c for c in sens.columns if not c.startswith(("e19", "e23"))]].to_excel(w, sheet_name="S1_2019_cutoff", index=False)
    pd.DataFrame({"Notes": [
        "Source: GBD 2023 Results Tool (IHME), China, 1990-2023, age-standardized rates per 100,000.",
        "EAPC = 100*(exp(b)-1), ln(ASR) ~ year OLS, 95% CI t-distribution. Point estimates only; GBD UI not propagated.",
        "AAPC (v2): Python joinpoint mirroring NCI Joinpoint defaults for annual data - log-linear, up to 5 joinpoints at integer years, >=2 obs from ends and between joinpoints, heteroscedastic weights 1/SE^2 with SE of ln(ASR) = (ln upper - ln lower)/3.92 from GBD 95% UI, BIC selection. SUPERSEDED 2026-11-08: Table 2 AAPC and segments now from official NCI Joinpoint 6.1.0 (log-linear, 0-5 joinpoints, WBIC, SE = (upper-lower)/3.92 heteroscedastic weights, Empirical Quantile CI 5001). S1 1990-2019 AAPC remains in-house.",
        "Crossover = first year from which HCV stays above HBV through 2023.",
        "S1: sensitivity analysis restricting to 1990-2019 (pre-COVID).",
        "MASH = GBD cause 'Liver cancer due to NASH'."]}).to_excel(w, sheet_name="Notes", index=False)
t2.to_csv(OUT + "table2_china.csv", index=False)

# ---------- Fig 1: etiology share of counts, China, by sex ----------
fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True, sharey=True)
for i, meas in enumerate(["Incidence", "Deaths"]):
    for j, sex in enumerate(["Male", "Female"]):
        ax = axes[i, j]
        tot = get("China", "All", meas, sex, "All ages", "Number").set_index("year").val
        st = [get("China", c, meas, sex, "All ages", "Number").set_index("year").val / tot * 100 for c in ETI]
        ax.stackplot(tot.index, st, colors=[COL[c] for c in ETI], labels=ETI, linewidth=0)
        ax.set_title(f"{meas}, {sex.lower()}s", loc="left"); ax.set_ylim(0, 100); ax.set_xlim(1990, 2023)
        ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
        for yr in [1990, 2023]:
            hb = (get("China", "HBV", meas, sex, "All ages", "Number").set_index("year").val / tot * 100).loc[yr]
            ax.text(yr + (0.4 if yr == 1990 else -0.4), hb / 2, f"HBV {hb:.0f}%", ha="left" if yr == 1990 else "right",
                    va="center", fontsize=8.5, color=INK)
fig.legend(*axes[0, 0].get_legend_handles_labels(), ncol=6, loc="upper left", bbox_to_anchor=(0.055, 0.945), frameon=False, fontsize=9)
fig.suptitle("HBV dominates male liver cancer in China, while HCV now rivals HBV in women", x=0.06, ha="left",
             fontsize=13, fontweight="bold", y=1.02)
fig.text(0.06, 0.975, "Share of liver cancer cases and deaths by etiology, China, 1990-2023 (GBD 2023)", fontsize=10,
         color="#1d1d1f99")
fig.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.06, wspace=0.08, hspace=0.25); fig.savefig(OUT + "Fig1_China_etiology_share_by_sex.png", dpi=200, bbox_inches="tight"); plt.close()

# ---------- Fig 2: ASR trends by sex, crossover marked ----------
fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
for i, meas in enumerate(["Incidence", "Deaths"]):
    for j, sex in enumerate(["Male", "Female"]):
        ax = axes[i, j]
        for c in ["HBV", "HCV", "Alcohol", "MASH", "Other"]:
            g = get("China", c, meas, sex)
            ax.plot(g.year, g.val, color=COL[c], lw=2 if c in ("HBV", "HCV", "MASH") else 1.5, label=c,
                    solid_capstyle="round")
        ax.grid(axis="y", color="#1d1d1f12"); ax.set_xlim(1990, 2023)
        if sex == "Male":
            ax.set_yscale("log"); ax.set_ylim(0.3, 15); ax.set_yticks([0.5, 1, 2, 5, 10]); ax.set_yticklabels(["0.5", "1", "2", "5", "10"]); ax.minorticks_off()
        else: ax.set_ylim(bottom=0)
        ax.set_title(f"{meas}, {sex.lower()}s (per 100,000{', log scale' if sex == 'Male' else ''})", loc="left")
        r = cross[(cross.measure == meas) & (cross.sex == sex) & (cross.metric == "ASR")].iloc[0]
        if r.HCV_exceeds_HBV_from != "not crossed":
            yr = int(r.HCV_exceeds_HBV_from); v = get("China", "HCV", meas, sex).set_index("year").val.loc[yr]
            ax.plot(yr, v, "o", color="#607868", ms=6); ax.annotate(f"HCV > HBV from {yr}", (yr, v), xytext=(6, 12),
                                                                     textcoords="offset points", fontsize=9, color=INK)
        if False:
            g = get("China", "MASH", meas, sex)
            e = eapc(g.year.values.astype(float), g.val.values)
            ax.annotate(f"MASH EAPC {e[0]:+.2f}%", (2023, g.val.values[-1]), xytext=(-2, -16), textcoords="offset points",
                        ha="right", fontsize=9, color="#6f6782")
fig.legend(*axes[0, 0].get_legend_handles_labels(), ncol=5, loc="upper left", bbox_to_anchor=(0.055, 0.945), frameon=False, fontsize=9)
fig.suptitle("Etiology-specific liver cancer rates diverge by sex in China", x=0.06, ha="left", fontsize=13,
             fontweight="bold", y=1.02)
fig.text(0.06, 0.975, "Age-standardized incidence and death rates per 100,000, 1990-2023 (GBD 2023)", fontsize=10,
         color="#1d1d1f99")
fig.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.06, wspace=0.15, hspace=0.25); fig.savefig(OUT + "Fig2_China_ASR_by_sex_crossover.png", dpi=200, bbox_inches="tight"); plt.close()

# ---------- Fig 3: AAPC forest by sex x etiology (China) ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 6), sharey=True)
order = ["All", "HBV", "HCV", "Alcohol", "MASH", "Other"]
for ax, meas in zip(axes, ["Incidence", "Deaths"]):
    for k, sex in enumerate(["Male", "Female", "Both"]):
        sub = t2[(t2.measure == meas) & (t2.sex == sex)].set_index("etiology").loc[order]
        yy = np.arange(len(order)) * 1.0 + (k - 1) * 0.25
        ax.errorbar(sub.AAPC, yy, xerr=[sub.AAPC - sub.AAPC_lo, sub.AAPC_hi - sub.AAPC], fmt="o", color=SEXC[sex],
                    ecolor=SEXC[sex], ms=5, lw=1.4, capsize=0, label=sex)
    ax.axvline(0, color="#1d1d1f24", ls="--", lw=1)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(["All liver cancer" if o == "All" else o for o in order], fontsize=11)
    ax.invert_yaxis(); ax.grid(axis="x", color="#1d1d1f12")
    ax.set_title(f"{meas} ASR", loc="left"); ax.set_xlabel("AAPC 1990-2023, % per year (95% CI)")
fig.legend(*axes[0].get_legend_handles_labels(), ncol=3, frameon=False, loc="upper left", bbox_to_anchor=(0.06, 0.94), fontsize=9)
fig.suptitle("In Chinese men, MASH- and alcohol-related incidence showed no significant decline, unlike viral and other etiologies", x=0.06, ha="left", fontsize=13,
             fontweight="bold", y=1.03)
fig.text(0.06, 0.965, "Joinpoint average annual percent change of age-standardized rates, by sex (NCI Joinpoint 6.1.0)", fontsize=10,
         color="#1d1d1f99")
fig.subplots_adjust(left=0.1, right=0.98, top=0.82, bottom=0.1, wspace=0.08); fig.savefig(OUT + "Fig3_China_AAPC_by_sex.png", dpi=200, bbox_inches="tight"); plt.close()

# ---------- Fig 4: China vs Japan/Korea/USA mirrors ----------
locs = ["China", "Japan", "Korea", "USA"]
fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]; xt, xl = [], []
for li, L in enumerate(locs):
    for k, (sex, yr) in enumerate([("Male", 1990), ("Male", 2023), ("Female", 1990), ("Female", 2023)]):
        xpos = li * 5 + k + (0.3 if k >= 2 else 0)
        tot = get(L, "All", "Deaths", sex, "All ages", "Number").set_index("year").val.loc[yr]
        bottom = 0
        for c in ETI[:5]:
            v = get(L, c, "Deaths", sex, "All ages", "Number").set_index("year").val.loc[yr] / tot * 100
            ax.bar(xpos, v, bottom=bottom, color=COL[c], width=0.85, label=c if (li == 0 and k == 0) else None)
            bottom += v
        xt.append(xpos); xl.append(f"{sex[0]}\n{str(yr)[2:]}")
    ax.text(li * 5 + 1.65, 104, L, ha="center", fontsize=11, fontweight="bold", color=INK)
ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=8); ax.set_ylim(0, 110)
ax.set_yticks([0, 25, 50, 75, 100]); ax.yaxis.set_major_formatter(lambda v, p: f"{v:.0f}%")
ax.set_title("Share of liver cancer deaths by etiology", loc="left", pad=18)
ax.legend(ncol=5, frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0, -0.13))
ax.text(0, -0.24, "M = males, F = females; 90 = 1990, 23 = 2023; hepatoblastoma (<1.5%) omitted", transform=ax.transAxes, fontsize=8,
        color="#1d1d1f99", va="top", ha="left")
ax = axes[1]
_rows = []
for _L in ["China", "Japan", "Korea", "USA"]:
    for _s in ["Male", "Female"]:
        for _c in ["All", "HBV", "HCV", "Alcohol", "MASH"]:
            try: _r = nci.result(_L, "Deaths", _s, _c)[1]; _src = "NCI Joinpoint 6.1.0"
            except IndexError:
                _g = get(_L, _c, "Deaths", _s); _r = jp(_g.year.values.astype(float), np.log(_g.val.values), SE(_g))[1]
                _src = "in-house jp2 (NCI rerun pending)"
            _rows.append(dict(location=_L, sex=_s, cause="All liver cancer" if _c == "All" else _c, AAPC=_r[0], AAPC_lo=_r[1], AAPC_hi=_r[2], source=_src))
jk = pd.DataFrame(_rows); jk.to_csv(OUT + "fig4_JKUS_deaths_AAPC_by_sex.csv", index=False)
LC = {"China": "#1d1d1f", "Japan": "#b1bcd2", "Korea": "#dbb691", "USA": "#a79fb8"}
causes = ["All liver cancer", "HBV", "HCV", "Alcohol", "MASH"]
for k, L in enumerate(locs):
    for _s, _mk in [("Male", "o"), ("Female", "D")]:
        sub = jk[(jk.location == L) & (jk.sex == _s)].set_index("cause").reindex(causes)
        yy = np.arange(len(causes)) * 2 + (0 if _s == "Male" else 0.9) + (k - 1.5) * 0.17
        ax.errorbar(sub.AAPC, yy, xerr=[sub.AAPC - sub.AAPC_lo, sub.AAPC_hi - sub.AAPC], fmt=_mk, color=LC[L], mfc=LC[L] if _s == "Male" else "none",
                    ecolor=LC[L], ms=4.5, lw=1.2, capsize=0, label=L if _s == "Male" else None)
ax.axvline(0, color="#1d1d1f24", ls="--", lw=1)
ax.set_yticks([i * 2 + d_ for i in range(len(causes)) for d_ in (0, 0.9)])
ax.set_yticklabels([f"{c_}  M" if d_ == 0 else "F" for c_ in ["All liver cancer"] + causes[1:] for d_ in (0, 0.9)], fontsize=8.5)
ax.invert_yaxis(); ax.grid(axis="x", color="#1d1d1f12")
ax.set_title("Death ASR, AAPC 1990-2023, by sex", loc="left", pad=18); ax.set_xlabel("% per year (95% CI); filled = male, open = female")
ax.legend(frameon=False, fontsize=9, loc="lower left")
fig.suptitle("China, Japan, Korea and the USA follow distinct etiologic paths of liver cancer deaths",
             x=0.04, ha="left", fontsize=12.5, fontweight="bold", y=1.04)
fig.text(0.04, 0.975, "Liver cancer deaths by sex, China vs Japan, Korea and USA (GBD 2023)", fontsize=10,
         color="#1d1d1f99")
fig.tight_layout(); fig.savefig(OUT + "Fig4_China_vs_JKUS_mirror.png", dpi=200, bbox_inches="tight"); plt.close()

# ---------- Fig S1: 2019 cutoff sensitivity ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 6), sharey=True)
order = ["All", "HBV", "HCV", "Alcohol", "MASH", "Other"]
for ax, meas in zip(axes, ["Incidence", "Deaths"]):
    for k, sex in enumerate(["Male", "Female"]):
        sub = sens[(sens.measure == meas) & (sens.sex == sex)].set_index("etiology").loc[order]
        base = np.arange(len(order)) + (k - 0.5) * 0.36
        ax.errorbar(sub.e23, base - 0.07, xerr=[sub.e23 - sub.e23lo, sub.e23hi - sub.e23], fmt="o", color=SEXC[sex],
                    ms=5, lw=1.3, capsize=0, label=f"{sex}, 1990-2023")
        ax.errorbar(sub.e19, base + 0.07, xerr=[sub.e19 - sub.e19lo, sub.e19hi - sub.e19], fmt="o", mfc="none",
                    color=SEXC[sex], ms=5, lw=1.0, capsize=0, label=f"{sex}, 1990-2019")
    ax.axvline(0, color="#1d1d1f24", ls="--", lw=1); ax.set_yticks(range(len(order)))
    ax.set_yticklabels(["All liver cancer" if o == "All" else o for o in order], fontsize=11); ax.invert_yaxis()
    ax.grid(axis="x", color="#1d1d1f12"); ax.set_title(f"{meas} ASR", loc="left"); ax.set_xlabel("EAPC, % per year (95% CI)")
fig.legend(*axes[0].get_legend_handles_labels(), ncol=4, frameon=False, loc="upper left", bbox_to_anchor=(0.06, 0.94), fontsize=9)
fig.suptitle("Trend estimates barely change when COVID-era years are excluded", x=0.06, ha="left", fontsize=13,
             fontweight="bold", y=1.03)
fig.text(0.06, 0.965, "Sensitivity analysis: EAPC 1990-2019 vs 1990-2023, China, by sex", fontsize=10, color="#1d1d1f99")
fig.subplots_adjust(left=0.1, right=0.98, top=0.82, bottom=0.1, wspace=0.08); fig.savefig(OUT + "FigS1_2019_cutoff_sensitivity.png", dpi=200, bbox_inches="tight"); plt.close()

pd.set_option("display.width", 250)
print(t2[["measure", "sex", "etiology", "pct_change", "EAPC_95CI", "AAPC_95CI", "joinpoints"]].to_string())
print(cross.to_string()); print(sens[["measure", "sex", "etiology", "EAPC_1990_2019", "EAPC_1990_2023", "change_2019_2023_pct"]].to_string())
