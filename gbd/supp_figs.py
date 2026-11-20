import matplotlib.ticker
"""Supplementary figures S2-S5 + integrated master dataset. Reuses header of main_figs.py."""
src = open("gbd/main_figs.py").read()
exec(src.split("# ---------- Table 2")[0])
E5 = ["HBV", "HCV", "Alcohol", "MASH", "Other"]
t2 = pd.read_csv(OUT + "table2_china.csv")

def head(fig, title, sub, x=0.06):
    fig.suptitle(title, x=x, ha="left", fontsize=13, fontweight="bold", y=1.03)
    fig.text(x, 0.965, sub, fontsize=10, color="#1d1d1f99")

# ---- S2: China vs Global by sex ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]; pos = 0; ticks = []; labs = []
for loc in ["China", "Global"]:
    for sex in ["Male", "Female"]:
        for yr in [1990, 2023]:
            tot = get(loc, "All", "Incidence", sex, "All ages", "Number").set_index("year").val[yr]
            bot = 0
            for c in ETI:
                v = 100 * get(loc, c, "Incidence", sex, "All ages", "Number").set_index("year").val[yr] / tot
                ax.bar(pos, v, bottom=bot, color=COL[c], width=0.8, label=c if pos == 0 else None); bot += v
            ticks.append(pos); labs.append(f"{sex[0]}\n{str(yr)[2:]}"); pos += 1
        pos += 0.4
    ax.text(pos - 2.9 - (0.0), 104, loc, ha="center", fontweight="bold", fontsize=11)
    pos += 0.8
ax.set_xticks(ticks); ax.set_xticklabels(labs, fontsize=8.5); ax.set_ylim(0, 110)
ax.set_yticks([0, 25, 50, 75, 100]); ax.set_yticklabels([f"{v}%" for v in [0, 25, 50, 75, 100]])
ax.set_title("Share of incident cases by etiology", loc="left", pad=18)
ax.legend(ncol=6, frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0, -0.12))
ax = axes[1]; rows = []
for loc in ["China", "Global"]:
    for sex in ["Male", "Female"]:
        for c in ["All"] + E5:
            g = get(loc, c, "Incidence", sex); x = g.year.values.astype(float); y = g.val.values
            _, A, _ = jp(x, np.log(y), SE(g)); rows.append((loc, sex, c, *A))
S2 = pd.DataFrame(rows, columns=["loc", "sex", "etiology", "AAPC", "lo", "hi"])
style = {("China", "Male"): ("#1d1d1f", "o"), ("China", "Female"): ("#c9a3a8", "o"),
         ("Global", "Male"): ("#1d1d1f", "s"), ("Global", "Female"): ("#c9a3a8", "s")}
off = {("China", "Male"): -0.27, ("China", "Female"): -0.09, ("Global", "Male"): 0.09, ("Global", "Female"): 0.27}
for (loc, sex), g in S2.groupby(["loc", "sex"]):
    cc, mk = style[(loc, sex)]
    for i, c in enumerate(["All"] + E5):
        r = g[g.etiology == c].iloc[0]
        ax.errorbar(r.AAPC, i + off[(loc, sex)], xerr=[[r.AAPC - r.lo], [r.hi - r.AAPC]], fmt=mk, color=cc,
                    mfc=cc if loc == "China" else "none", ms=6, lw=1.6, label=f"{loc}, {sex.lower()}" if i == 0 else None)
ax.set_yticks(range(6)); ax.set_yticklabels(["All causes"] + E5); ax.invert_yaxis()
ax.axvline(0, color="#1d1d1f55", ls="--", lw=1); ax.grid(axis="x", color="#1d1d1f12")
ax.set_xlabel("AAPC 1990-2023, % per year (95% CI)"); ax.set_title("Incidence ASR trend", loc="left", pad=18)
ax.legend(frameon=False, fontsize=8.5, loc="lower right"); ax.set_xlim(right=1.6)
head(fig, "Unlike the global pattern, liver cancer incidence is falling for every etiology in Chinese women",
     "China vs global, incidence by sex (GBD 2023). Filled = China, open = global. M = male, F = female; 90 = 1990, 23 = 2023", x=0.04)
fig.tight_layout(); fig.savefig(OUT + "FigS2_China_vs_Global_by_sex.png", dpi=200, bbox_inches="tight"); plt.close()

# ---- S3: male-to-female ASR ratio ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
S3 = []
for ax, meas in zip(axes, ["Incidence", "Deaths"]):
    for c in ["All"] + E5:
        m = get("China", c, meas, "Male").set_index("year").val; f = get("China", c, meas, "Female").set_index("year").val
        r = m / f; S3.append(pd.DataFrame({"measure": meas, "etiology": c, "year": r.index, "MF_ratio": r.values.round(3)}))
        ax.plot(r.index, r.values, color=COL[c], lw=2.6 if c in ("All", "MASH") else 1.8, label="All causes" if c == "All" else c)
        if c != "Other": ax.text(2023.4, r.values[-1] * (1.03 if c == "HCV" else 1), f"{r.values[-1]:.1f}" + (" HCV/other" if c == "HCV" else ""), fontsize=8.5, va="center", color=COL[c] if c != "HCV" else INK)
    ax.axhline(1, color="#1d1d1f55", ls="--", lw=1); ax.set_yscale("log")
    ax.set_yticks([0.5, 1, 2, 4, 8]); ax.set_yticklabels(["0.5", "1", "2", "4", "8"]); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlim(1990, 2028)
    ax.grid(axis="y", color="#1d1d1f12"); ax.set_title(f"{meas} ASR, male / female (log scale)", loc="left")
fig.legend(*axes[0].get_legend_handles_labels(), ncol=6, frameon=False, loc="upper left", bbox_to_anchor=(0.06, 0.93), fontsize=9)
head(fig, "Faster declines in women widened the male excess for every etiology; MASH flipped from a female to a male excess",
     "Male-to-female ratio of age-standardized rates by etiology, China, 1990-2023")
fig.subplots_adjust(left=0.06, right=0.97, top=0.8, bottom=0.08, wspace=0.08)
fig.savefig(OUT + "FigS3_China_male_female_ratio.png", dpi=200, bbox_inches="tight"); plt.close(); S3 = pd.concat(S3)

# ---- S4: mortality-to-incidence ratio ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
S4 = []
for ax, sex in zip(axes, ["Male", "Female"]):
    for c in ["All"] + E5:
        dd = get("China", c, "Deaths", sex).set_index("year").val; ii = get("China", c, "Incidence", sex).set_index("year").val
        r = dd / ii; S4.append(pd.DataFrame({"sex": sex, "etiology": c, "year": r.index, "MIR": r.values.round(3)}))
        ax.plot(r.index, r.values, color=COL[c], lw=2.6 if c == "All" else 1.8, label="All causes" if c == "All" else c)
    ax.set_xlim(1990, 2023); ax.grid(axis="y", color="#1d1d1f12"); ax.set_title(f"{sex}s", loc="left")
axes[0].set_ylabel("Death ASR / incidence ASR")
fig.legend(*axes[0].get_legend_handles_labels(), ncol=6, frameon=False, loc="upper left", bbox_to_anchor=(0.06, 0.93), fontsize=9)
head(fig, "Mortality-to-incidence ratios have fallen for every etiology but remain high",
     "Age-standardized mortality-to-incidence ratio (MIR) by etiology and sex, China, 1990-2023; a crude proxy of case fatality")
fig.subplots_adjust(left=0.07, right=0.98, top=0.8, bottom=0.08, wspace=0.06)
fig.savefig(OUT + "excluded_MIR_by_sex.png", dpi=200, bbox_inches="tight"); plt.close(); S4 = pd.concat(S4)

# ---- S5: joinpoint segment fits ----
for meas in ["Incidence", "Deaths"]:
    fig, axes = plt.subplots(2, 5, figsize=(16, 6.4), sharex=True)
    for i, sex in enumerate(["Male", "Female"]):
        for j, c in enumerate(E5):
            ax = axes[i, j]; g = get("China", c, meas, sex); x = g.year.values.astype(float); y = g.val.values
            ks, A, segs, beta = jp(x, np.log(y), SE(g), full=True)
            fitv = beta if np.shape(beta) == np.shape(x) else np.exp(design(x, ks) @ beta)
            ax.scatter(x, y, s=10, color=COL[c], zorder=3); ax.plot(x, fitv, color=INK, lw=1.4)
            for k in ks: ax.axvline(k, color="#1d1d1f40", ls=":", lw=1)
            txt = "\n".join(f"{s[0]}-{s[1]}: {s[2]:+.1f}" for s in segs)
            ax.text(0.02, 0.03, txt, transform=ax.transAxes, fontsize=7, va="bottom", color="#1d1d1fcc")
            ax.set_title(f"{c}, {sex.lower()}s", loc="left", fontsize=10); ax.grid(axis="y", color="#1d1d1f12")
            ax.text(0.98, 0.95, f"AAPC {A[0]:+.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=8, fontweight="bold")
    n = "S4" if meas == "Incidence" else "S5"
    head(fig, f"Joinpoint segments of liver cancer {meas.lower()} rates by etiology and sex, China",
         f"Age-standardized {meas.lower()} rate per 100,000 (points) with fitted NCI Joinpoint 6.1.0 model; segment APC (% per year) at bottom left", x=0.04)
    fig.tight_layout(); fig.savefig(OUT + f"Fig{n}_China_joinpoint_{meas.lower()}_by_sex.png", dpi=200, bbox_inches="tight"); plt.close()

# ---- integrated master dataset ----
m = d[["loc", "measure", "cause", "sex_name", "age_name", "metric_name", "year", "val", "lower", "upper"]].rename(
    columns={"loc": "location", "cause": "etiology", "sex_name": "sex", "age_name": "age", "metric_name": "metric"})
m = m.sort_values(["location", "measure", "etiology", "sex", "age", "metric", "year"])
m.to_csv(OUT + "master_dataset_long.csv", index=False)
inv = m.groupby(["location", "measure", "sex"]).agg(etiologies=("etiology", lambda s: ", ".join(sorted(s.unique()))),
      years=("year", lambda s: f"{s.min()}-{s.max()}"), rows=("val", "size")).reset_index()
missing = pd.DataFrame({"Dataset still needed": [
    "Japan, Korea, USA: males (all measures) and incidence (all sexes)",
    "China (and JKUS): 5-year age groups, 1990-2023, by sex and etiology - for age-period-cohort (core) and decomposition",
    "Population by age and sex, China 1990-2023 (+ IHME forecast) - for decomposition and BAPC",
    "Risk-factor attributable burden (high BMI, high FPG, alcohol, smoking, drug use) for all years",
    "SDI values (GHDx)"], "Used for": ["Fig4 complete", "Fig5, Fig6", "Fig6, Fig7", "Supplementary", "Optional"]})
with pd.ExcelWriter(OUT + "GBD2023_integrated_data.xlsx") as w:
    pd.read_excel(OUT + "GBD2023_Table1_China_burden_by_sex_etiology.xlsx", sheet_name=None)
    for sh, df in pd.read_excel(OUT + "GBD2023_Table1_China_burden_by_sex_etiology.xlsx", sheet_name=None).items():
        if sh.startswith("Table1"): df.to_excel(w, sheet_name=sh, index=False)
    for sh, df in pd.read_excel(OUT + "GBD2023_China_sex_etiology_tables.xlsx", sheet_name=None).items():
        if sh != "Notes": df.to_excel(w, sheet_name=sh[:31], index=False)
    S2.round(2).to_excel(w, sheet_name="S2_China_vs_Global_AAPC", index=False)
    S3.to_excel(w, sheet_name="S3_MF_ratio", index=False)
    S4.to_excel(w, sheet_name="S4_MIR", index=False)
    inv.to_excel(w, sheet_name="Data_inventory", index=False)
    missing.to_excel(w, sheet_name="Data_still_needed", index=False)
from openpyxl import load_workbook
wb = load_workbook(OUT + "GBD2023_integrated_data.xlsx")
for ws in wb.worksheets:
    for col in ws.columns: ws.column_dimensions[col[0].column_letter].width = min(40, max(len(str(c.value or "")) for c in col[:60]) + 2)
wb.save(OUT + "GBD2023_integrated_data.xlsx")
print(S2.round(2).to_string()); print(S3[S3.year.isin([1990, 2023])].to_string()); print(S4[S4.year.isin([1990, 2023])].to_string())
print(len(m), "master rows")
