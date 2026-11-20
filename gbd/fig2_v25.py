"""Fig 2 v2.5 (2026-11-11): adds 95% UI bands for female HBV/HCV (reviewer #4) and a bottom row of
incidence indexed to 1990=100 on a common linear scale for both sexes (reviewer Fig 2 note)."""
src = open("gbd/main_figs.py").read().split("# ---------- joinpoint")[0]
exec(src)
CROSS = {("Incidence", "Female"): 2012, ("Deaths", "Female"): 2006}
fig, axes = plt.subplots(3, 2, figsize=(11, 10.5), sharex=True)
for i, meas in enumerate(["Incidence", "Deaths"]):
    for j, sex in enumerate(["Male", "Female"]):
        ax = axes[i, j]
        for c in ["HBV", "HCV", "Alcohol", "MASH", "Other"]:
            g = get("China", c, meas, sex)
            if sex == "Female" and c in ("HBV", "HCV"):
                ax.fill_between(g.year, g.lower, g.upper, color=COL[c], alpha=0.28, lw=0)
            ax.plot(g.year, g.val, color=COL[c], lw=2 if c in ("HBV", "HCV", "MASH") else 1.5, label=c, solid_capstyle="round")
        ax.grid(axis="y", color="#1d1d1f12"); ax.set_xlim(1990, 2023)
        if sex == "Male":
            ax.set_yscale("log"); ax.set_ylim(0.3, 15); ax.set_yticks([0.5, 1, 2, 5, 10]); ax.set_yticklabels(["0.5", "1", "2", "5", "10"]); ax.minorticks_off()
        else: ax.set_ylim(bottom=0)
        ax.set_title(f"{'ABCD'[2*i+j]}  {meas}, {sex.lower()}s (per 100,000{', log scale' if sex == 'Male' else ''})", loc="left")
        if (meas, sex) in CROSS:
            yr = CROSS[(meas, sex)]; v = get("China", "HCV", meas, sex).set_index("year").val.loc[yr]
            ax.plot(yr, v, "o", color="#607868", ms=6)
            ax.annotate(f"HCV > HBV point estimate\nfrom {yr} (UIs overlap)", (yr, v), xytext=(8, 26), textcoords="offset points", fontsize=8.5, color=INK)
for j, sex in enumerate(["Male", "Female"]):
    ax = axes[2, j]
    for c in ["HBV", "HCV", "Alcohol", "MASH", "Other"]:
        g = get("China", "%s" % c, "Incidence", sex)
        ax.plot(g.year, 100 * g.val.values / g.val.values[0], color=COL[c], lw=2 if c in ("HBV", "HCV", "MASH") else 1.5)
    ax.axhline(100, color="#1d1d1f55", lw=0.8, ls="--")
    ax.set_ylim(35, 125); ax.grid(axis="y", color="#1d1d1f12"); ax.set_xlim(1990, 2023)
    ax.set_title(f"{'E' if j==0 else 'F'}  Incidence, {sex.lower()}s, 1990 = 100", loc="left")
fig.legend(*axes[0, 0].get_legend_handles_labels(), ncol=5, loc="upper left", bbox_to_anchor=(0.055, 0.952), frameon=False, fontsize=9)
fig.suptitle("Etiology-specific liver cancer rates by sex in China", x=0.06, ha="left", fontsize=13, fontweight="bold", y=1.0)
fig.text(0.06, 0.97, "Age-standardized rates, 1990-2023 (GBD 2023). Shaded bands: 95% UI for female HBV and HCV", fontsize=10, color="#1d1d1f99")
fig.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.05, wspace=0.15, hspace=0.3)
fig.savefig(BASE + "out/main/Fig2_China_ASR_by_sex_v2.5.png", dpi=200, bbox_inches="tight"); plt.close()
for sex in ["Male", "Female"]:
    for c in ["HBV", "HCV", "Alcohol", "MASH", "Other"]:
        g = get("China", c, "Incidence", sex); print(sex, c, round(100 * g.val.values[-1] / g.val.values[0], 1))
