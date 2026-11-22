"""Load official NCI Joinpoint 6.1.0 results (198 series, run 2026-11-08 by secretary).
Mapping id->series in id_map.json was independently verified by SE fingerprint (0 mismatches)."""
import json, numpy as np, pandas as pd
B = "out/nci_joinpoint/rerun_output/02_raw_exports/"
_m = {int(k): v for k, v in json.load(open(B + "id_map.json")).items()}
def _rd(f):
    d = pd.read_csv(B + f, sep="\t"); d["name"] = d.series.map(_m); return d
AAPC = _rd("all.Export.AAPC.txt"); APC = _rd("all.Export.APC.txt"); EST = _rd("all.Export.Model.Estimates.txt")
# ---- batch 2 (2026-11-10, secretary; rest-54 replaced 19:17 by v2 rerun with max 5 JP, v1 had max 6 -> archived *_v1_max6_SUPERSEDED): 72 JKUS series from co-author GBD 2023 download (male deaths, M/F/Both incidence); labelled exports
B2 = "out/nci_joinpoint/rerun2_output/"
def _rd2(kind):
    out = []
    for sub, stem in [("01_male_deaths_18_series/", "jkus_male_deaths"), ("02_other_54_series/", "jkus_rest_v2")]:
        d = pd.read_csv(B2 + sub + f"{stem}.Export.{kind}.txt", sep="\t")
        d["name"] = d.location + "_" + d.measure + "_" + d.sex + "_" + d.etiology; d["batch"] = 2; out.append(d)
    d = pd.concat(out, ignore_index=True)
    if "AAPC Index" in d: d["AAPC Index"] = d["AAPC Index"].map({"Full Range": 0, "Last 5": 1, "Last 10": 2})
    return d
_A2, _P2, _E2 = _rd2("AAPC"), _rd2("APC"), _rd2("Model.Estimates")
assert not set(_A2.name) & set(AAPC.name), "batch overlap"
AAPC = pd.concat([AAPC, _A2], ignore_index=True); APC = pd.concat([APC, _P2], ignore_index=True); EST = pd.concat([EST, _E2], ignore_index=True)
def key(loc, meas, sex, cause): return f"{loc}_{meas}_{sex}_{cause}"
def result(loc, meas, sex, cause, years=None):
    n = key(loc, meas, sex, cause)
    a = AAPC[(AAPC.name == n) & (AAPC["AAPC Index"] == 0)].iloc[0]
    A = (a["AAPC"], a["AAPC C.I. Low"], a["AAPC C.I. High"])
    s = APC[APC.name == n].sort_values("Segment")
    segs = [(int(r["Segment Start"]), int(r["Segment End"]), r["APC"], r["APC 95% LCL"], r["APC 95% UCL"]) for _, r in s.iterrows()]
    ks = [sg[0] for sg in segs[1:]]
    e = EST[EST.name == n].sort_values("Segment")
    fit = None
    if years is not None:
        years = np.asarray(years, float); fit = np.empty_like(years)
        bounds = [segs[0][0]] + ks + [segs[-1][1] + 1]
        for i, (_, r) in enumerate(e.iterrows()):
            msk = (years >= bounds[i]) & (years < bounds[i + 1]) if i < len(segs) - 1 else (years >= bounds[i])
            fit[msk] = np.exp(r["Intercept Estimate"] + r["Slope Estimate"] * years[msk])
    return ks, A, segs, fit, float(a["P-Value"])
