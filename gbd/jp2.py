"""Joinpoint v2, closer to NCI Joinpoint defaults for annual data (n=34):
max 5 joinpoints, >=2 obs from ends and between joinpoints (joinpoints at integer years),
log-linear, heteroscedastic WLS with SE of ln(ASR) from GBD 95% UI, BIC model selection.
AAPC = length-weighted sum of segment slopes; CI via t (df = n - params)."""
import itertools, numpy as np
from scipy import stats

def design(x, ks):
    return np.column_stack([np.ones_like(x), x] + [np.clip(x - k, 0, None) for k in ks])

def jp(x, y, se=None, maxk=5, gap=2):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    w = np.ones(n) if se is None else 1 / np.asarray(se, float) ** 2
    sw = np.sqrt(w); yw = y * sw
    cand = x[gap:n - gap]  # interior years allowed
    best = None
    for nk in range(maxk + 1):
        for ks in itertools.combinations(cand, nk):
            if any(ks[i + 1] - ks[i] < gap for i in range(nk - 1)): continue
            X = design(x, ks); Xw = X * sw[:, None]
            beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
            rss = float(((yw - Xw @ beta) ** 2).sum())
            p = 2 + 2 * nk
            bic = np.log(rss / n) + p * np.log(n) / n
            if best is None or bic < best[0]: best = (bic, ks, beta, rss, X, Xw)
    bic, ks, beta, rss, X, Xw = best
    p = X.shape[1]; df = n - p - len(ks); s2 = rss / df
    cov = s2 * np.linalg.inv(Xw.T @ Xw); tq = stats.t.ppf(0.975, df)
    f = lambda v: 100 * (np.exp(v) - 1)
    bounds = [x[0], *ks, x[-1]]; segs = []; L = []
    for i in range(len(bounds) - 1):
        c = np.zeros(p); c[1:2 + i] = 1; sl = c @ beta; s = np.sqrt(c @ cov @ c)
        segs.append((int(bounds[i]), int(bounds[i + 1]), f(sl), f(sl - tq * s), f(sl + tq * s))); L.append((bounds[i + 1] - bounds[i], c))
    tot = sum(l for l, _ in L); cA = sum(l / tot * c for l, c in L); sl = cA @ beta; s = np.sqrt(cA @ cov @ cA)
    return tuple(int(k) for k in ks), (f(sl), f(sl - tq * s), f(sl + tq * s)), segs, beta
