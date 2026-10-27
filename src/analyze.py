"""
analyze.py — Pre-specified analyses (protocol locked 2026-10-26).
Input : out/analytic_dataset.csv (from build_data.py)
Output: out/*.csv tables, out/results.json
"""
import os, json, warnings, itertools
import numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import median_survival_times
import statsmodels.formula.api as smf
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

OUT = os.environ.get("AML_OUT", "out")
COHORTS = ["TCGA", "BeatAML", "AMLSG"]
GROUPS = ["NPM1/ITD/epi-wt", "NPM1/ITD/epi-mut", "NPM1/TKD/epi-wt", "NPM1/TKD/epi-mut"]
GLAB = {"NPM1/ITD/epi-wt": "G1", "NPM1/ITD/epi-mut": "G2", "NPM1/TKD/epi-wt": "G3", "NPM1/TKD/epi-mut": "G4"}
MIN_EVENTS_PER_CELL = 3
warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261026)


class AnalysisError(RuntimeError):
    pass


def load():
    df = pd.read_csv(f"{OUT}/analytic_dataset.csv")
    for c in ["npm1", "itd", "tkd", "epi", "dnmt3a", "tet2", "epi_restr", "dnmt3a_restr", "tet2_restr", "primary",
              "itd_strict", "tkd_ext", "dnmt3a_r882"]:
        df[c] = df[c].astype(str).str.lower().eq("true")
    p = df[df.primary].copy()
    if p.empty:
        raise AnalysisError("empty primary set")
    unexpected = set(p.group4) - set(GROUPS)
    if unexpected:
        raise AnalysisError(f"unexpected group labels {unexpected}")
    return df, p


def fmt_p(p):
    return "NA" if p is None or np.isnan(p) else (f"{p:.3f}" if p >= 0.001 else "<0.001")


def iqr(x):
    x = pd.Series(x).dropna()
    return (np.nan, np.nan, np.nan, 0) if x.empty else (x.median(), x.quantile(.25), x.quantile(.75), len(x))


# =========================================================================== Table 1
def mc_chi2(groups, cats, B=20000):
    """Monte-Carlo permutation p-value for an r x c chi-square (small expected counts)."""
    g = np.asarray(groups); c = np.asarray(cats)
    ok = ~pd.isna(c); g, c = g[ok], c[ok]
    if len(np.unique(c)) < 2:
        return np.nan
    tab = pd.crosstab(g, c).values
    obs = stats.chi2_contingency(tab, correction=False)[0]
    cnt = 0
    for _ in range(B // 10):
        perm = rng.permutation(c)
        if stats.chi2_contingency(pd.crosstab(g, perm).values, correction=False)[0] >= obs - 1e-12:
            cnt += 1
    return (cnt + 1) / (B // 10 + 1)


def table1(p):
    rows = []
    def add(label, fn, test=None):
        r = {"characteristic": label}
        for gname in GROUPS:
            r[GLAB[gname] + " " + gname] = fn(p[p.group4 == gname])
        r["all"] = fn(p)
        r["p"] = test
        rows.append(r)
    add("n", lambda d: str(len(d)))
    for c in COHORTS:
        add(f"  cohort {c}, n (%)", lambda d, c=c: f"{(d.cohort == c).sum()} ({(d.cohort == c).mean()*100:.0f}%)")
    def med(col, dig=1):
        def f(d):
            m, q1, q3, n = iqr(d[col])
            return "NA" if n == 0 else f"{m:.{dig}f} ({q1:.{dig}f}-{q3:.{dig}f}); n={n}"
        return f
    kw = lambda col: fmt_p(stats.kruskal(*[p.loc[p.group4 == g, col].dropna() for g in GROUPS])[1])
    add("Age, y, median (IQR)", med("age", 0), kw("age"))
    add("Male, n/N (%)", lambda d: f"{int(d.male.sum())}/{d.male.notna().sum()} ({d.male.mean()*100:.0f}%)",
        fmt_p(mc_chi2(p.group4, p.male)))
    add("WBC x10^9/L, median (IQR)", med("wbc"), kw("wbc"))
    add("BM blasts %, median (IQR)", med("bm_blast", 0), kw("bm_blast"))
    add("PB blasts %, median (IQR)", med("pb_blast", 0), kw("pb_blast"))
    add("Normal karyotype, n/N (%)", lambda d: f"{int(d.nk.sum())}/{d.nk.notna().sum()} ({d.nk.mean()*100:.0f}%)",
        fmt_p(mc_chi2(p.group4, p.nk)))
    add("Adverse cytogenetics (TCGA, AMLSG), n/N", lambda d: f"{int(d.cyto_adverse.sum())}/{d.cyto_adverse.notna().sum()}",
        fmt_p(mc_chi2(p.group4, p.cyto_adverse)))
    add("De novo AML, n/N (%)", lambda d: f"{int(d.de_novo.sum())}/{d.de_novo.notna().sum()} ({d.de_novo.mean()*100:.0f}%)",
        fmt_p(mc_chi2(p.group4, p.de_novo)))
    add("Any HSCT recorded (TCGA, AMLSG), n/N", lambda d: f"{int(d.hsct.sum())}/{d.hsct.notna().sum()}",
        fmt_p(mc_chi2(p.group4, p.hsct)))
    add("DNMT3A mutated, n", lambda d: str(int(d.dnmt3a.sum())))
    add("  DNMT3A R882, n", lambda d: str(int(d.dnmt3a_r882.sum())))
    add("TET2 mutated, n", lambda d: str(int(d.tet2.sum())))
    add("DNMT3A and TET2 mutated, n", lambda d: str(int((d.dnmt3a & d.tet2).sum())))
    add("Deaths, n", lambda d: str(int(d.os_event.sum())))
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/Table1_baseline.csv", index=False)
    return t


# =========================================================================== Aim 1 co-mutation
def or_ci(tab):
    r = stats.contingency.odds_ratio(tab, kind="conditional")
    ci = r.confidence_interval(0.95)
    return r.statistic, ci.low, ci.high, stats.fisher_exact(tab)[1]


def aim1(p):
    rows = []
    for gene, col in [("DNMT3A", "dnmt3a"), ("TET2", "tet2"), ("DNMT3A and/or TET2", "epi")]:
        tabs = []
        for c in COHORTS + ["pooled"]:
            d = p if c == "pooled" else p[p.cohort == c]
            i, t = d[d.flt3_group == "ITD"], d[d.flt3_group == "TKD"]
            tab = np.array([[i[col].sum(), (~i[col]).sum()], [t[col].sum(), (~t[col]).sum()]], dtype=int)
            if c != "pooled":
                tabs.append(tab)
                o, lo, hi, pv = or_ci(tab)
                rows.append(dict(gene=gene, cohort=c, itd_mut=tab[0, 0], itd_n=tab[0].sum(),
                                 itd_pct=100 * tab[0, 0] / tab[0].sum(), tkd_mut=tab[1, 0], tkd_n=tab[1].sum(),
                                 tkd_pct=100 * tab[1, 0] / tab[1].sum(), method="Fisher exact (conditional MLE OR)",
                                 OR=o, lo=lo, hi=hi, p=pv))
            else:
                st = StratifiedTable([t_.astype(float) + (0.5 if (t_ == 0).any() else 0) for t_ in tabs])
                lo, hi = st.oddsratio_pooled_confint()
                rows.append(dict(gene=gene, cohort="pooled (Mantel-Haenszel, cohort-stratified)",
                                 itd_mut=tab[0, 0], itd_n=tab[0].sum(), itd_pct=100 * tab[0, 0] / tab[0].sum(),
                                 tkd_mut=tab[1, 0], tkd_n=tab[1].sum(), tkd_pct=100 * tab[1, 0] / tab[1].sum(),
                                 method="Mantel-Haenszel (0.5 added to tables with a zero cell)",
                                 OR=st.oddsratio_pooled, lo=lo, hi=hi, p=st.test_null_odds(correction=True).pvalue))
    t = pd.DataFrame(rows)
    pooled = t.cohort.str.startswith("pooled")
    t.loc[pooled, "p_BH"] = multipletests(t.loc[pooled, "p"], method="fdr_bh")[1]
    t.to_csv(f"{OUT}/Table2_comutation.csv", index=False)
    return t


# =========================================================================== Aim 2 phenotype
def dunn(values, groups):
    d = pd.DataFrame({"v": values, "g": groups}).dropna()
    d["r"] = d.v.rank()
    N = len(d); ties = d.v.value_counts(); T = ((ties ** 3 - ties).sum()) / (12 * (N - 1))
    out = []
    for a, b in itertools.combinations(GROUPS, 2):
        ra, rb = d.loc[d.g == a, "r"], d.loc[d.g == b, "r"]
        if len(ra) == 0 or len(rb) == 0:
            continue
        se = np.sqrt((N * (N + 1) / 12 - T) * (1 / len(ra) + 1 / len(rb)))
        z = (ra.mean() - rb.mean()) / se
        out.append((GLAB[a] + " vs " + GLAB[b], z, 2 * stats.norm.sf(abs(z))))
    r = pd.DataFrame(out, columns=["contrast", "z", "p"])
    if len(r):
        r["p_BH"] = multipletests(r.p, method="fdr_bh")[1]
    return r


def aim2(p):
    desc, posthoc, models = [], [], []
    for var in ["wbc", "bm_blast", "pb_blast"]:
        for c in COHORTS + ["pooled"]:
            d = p if c == "pooled" else p[p.cohort == c]
            for g in GROUPS:
                m, q1, q3, n = iqr(d.loc[d.group4 == g, var])
                desc.append(dict(variable=var, cohort=c, group=g, n=n, median=m, q1=q1, q3=q3))
            if c != "pooled":
                samples = [d.loc[d.group4 == g, var].dropna() for g in GROUPS]
                samples = [s for s in samples if len(s) > 0]
                kwp = stats.kruskal(*samples)[1] if len(samples) > 1 else np.nan
                ph = dunn(d[var], d.group4); ph["variable"] = var; ph["cohort"] = c; ph["kruskal_p"] = kwp
                posthoc.append(ph)
        # pooled linear model with cohort fixed effects, HC3
        d = p[[var, "itd", "epi", "cohort"]].dropna().copy()
        d["y"] = np.log(d[var]) if var == "wbc" else d[var]
        if var == "wbc" and (d[var] <= 0).any():
            raise AnalysisError("non-positive WBC")
        d["itd"] = d.itd.astype(int); d["epi"] = d.epi.astype(int)
        fit = smf.ols("y ~ itd * epi + C(cohort)", data=d).fit(cov_type="HC3")
        b, V = fit.params, fit.cov_params()
        L = {"TKD: epi-mut vs epi-wt (G4 vs G3)": {"epi": 1},
             "ITD: epi-mut vs epi-wt (G2 vs G1)": {"epi": 1, "itd:epi": 1},
             "epi-wt: ITD vs TKD (G1 vs G3)": {"itd": 1},
             "epi-mut: ITD vs TKD (G2 vs G4)": {"itd": 1, "itd:epi": 1},
             "G4 (TKD epi-mut) vs G1 (ITD epi-wt)": {"epi": 1, "itd": -1},
             "interaction (difference of differences)": {"itd:epi": 1}}
        for name, w in L.items():
            vec = pd.Series(0.0, index=b.index); [vec.__setitem__(k, v) for k, v in w.items()]
            est = float(vec @ b); se = float(np.sqrt(vec @ V @ vec)); z = est / se
            lo, hi = est - 1.96 * se, est + 1.96 * se
            if var == "wbc":
                est, lo, hi = np.exp(est), np.exp(lo), np.exp(hi)
            models.append(dict(variable=var, contrast=name, scale="geometric-mean ratio" if var == "wbc" else "mean difference (%)",
                               estimate=est, lo=lo, hi=hi, p=2 * stats.norm.sf(abs(z)), n=int(fit.nobs)))
    D, PH, MO = pd.DataFrame(desc), pd.concat(posthoc), pd.DataFrame(models)
    D.to_csv(f"{OUT}/phenotype_descriptives.csv", index=False)
    PH.to_csv(f"{OUT}/phenotype_kruskal_dunn.csv", index=False)
    MO.to_csv(f"{OUT}/phenotype_models.csv", index=False)
    return D, PH, MO


# =========================================================================== Aim 3 KM / log-rank
def km_table(p, time="os_months", event="os_event"):
    rows = []
    for c in COHORTS + ["pooled (descriptive)"]:
        d = p if c.startswith("pooled") else p[p.cohort == c]
        for g in GROUPS:
            s = d[d.group4 == g]
            if len(s) == 0:
                continue
            k = KaplanMeierFitter().fit(s[time], s[event])
            med = k.median_survival_time_
            ci = median_survival_times(k.confidence_interval_)
            sf = k.survival_function_at_times(36).iloc[0]
            cis = k.confidence_interval_survival_function_
            idx = cis.index[cis.index <= 36]
            lo36, hi36 = (cis.loc[idx[-1]].values if len(idx) else (np.nan, np.nan))
            rows.append(dict(cohort=c, group=g, n=len(s), events=int(s[event].sum()),
                             median_os=med, median_lo=ci.iloc[0, 0], median_hi=ci.iloc[0, 1],
                             os36=sf, os36_lo=lo36, os36_hi=hi36,
                             max_followup=s[time].max()))
    return pd.DataFrame(rows)


def stratified_logrank(d, grp_a, grp_b, time="os_months", event="os_event"):
    d = d[d.group4.isin([grp_a, grp_b])]
    O = E = V = 0.0
    for _, s in d.groupby("cohort"):
        x = (s.group4 == grp_a).values.astype(float); t = s[time].values; e = s[event].values
        for tt in np.unique(t[e == 1]):
            at = t >= tt; n = at.sum(); n1 = (at & (x == 1)).sum()
            dd = ((t == tt) & (e == 1)).sum(); d1 = ((t == tt) & (e == 1) & (x == 1)).sum()
            if n < 2:
                continue
            O += d1; E += dd * n1 / n; V += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    chi = (O - E) ** 2 / V if V > 0 else np.nan
    return O, E, chi, stats.chi2.sf(chi, 1)


def aim3(p):
    km = km_table(p); km.to_csv(f"{OUT}/Table3a_km.csv", index=False)
    rows = []
    for a, b, lab in [("NPM1/ITD/epi-wt", "NPM1/TKD/epi-wt", "ITD vs TKD within epi-wt (G1 vs G3)"),
                      ("NPM1/ITD/epi-mut", "NPM1/TKD/epi-mut", "ITD vs TKD within epi-mut (G2 vs G4)"),
                      ("NPM1/TKD/epi-mut", "NPM1/TKD/epi-wt", "epi-mut vs epi-wt within TKD (G4 vs G3)"),
                      ("NPM1/ITD/epi-mut", "NPM1/ITD/epi-wt", "epi-mut vs epi-wt within ITD (G2 vs G1)"),
                      ("NPM1/TKD/epi-mut", "NPM1/ITD/epi-wt", "G4 vs G1 (draft 'equalization' contrast; descriptive)")]:
        for c in COHORTS + ["stratified"]:
            d = p if c == "stratified" else p[p.cohort == c]
            O, E, chi, pv = stratified_logrank(d, a, b)
            rows.append(dict(contrast=lab, cohort=c, observed_first=O, expected_first=E, chi2=chi, p=pv))
    lr = pd.DataFrame(rows); lr.to_csv(f"{OUT}/Table3b_logrank.csv", index=False)
    return km, lr


# =========================================================================== Cox models
CONTRASTS = {"HR ITD vs TKD | epi-wt": {"itd": 1},
             "HR ITD vs TKD | epi-mut": {"itd": 1, "itd_x_epi": 1},
             "HR epi-mut vs wt | TKD": {"epi": 1},
             "HR epi-mut vs wt | ITD": {"epi": 1, "itd_x_epi": 1},
             "Ratio of HRs (interaction)": {"itd_x_epi": 1}}


def design(d, epi_col="epi", itd_col="itd", extra=()):
    x = pd.DataFrame({"T": d.os_months.values, "E": d.os_event.values.astype(int),
                      "itd": d[itd_col].astype(int).values, "epi": d[epi_col].astype(int).values,
                      "age10": d.age10.values, "male": d.male.values, "cohort": d.cohort.values})
    x["itd_x_epi"] = x.itd * x.epi
    for e in extra:
        x[e] = d[e].values
    return x.dropna()


def fit_cox(x, strata=True, covars=("age10", "male"), time="T", event="E"):
    cols = [time, event, "itd", "epi", "itd_x_epi", *covars] + (["cohort"] if strata else [])
    cph = CoxPHFitter()
    cph.fit(x[cols], duration_col=time, event_col=event, strata=["cohort"] if strata else None)
    return cph


def contrasts(cph):
    b, V = cph.params_, cph.variance_matrix_
    out = {}
    for name, w in CONTRASTS.items():
        vec = pd.Series(0.0, index=b.index)
        for k, v in w.items():
            vec[k] = v
        est = float(vec @ b); se = float(np.sqrt(vec @ V.values @ vec))
        out[name] = dict(HR=np.exp(est), lo=np.exp(est - 1.96 * se), hi=np.exp(est + 1.96 * se),
                         p=2 * stats.norm.sf(abs(est / se)), logHR=est, se=se)
    return out


def cell_events(x):
    return x.groupby(["itd", "epi"]).E.sum().to_dict()


def cox_row(label, x, strata=True, covars=("age10", "male")):
    ce = cell_events(x)
    cph = fit_cox(x, strata, covars)
    c = contrasts(cph)
    rows = []
    for k, v in c.items():
        rows.append(dict(model=label, contrast=k, n=len(x), events=int(x.E.sum()),
                         min_cell_events=int(min(ce.values())) if len(ce) == 4 else 0,
                         sparse_flag=(len(ce) < 4 or min(ce.values()) < MIN_EVENTS_PER_CELL), **v))
    return rows, cph


def pool_iv(est, se):
    est, se = np.asarray(est), np.asarray(se); w = 1 / se ** 2
    fe = np.sum(w * est) / w.sum(); fe_se = np.sqrt(1 / w.sum())
    Q = np.sum(w * (est - fe) ** 2); dfree = len(est) - 1
    tau2 = max(0.0, (Q - dfree) / (w.sum() - np.sum(w ** 2) / w.sum())) if dfree > 0 else 0.0
    wr = 1 / (se ** 2 + tau2); re = np.sum(wr * est) / wr.sum(); re_se = np.sqrt(1 / wr.sum())
    I2 = max(0.0, (Q - dfree) / Q) * 100 if Q > 0 else 0.0
    return dict(fe=fe, fe_se=fe_se, re=re, re_se=re_se, Q=Q, Q_p=stats.chi2.sf(Q, dfree) if dfree > 0 else np.nan,
                I2=I2, tau2=tau2)


def bootstrap_interaction(x, B=2000):
    vals, fails = [], 0
    groups = {c: g.index.values for c, g in x.groupby("cohort")}
    for _ in range(B):
        idx = np.concatenate([rng.choice(ix, size=len(ix), replace=True) for ix in groups.values()])
        xb = x.loc[idx].reset_index(drop=True)
        try:
            cph = fit_cox(xb)
            v = cph.params_["itd_x_epi"]
            if np.isfinite(v) and abs(v) < 15:
                vals.append(v)
            else:
                fails += 1
        except Exception:
            fails += 1
    vals = np.array(vals)
    return dict(B=B, successful=len(vals), failed_or_divergent=fails,
                ratio_median=float(np.exp(np.median(vals))),
                lo=float(np.exp(np.percentile(vals, 2.5))), hi=float(np.exp(np.percentile(vals, 97.5))))


def primary_and_sensitivity(df, p):
    rows, extra = [], {}
    x = design(p)
    r, cph = cox_row("PRIMARY: pooled, cohort-stratified, adj. age+sex", x); rows += r
    # PH diagnostics
    ph = proportional_hazard_test(cph, x[["T", "E", "itd", "epi", "itd_x_epi", "age10", "male", "cohort"]],
                                  time_transform="rank").summary.reset_index().rename(columns={"index": "term"})
    ph.to_csv(f"{OUT}/cox_ph_schoenfeld.csv", index=False)
    extra["ph_global_min_p"] = float(ph["p"].min())
    extra["primary_coefficients"] = cph.summary[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].reset_index().to_dict("records")
    extra["events_per_parameter"] = float(x.E.sum() / 5)
    extra["cell_events"] = {f"itd={k[0]},epi={k[1]}": int(v) for k, v in cell_events(x).items()}
    # unadjusted
    r, _ = cox_row("pooled, cohort-stratified, unadjusted", x, covars=()); rows += r
    # cohort-specific
    per = {}
    for c in COHORTS:
        xc = x[x.cohort == c]
        try:
            r, cc = cox_row(f"cohort-specific: {c}, adj. age+sex", xc, strata=False); rows += r
            per[c] = contrasts(cc)
        except Exception as e:
            rows.append(dict(model=f"cohort-specific: {c}", contrast="FAILED", note=str(e)[:120]))
    # IV pooling of cohort-specific interactions
    ks = [c for c in COHORTS if c in per]
    pool = {}
    for k in CONTRASTS:
        pl = pool_iv([per[c][k]["logHR"] for c in ks], [per[c][k]["se"] for c in ks])
        pool[k] = pl
        for kind in ["fe", "re"]:
            rows.append(dict(model=f"inverse-variance pooling of cohort-specific estimates ({'fixed' if kind=='fe' else 'DerSimonian-Laird random'})",
                             contrast=k, n=len(x), events=int(x.E.sum()), HR=np.exp(pl[kind]),
                             lo=np.exp(pl[kind] - 1.96 * pl[kind + "_se"]), hi=np.exp(pl[kind] + 1.96 * pl[kind + "_se"]),
                             p=2 * stats.norm.sf(abs(pl[kind] / pl[kind + "_se"])), I2=pl["I2"], Q_p=pl["Q_p"]))
    extra["heterogeneity"] = {k: dict(I2=v["I2"], Q_p=v["Q_p"]) for k, v in pool.items()}
    # S1 restricted definitions
    r, _ = cox_row("S1: restricted variant definition", design(p, epi_col="epi_restr")); rows += r
    # S2 DNMT3A-only modifier
    r, _ = cox_row("S2: modifier = DNMT3A (any), irrespective of TET2", design(p, epi_col="dnmt3a")); rows += r
    # S3 + normal karyotype
    xs = design(p, extra=("nk",))
    r, _ = cox_row("S3: + normal karyotype (complete cases)", xs, covars=("age10", "male", "nk")); rows += r
    # S4 leave-one-cohort-out
    for c in COHORTS:
        r, _ = cox_row(f"S4: leave out {c}", x[x.cohort != c]); rows += r
    # S5 TCGA strict ITD / extended TKD
    alt = df[df.npm1].copy()
    alt["itd2"] = np.where(alt.cohort == "TCGA", alt.itd_strict, alt.itd)
    alt["tkd2"] = alt.tkd_ext
    alt["fg2"] = np.select([alt.itd2 & alt.tkd2, alt.itd2, alt.tkd2], ["ITD+TKD", "ITD", "TKD"], "none")
    s5 = alt[alt.fg2.isin(["ITD", "TKD"]) & alt.os_months.notna() & alt.os_event.notna()].copy()
    s5["itd"] = s5.fg2 == "ITD"
    r, _ = cox_row("S5: TCGA strict ITD + extended TKD (all cohorts)", design(s5)); rows += r
    extra["S5_n_by_group"] = s5.groupby(["cohort", "fg2"]).size().astype(int).to_dict().__repr__()
    # S6 ITD+TKD counted as ITD
    s6 = df[df.npm1 & df.flt3_group.isin(["ITD", "TKD", "ITD+TKD"]) & df.os_months.notna() & df.os_event.notna()].copy()
    s6["itd"] = s6.flt3_group.isin(["ITD", "ITD+TKD"])
    r, _ = cox_row("S6: ITD+TKD patients counted as ITD", design(s6)); rows += r
    # S8 EFS (AMLSG)
    a = p[p.cohort == "AMLSG"].copy()
    xe = design(a); xe["T"] = a.efs_months.clip(lower=0.01).values; xe["E"] = a.efs_event.astype(int).values
    r, _ = cox_row("S8: EFS, AMLSG only, adj. age+sex", xe, strata=False); rows += r
    # S9 de novo only
    r, _ = cox_row("S9: de novo AML only", design(p[p.de_novo == 1])); rows += r
    # S7 bootstrap
    extra["S7_bootstrap_interaction"] = bootstrap_interaction(x)
    T = pd.DataFrame(rows); T.to_csv(f"{OUT}/Table3c_cox_models.csv", index=False)
    return T, extra


def vaf_summary(p):
    rows = []
    for g in ["dnmt3a", "tet2"]:
        for c in COHORTS:
            for f in ["ITD", "TKD"]:
                v = p.loc[(p.cohort == c) & (p.flt3_group == f) & p[g], g + "_vaf"].dropna()
                m, q1, q3, n = iqr(v)
                rows.append(dict(gene=g.upper(), cohort=c, flt3=f, n_with_vaf=n, median_vaf=m, q1=q1, q3=q3))
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/Supp_vaf.csv", index=False)
    return t


def main():
    df, p = load()
    res = {"n_primary": len(p), "deaths": int(p.os_event.sum()),
           "n_by_cohort": p.cohort.value_counts().to_dict(),
           "group_counts": p.groupby(["cohort", "group4"]).size().astype(int).to_dict().__repr__()}
    t1 = table1(p); t2 = aim1(p); D, PH, MO = aim2(p); km, lr = aim3(p)
    T, extra = primary_and_sensitivity(df, p); res.update(extra); vaf_summary(p)
    res["followup_reverse_km_median_months"] = {c: float(KaplanMeierFitter().fit(g.os_months, 1 - g.os_event).median_survival_time_)
                                                for c, g in p.groupby("cohort")}
    json.dump(res, open(f"{OUT}/results.json", "w"), indent=2, default=str)
    return res


if __name__ == "__main__":
    r = main()
    print(json.dumps({k: r[k] for k in ["n_primary", "deaths", "n_by_cohort", "events_per_parameter",
                                         "cell_events", "ph_global_min_p", "S7_bootstrap_interaction"]}, indent=1, default=str))
