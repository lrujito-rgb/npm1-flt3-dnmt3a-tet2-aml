"""
figures.py — Publication figures from out/ tables (no values typed by hand).
Fig1 cohort flow + mutation landscape; Fig2 haematological phenotype; Fig3 Kaplan-Meier; Fig4 forest plots.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "publication.mplstyle"))
from matplotlib.patches import FancyBboxPatch, Patch
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts

OUT = os.environ.get("AML_OUT", "out")
COHORTS = ["TCGA", "BeatAML", "AMLSG"]
CLAB = {"TCGA": "TCGA (cBioPortal)", "BeatAML": "Beat AML (cBioPortal)", "AMLSG": "AMLSG"}
GROUPS = ["NPM1/ITD/epi-wt", "NPM1/ITD/epi-mut", "NPM1/TKD/epi-wt", "NPM1/TKD/epi-mut"]
COL = {"NPM1/ITD/epi-wt": "#F4A261", "NPM1/ITD/epi-mut": "#B5470F",
       "NPM1/TKD/epi-wt": "#8DB8DD", "NPM1/TKD/epi-mut": "#1D4E89"}
SHORT = {"NPM1/ITD/epi-wt": "ITD, DNMT3A/TET2 wt", "NPM1/ITD/epi-mut": "ITD, DNMT3A/TET2 mut",
         "NPM1/TKD/epi-wt": "TKD, DNMT3A/TET2 wt", "NPM1/TKD/epi-mut": "TKD, DNMT3A/TET2 mut"}
MARK = {"TCGA": "^", "BeatAML": "s", "AMLSG": "o"}


def load():
    df = pd.read_csv(f"{OUT}/analytic_dataset.csv")
    for c in ["primary", "itd", "tkd", "epi", "dnmt3a", "tet2", "dnmt3a_r882"]:
        df[c] = df[c].astype(str).str.lower().eq("true")
    return df, df[df.primary].copy()


def letter(ax, s, x=-0.12, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom", ha="left")


# ---------------------------------------------------------------- Figure 1
def fig1(df, p):
    fl = pd.read_csv(f"{OUT}/cohort_flow.csv")
    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.0], hspace=0.12)
    ax = fig.add_subplot(gs[0]); ax.set_axis_off(); ax.set_xlim(0, 3); ax.set_ylim(0, 1)
    letter(ax, "a", x=0.0, y=1.0)
    for i, c in enumerate(COHORTS):
        steps = fl[fl.cohort == c].reset_index(drop=True)
        keep = [s for s in steps.step if not s.strip().startswith("FLT3 status")]
        steps_main = steps[steps.step.isin(keep)]
        fs = steps[steps.step.str.strip().str.startswith("FLT3 status")]
        ax.text(i + 0.5, 0.985, CLAB[c], ha="center", va="top", fontweight="bold", fontsize=8)
        nmax = max(len([s for s in fl[fl.cohort == cc].step if not s.strip().startswith("FLT3 status")]) for cc in COHORTS)
        y = 0.93; dy = 0.92 / nmax
        for _, r in steps_main.iterrows():
            txt = f"{r.step.replace('OS available (primary analysis set)', 'Primary set (OS known)').replace('one diagnostic sample per patient', 'one sample per patient')}: n = {r.n}"
            if r.step == "NPM1 mutated":
                txt += "\n" + "; ".join(f"{s.strip().replace('FLT3 status: ', '')} {n}" for s, n in zip(fs.step, fs.n))
            box = FancyBboxPatch((i + 0.06, y - dy * 0.82), 0.88, dy * 0.78, boxstyle="round,pad=0.004",
                                 fc="#F2F2F2" if "primary" not in r.step else "#DDE8F3", ec="#7F7F7F", lw=0.6)
            ax.add_patch(box)
            ax.text(i + 0.5, y - dy * 0.43, txt, ha="center", va="center", fontsize=5.5)
            y -= dy
    # mutation landscape
    ax2 = fig.add_subplot(gs[1]); letter(ax2, "b", x=-0.08)
    q = p.copy()
    q["gi"] = q.group4.map({g: i for i, g in enumerate(GROUPS)}); q["ci"] = q.cohort.map({c: i for i, c in enumerate(COHORTS)})
    q = q.sort_values(["gi", "ci", "dnmt3a_r882", "dnmt3a", "tet2"], ascending=[True, True, False, False, False]).reset_index(drop=True)
    rows = [("Group", None), ("Cohort", None), ("FLT3-ITD", "itd"), ("FLT3-TKD (D835/I836)", "tkd"),
            ("DNMT3A", "dnmt3a"), ("TET2", "tet2")]
    ccol = {"TCGA": "#BDBDBD", "BeatAML": "#737373", "AMLSG": "#252525"}
    for j, (lab, col) in enumerate(rows):
        yy = len(rows) - 1 - j
        for k, r in q.iterrows():
            if lab == "Group":
                fc = COL[r.group4]
            elif lab == "Cohort":
                fc = ccol[r.cohort]
            elif col == "dnmt3a":
                fc = ("#5E3C99" if r.dnmt3a_r882 else "#B2ABD2") if r.dnmt3a else "#F0F0F0"
            elif col == "tet2":
                fc = "#5E3C99" if r.tet2 else "#F0F0F0"
            else:
                fc = "#404040" if r[col] else "#F0F0F0"
            ax2.add_patch(plt.Rectangle((k, yy + 0.08), 0.92, 0.84, fc=fc, lw=0))
    ax2.set_xlim(0, len(q)); ax2.set_ylim(0, len(rows))
    ax2.set_yticks([len(rows) - 1 - j + 0.5 for j in range(len(rows))]); ax2.set_yticklabels([r[0] for r in rows])
    for lab in ax2.get_yticklabels():
        if lab.get_text() in ("DNMT3A", "TET2"):
            lab.set_fontstyle("italic")
    ax2.set_xticks([]); [s.set_visible(False) for s in ax2.spines.values()]
    ax2.set_xlabel(f"Patients in the primary analysis set (n = {len(q)}), one column each")
    leg = [Patch(fc=COL[g], label=SHORT[g]) for g in GROUPS] + \
          [Patch(fc=ccol[c], label=CLAB[c]) for c in COHORTS] + \
          [Patch(fc="#5E3C99", label="DNMT3A R882 / TET2 mutated"), Patch(fc="#B2ABD2", label="DNMT3A non-R882"),
           Patch(fc="#404040", label="FLT3 alteration present"), Patch(fc="#F0F0F0", label="not detected")]
    ax2.legend(handles=leg, ncol=4, frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.12),
               handlelength=1.0, columnspacing=0.8)
    fig.savefig(f"{OUT}/Figure1_cohort_landscape.png", dpi=300, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- Figure 2
def fig2(p):
    mo = pd.read_csv(f"{OUT}/phenotype_models.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for ax, var, ylab, logy, lt in [(axes[0], "wbc", "White-cell count (×10⁹/L)", True, "a"),
                                     (axes[1], "bm_blast", "Bone-marrow blasts (%)", False, "b")]:
        for i, g in enumerate(GROUPS):
            s = p[(p.group4 == g) & p[var].notna()]
            bp = ax.boxplot(s[var], positions=[i], widths=0.55, showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.2), boxprops=dict(fc="none", ec="#404040", lw=0.7),
                            whiskerprops=dict(color="#404040", lw=0.7), capprops=dict(color="#404040", lw=0.7))
            jit = np.random.default_rng(i).uniform(-0.18, 0.18, len(s))
            for c in COHORTS:
                m = (s.cohort == c).values
                ax.scatter(i + jit[m], s[var].values[m], s=9, marker=MARK[c], color=COL[g], alpha=0.75,
                           edgecolor="none", zorder=3)
            ax.text(i, 1.0, f"n={len(s)}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6)
        if logy:
            ax.set_yscale("log"); ax.set_yticks([1, 10, 100]); ax.set_yticklabels(["1", "10", "100"])
        ax.set_xticks(range(4)); ax.set_xticklabels(["ITD\nwt", "ITD\nmut", "TKD\nwt", "TKD\nmut"], fontsize=6)
        ax.set_xlabel("FLT3 subtype / DNMT3A-TET2 status", fontsize=7)
        ax.set_ylabel(ylab); ax.margins(y=0.04); letter(ax, lt, y=1.08)
        r = mo[(mo.variable == var) & (mo.contrast == "TKD: epi-mut vs epi-wt (G4 vs G3)")].iloc[0]
        unit = "GMR" if var == "wbc" else "Δ"
        ax.set_title(f"TKD mut vs wt: {unit} {r.estimate:.2f} ({r.lo:.2f}–{r.hi:.2f})" if var == "wbc"
                     else f"TKD mut vs wt: Δ {r.estimate:.1f} pts ({r.lo:.1f} to {r.hi:.1f})", fontsize=7, loc="left", pad=12)
    h = [plt.Line2D([], [], marker=MARK[c], ls="", color="#606060", ms=4, label=CLAB[c]) for c in COHORTS]
    fig.legend(handles=h, frameon=False, fontsize=6, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(f"{OUT}/Figure2_phenotype.png", dpi=300, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- Figure 3
def fig3(p):
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.9), sharey=True)
    for ax, c, lt in zip(axes, COHORTS, "abc"):
        d = p[p.cohort == c]; fits = []
        for g in GROUPS:
            s = d[d.group4 == g]
            k = KaplanMeierFitter(label=SHORT[g].replace("DNMT3A/TET2 ", "")).fit(s.os_months, s.os_event)
            k.plot_survival_function(ax=ax, ci_show=False, color=COL[g], lw=1.3, show_censors=True,
                                     censor_styles={"ms": 3, "marker": "|"})
            fits.append(k)
        ax.set_xlim(0, 96); ax.set_xticks([0, 24, 48, 72, 96]); ax.set_ylim(0, 1.02)
        ax.set_title(f"{CLAB[c]} (n={len(d)})", fontsize=8, loc="left")
        ax.set_xlabel("Months from diagnosis"); letter(ax, lt, x=-0.08)
        ax.get_legend().remove() if ax.get_legend() else None
        add_at_risk_counts(*fits, ax=ax, rows_to_show=["At risk"], xticks=[0, 24, 48, 72, 96], fontsize=5)
    axes[0].set_ylabel("Overall survival probability")
    h = [plt.Line2D([], [], color=COL[g], lw=1.5, label=SHORT[g]) for g in GROUPS]
    fig.legend(handles=h, ncol=4, frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(f"{OUT}/Figure3_KM_OS.png", dpi=300, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- Figure 4
def fig4():
    T = pd.read_csv(f"{OUT}/Table3c_cox_models.csv")
    import json
    res = json.load(open(f"{OUT}/results.json"))
    fig, axes = plt.subplots(2, 1, figsize=(6.0, 8.6), gridspec_kw=dict(height_ratios=[1.0, 1.05], hspace=0.32))
    # (a) stratum-specific HRs
    ax = axes[0]
    models = [("cohort-specific: TCGA, adj. age+sex", "TCGA"), ("cohort-specific: BeatAML, adj. age+sex", "Beat AML"),
              ("cohort-specific: AMLSG, adj. age+sex", "AMLSG"), ("PRIMARY: pooled, cohort-stratified, adj. age+sex", "Pooled (stratified)")]
    rows = []
    for con, lab, col in [("HR ITD vs TKD | epi-wt", "ITD vs TKD in DNMT3A/TET2 wt", COL["NPM1/TKD/epi-wt"]),
                          ("HR ITD vs TKD | epi-mut", "ITD vs TKD in DNMT3A/TET2 mut", COL["NPM1/TKD/epi-mut"]),
                          ("Ratio of HRs (interaction)", "Ratio of HRs (mut / wt)", "#404040")]:
        rows.append(("header", lab, None, None))
        for m, ml in models:
            r = T[(T.model == m) & (T.contrast == con)].iloc[0]
            rows.append((ml, r, col, m.startswith("PRIMARY")))
    y = len(rows)
    for item in rows:
        y -= 1
        if item[0] == "header":
            ax.text(0.02, y, item[1], transform=ax.get_yaxis_transform(), fontsize=6.5, fontweight="bold", va="center")
            continue
        ml, r, col, prim = item
        ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.0)
        ax.plot(r.HR, y, marker="D" if prim else "s", color=col, ms=5 if prim else 3.5)
        ax.text(1.02, y, f"{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f})", transform=ax.get_yaxis_transform(), fontsize=5.5, va="center")
        ax.text(-0.02, y, ml, transform=ax.get_yaxis_transform(), fontsize=5.5, va="center", ha="right")
    ax.axvline(1, color="#7F7F7F", lw=0.6, ls="--"); ax.set_xscale("log")
    ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10"]); ax.set_xlim(0.08, 15)
    ax.set_yticks([]); ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("Hazard ratio for death (95% CI), log scale"); letter(ax, "a", x=-0.45, y=1.0)
    # (b) interaction across sensitivity analyses
    ax = axes[1]
    sens = T[(T.contrast == "Ratio of HRs (interaction)") & (T.model.str.startswith(("PRIMARY", "S", "inverse", "pooled, cohort-stratified, unadj")))]
    labs = []
    for i, (_, r) in enumerate(sens.iloc[::-1].iterrows()):
        prim = r.model.startswith("PRIMARY")
        ax.plot([r.lo, r.hi], [i, i], color="#404040", lw=1.0)
        ax.plot(r.HR, i, marker="D" if prim else "s", color="#1D4E89" if prim else "#404040", ms=5 if prim else 3.5)
        labs.append(r.model.replace("PRIMARY: ", "Primary: ").replace("inverse-variance pooling of cohort-specific estimates", "IV pooling").replace(", adj. age+sex", ""))
        ax.text(1.02, i, f"{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f})", transform=ax.get_yaxis_transform(), fontsize=5.5, va="center")
    b = res["S7_bootstrap_interaction"]; i = len(labs)
    ax.plot([b["lo"], b["hi"]], [i, i], color="#404040", lw=1.0); ax.plot(b["ratio_median"], i, "s", color="#404040", ms=3.5)
    ax.text(1.02, i, f"{b['ratio_median']:.2f} ({b['lo']:.2f}–{b['hi']:.2f})", transform=ax.get_yaxis_transform(), fontsize=5.5, va="center")
    labs.append("S7: bootstrap percentile CI (2,000 replicates)")
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=5.5)
    ax.axvline(1, color="#7F7F7F", lw=0.6, ls="--"); ax.set_xscale("log")
    ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10"]); ax.set_xlim(0.08, 15)
    ax.set_xlabel("Ratio of hazard ratios (95% CI), log scale"); letter(ax, "b", x=-0.62, y=1.02)
    ax.set_title("Interaction estimate is imprecise in every specification", fontsize=7, loc="left")
    fig.savefig(f"{OUT}/Figure4_forest.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def main():
    df, p = load()
    fig1(df, p); fig2(p); fig3(p); fig4()


if __name__ == "__main__":
    main()
