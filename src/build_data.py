"""
build_data.py — Data import, validation, harmonization and cohort assignment.

Cohorts
  TCGA      cBioPortal laml_tcga_pub         (sequence-derived FLT3-ITD; validated classifier)
  BeatAML   cBioPortal aml_ohsu_2022         (clinical PCR/CE FLT3_ITD; WES TKD/DNMT3A/TET2)
  AMLSG     github gerstung-lab/AML-multistage/data (clinical ITD/TKD; 111-gene panel)

Run via the repository entry point:  python run_all.py   (or set AML_RAW / AML_OUT and run this file directly)
Writes out/analytic_dataset.csv, out/cohort_flow.csv, out/data_dictionary.csv.
All column names below were verified against the downloaded files on 2026-10-26.
"""
import os, re, json, urllib.request, urllib.parse
import numpy as np, pandas as pd

RAW = os.environ.get("AML_RAW", "raw")
OUT = os.environ.get("AML_OUT", "out")
os.makedirs(RAW, exist_ok=True); os.makedirs(OUT, exist_ok=True)
API = "https://www.cbioportal.org/api"
GENES = {"NPM1": 4869, "FLT3": 2322, "DNMT3A": 1788, "TET2": 54790}
PROTEIN_ALTERING = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                    "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site"}
TRUNCATING = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation"}


class DataError(RuntimeError):
    pass


def require(df, cols, where):
    miss = [c for c in cols if c not in df.columns]
    if miss:
        raise DataError(f"{where}: required columns missing: {miss}")


# --------------------------------------------------------------------------- download helpers
def _get(path, params=None):
    u = API + path + ("?" + urllib.parse.urlencode(params) if params else "")
    with urllib.request.urlopen(urllib.request.Request(u, headers={"Accept": "application/json"}), timeout=300) as r:
        return json.loads(r.read())


def _post(path, body, params=None):
    u = API + path + ("?" + urllib.parse.urlencode(params) if params else "")
    req = urllib.request.Request(u, data=json.dumps(body).encode(),
                                 headers={"Accept": "application/json", "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.loads(r.read())


def fetch_cbioportal(study):
    """Download clinical (sample+patient), sequenced sample list and target-gene mutations (skips if cached)."""
    f_s = f"{RAW}/{study}_clinical_sample.tsv"
    if os.path.exists(f_s):
        return
    sc = pd.DataFrame(_get(f"/studies/{study}/clinical-data", {"clinicalDataType": "SAMPLE", "pageSize": 10**7}))
    pc = pd.DataFrame(_get(f"/studies/{study}/clinical-data", {"clinicalDataType": "PATIENT", "pageSize": 10**7}))
    sc.pivot_table(index=["sampleId", "patientId"], columns="clinicalAttributeId", values="value",
                   aggfunc="first").reset_index().to_csv(f_s, sep="\t", index=False)
    pc.pivot_table(index="patientId", columns="clinicalAttributeId", values="value",
                   aggfunc="first").reset_index().to_csv(f"{RAW}/{study}_clinical_patient.tsv", sep="\t", index=False)
    sl = pd.DataFrame(_get(f"/studies/{study}/sample-lists", {"projection": "DETAILED"}))
    seq = sl[sl.sampleListId == f"{study}_sequenced"].iloc[0]
    pd.Series(seq.sampleIds, name="sampleId").to_csv(f"{RAW}/{study}_sequenced_samples.tsv", sep="\t", index=False)
    mu = pd.DataFrame(_post(f"/molecular-profiles/{study}_mutations/mutations/fetch",
                            {"entrezGeneIds": list(GENES.values()), "sampleListId": f"{study}_sequenced"},
                            {"projection": "DETAILED"}))
    mu["hugo"] = mu["gene"].apply(lambda g: g["hugoGeneSymbol"])
    mu.drop(columns=["gene"]).to_csv(f"{RAW}/{study}_mutations_targetgenes.tsv", sep="\t", index=False)
    if study == "aml_ohsu_2022":
        meta = _get("/generic-assay-meta/aml_ohsu_2022_BLOOD_MEASUREMENTS")
        ids = [m["stableId"] for m in meta]
        gd = pd.DataFrame(_post("/generic_assay_data/aml_ohsu_2022_BLOOD_MEASUREMENTS/fetch",
                                {"genericAssayStableIds": ids, "sampleListId": "aml_ohsu_2022_all"}))
        gd.pivot_table(index=["sampleId", "patientId"], columns="stableId", values="value",
                       aggfunc="first").reset_index().to_csv(f"{RAW}/aml_ohsu_2022_blood_measurements.tsv",
                                                             sep="\t", index=False)


def fetch_amlsg():
    d = f"{RAW}/amlsg"; os.makedirs(d, exist_ok=True)
    base = f"https://raw.githubusercontent.com/gerstung-lab/AML-multistage/{os.environ.get('AMLSG_REF', 'a3f96ba727bfa3f30e47d9419a9faf1caab44d51')}/data/"
    for fn in ["AMLSG_Clinical_Anon.RData", "AMLSG_Genetic.txt", "README.md"]:
        if not os.path.exists(f"{d}/{fn}"):
            open(f"{d}/{fn}", "wb").write(urllib.request.urlopen(base + fn, timeout=300).read())


def fetch_flt3_reference():
    f = f"{RAW}/FLT3_P36888.fasta"
    if not os.path.exists(f):
        open(f, "w").write(urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/P36888.fasta", timeout=60).read().decode())
    return "".join(open(f).read().split("\n")[1:])


# --------------------------------------------------------------------------- FLT3 classifier
def make_flt3_classifier(flt3_seq):
    win = flt3_seq[559:650]                       # FLT3 residues 560-650 (JM + TKD1 start)

    def dup_evidence(pc):
        if pc.endswith("dup"):
            return True
        m = re.search(r"ins([A-Z\*]+)$", pc)
        ins = m.group(1) if m else ""
        for frag in ins.split("*"):
            for L in range(len(frag), 4, -1):     # >=5 consecutive native residues
                if any(frag[i:i + L] in win for i in range(len(frag) - L + 1)):
                    return True
        return False

    def classify(pc, mtype, pstart):
        pc = str(pc)
        if (("ins" in pc) or ("dup" in pc) or mtype == "In_Frame_Ins") and not pc.startswith("I836"):
            m = re.search(r"(\d+)", pc)
            pos = int(m.group(1)) if m else pstart
            if 570 <= pos <= 650:
                return "ITD_verified" if dup_evidence(pc) else "ITD_probable"
        if mtype in ("Splice_Site", "Splice_Region") and re.match(r"X(5[7-9]\d|6[0-4]\d)_splice", pc):
            return "ITD_JM_splice"
        if re.match(r"^D835[A-Z]$", pc) or pc.startswith("I836"):
            return "TKD_D835_I836"
        if re.match(r"^(D839|N841|Y842|S840|R834|M837|H838)", pc):
            return "TKD_loop_other"
        if mtype == "In_Frame_Del" and pstart is not None and 570 <= pstart <= 650:
            return "JM_del_unclassified"
        return "FLT3_other"
    return classify


def restricted_variant(hugo, mtype, pos):
    """Sensitivity definition: truncating, or missense/in-frame within catalytic domains."""
    if mtype in TRUNCATING:
        return True
    if pos is None or np.isnan(pos):
        return False
    if hugo == "DNMT3A":
        return 634 <= pos <= 912
    if hugo == "TET2":
        return (1104 <= pos <= 1481) or (1843 <= pos <= 2002)
    return False


def gene_flags(mu, samples, gene):
    x = mu[(mu.hugo == gene) & (mu.mutationType.isin(PROTEIN_ALTERING))].copy()
    x["vaf"] = x.tumorAltCount / (x.tumorAltCount + x.tumorRefCount)
    x["restr"] = [restricted_variant(gene, t, p) for t, p in zip(x.mutationType, x.proteinPosStart)]
    s = pd.Series(samples)
    any_ = s.isin(set(x.sampleId)).values
    restr = s.isin(set(x.loc[x.restr, "sampleId"])).values
    vaf = s.map(x.groupby("sampleId").vaf.max()).values
    r882 = s.isin(set(x.loc[x.proteinChange.astype(str).str.startswith("R882"), "sampleId"])).values
    return any_, restr, vaf, r882


def parse_bounded(v):
    """Harmonize censored lab strings: '>X' -> X ; '<1' -> 0.5 ; '<X' -> X/2. Returns float or NaN."""
    if pd.isna(v):
        return np.nan
    t = str(v).strip()
    if t.startswith(">"):
        return float(t[1:])
    if t.startswith("<"):
        return float(t[1:]) / 2.0
    return pd.to_numeric(t, errors="coerce")


NK_RE = re.compile(r"^\s*46\s*,\s*X[XY]\s*\[\s*\d+\s*\]\s*$")


def normal_karyotype(s):
    if pd.isna(s) or str(s).strip() == "":
        return np.nan
    return float(bool(NK_RE.match(str(s))))


# --------------------------------------------------------------------------- cohort builders
def build_tcga(classify, flow):
    st = "laml_tcga_pub"
    S = pd.read_csv(f"{RAW}/{st}_clinical_sample.tsv", sep="\t")
    P = pd.read_csv(f"{RAW}/{st}_clinical_patient.tsv", sep="\t")
    M = pd.read_csv(f"{RAW}/{st}_mutations_targetgenes.tsv", sep="\t", low_memory=False)
    seq = set(pd.read_csv(f"{RAW}/{st}_sequenced_samples.tsv", sep="\t").sampleId)
    require(P, ["AGE", "SEX", "WBC", "BM_BLAST_PERCENTAGE", "PB_BLAST_PERCENTAGE", "OS_MONTHS", "OS_STATUS",
                "CYTOGENETICS", "RISK_CYTO", "TRANSPLANT_TYPE"], st)
    require(M, ["sampleId", "hugo", "proteinChange", "mutationType", "proteinPosStart", "tumorAltCount", "tumorRefCount"], st)
    d = S.merge(P, on="patientId")
    if d.patientId.duplicated().any():
        raise DataError("TCGA: >1 sample per patient (unexpected)")
    flow.append(("TCGA", "all patients in study", len(d)))
    d = d[d.sampleId.isin(seq)]; flow.append(("TCGA", "sequenced", len(d)))
    d = d.assign(age=pd.to_numeric(d.AGE, errors="coerce"))
    d = d[d.age >= 18]; flow.append(("TCGA", "age >= 18 y", len(d)))
    f = M[M.hugo == "FLT3"].copy()
    f["cls"] = [classify(a, b, c) for a, b, c in zip(f.proteinChange, f.mutationType, f.proteinPosStart)]
    has = lambda cl: d.sampleId.isin(set(f.loc[f.cls.isin(cl), "sampleId"])).values
    out = pd.DataFrame({
        "cohort": "TCGA", "patient_id": d.patientId.values, "sample_id": d.sampleId.values, "age": d.age.values,
        "male": (d.SEX == "Male").astype(float).values,
        "npm1": d.sampleId.isin(set(M[(M.hugo == "NPM1") & M.mutationType.isin(PROTEIN_ALTERING)].sampleId)).values,
        "itd": has(["ITD_verified", "ITD_probable", "ITD_JM_splice"]),
        "itd_strict": has(["ITD_verified"]),
        "tkd": has(["TKD_D835_I836"]), "tkd_ext": has(["TKD_D835_I836", "TKD_loop_other"]),
        "itd_source": "sequence (WGS/WES)", "itd_allelic_ratio": np.nan,
        "wbc": pd.to_numeric(d.WBC, errors="coerce").values,
        "bm_blast": pd.to_numeric(d.BM_BLAST_PERCENTAGE, errors="coerce").values,
        "pb_blast": pd.to_numeric(d.PB_BLAST_PERCENTAGE, errors="coerce").values,
        "nk": [normal_karyotype(x) for x in d.CYTOGENETICS],
        "cyto_adverse": np.where(d.RISK_CYTO == "N.D.", np.nan, (d.RISK_CYTO == "Poor").astype(float)),
        "de_novo": 1.0,   # TCGA NEJM 2013 enrolled de novo AML (per source publication)
        "hsct": np.where(d.TRANSPLANT_TYPE.isna(), np.nan, (~d.TRANSPLANT_TYPE.astype(str).isin(["0"])).astype(float)),
        "os_months": pd.to_numeric(d.OS_MONTHS, errors="coerce").values,
        "os_event": d.OS_STATUS.map({"1:DECEASED": 1, "0:LIVING": 0}).values,
        "efs_months": np.nan, "efs_event": np.nan,
    })
    for g in ["DNMT3A", "TET2"]:
        a, r, v, r882 = gene_flags(M, out.sample_id, g)
        out[g.lower()], out[g.lower() + "_restr"], out[g.lower() + "_vaf"] = a, r, v
        if g == "DNMT3A":
            out["dnmt3a_r882"] = r882
    return out


def build_beataml(flow):
    st = "aml_ohsu_2022"
    S = pd.read_csv(f"{RAW}/{st}_clinical_sample.tsv", sep="\t", low_memory=False)
    P = pd.read_csv(f"{RAW}/{st}_clinical_patient.tsv", sep="\t")
    M = pd.read_csv(f"{RAW}/{st}_mutations_targetgenes.tsv", sep="\t", low_memory=False)
    B = pd.read_csv(f"{RAW}/{st}_blood_measurements.tsv", sep="\t").set_index("sampleId")
    seq = set(pd.read_csv(f"{RAW}/{st}_sequenced_samples.tsv", sep="\t").sampleId)
    require(S, ["CANCER_TYPE", "DISEASE_STAGE_AT_SPECIMEN_COLLECTION", "AGE_AT_DIAGNOSIS", "FLT3_ITD", "NPM1",
                "SPECIMEN_TYPE", "TIME_OF_SAMPLE_COLLECTION_RELATIVE_TO_INCLUSION", "KARYOTYPE", "IS_DENOVO",
                "ALLELIC_RATIO", "COHORT"], st)
    require(P, ["OS_MONTHS", "OS_STATUS", "SEX"], st)
    require(B, ["WBC_COUNT", "PERC_BLASTS_IN_BM", "PERC_BLASTS_IN_PB"], st)
    d = S.copy(); flow.append(("BeatAML", "all samples", len(d)))
    d = d[d.CANCER_TYPE == "Acute Myeloid Leukemia"]; flow.append(("BeatAML", "AML samples", len(d)))
    d = d[d.DISEASE_STAGE_AT_SPECIMEN_COLLECTION == "Initial Diagnosis"]; flow.append(("BeatAML", "initial-diagnosis specimens", len(d)))
    d = d[d.sampleId.isin(seq)]; flow.append(("BeatAML", "sequenced (WES)", len(d)))
    d = d.assign(age=pd.to_numeric(d.AGE_AT_DIAGNOSIS, errors="coerce"))
    d = d[d.age >= 18]; flow.append(("BeatAML", "age >= 18 y", len(d)))
    d = d.assign(t=pd.to_numeric(d.TIME_OF_SAMPLE_COLLECTION_RELATIVE_TO_INCLUSION, errors="coerce"),
                 site=d.SPECIMEN_TYPE.map({"Bone Marrow Aspirate": 0, "Peripheral Blood": 1, "Leukapheresis": 2}))
    d = d.sort_values(["patientId", "t", "site", "sampleId"]).drop_duplicates("patientId")
    flow.append(("BeatAML", "one diagnostic sample per patient", len(d)))
    pm = P.set_index("patientId")
    f = M[M.hugo == "FLT3"]
    tkd = lambda pat: d.sampleId.isin(set(f.loc[f.proteinChange.astype(str).str.match(pat), "sampleId"])).values
    out = pd.DataFrame({
        "cohort": "BeatAML", "patient_id": d.patientId.values, "sample_id": d.sampleId.values, "age": d.age.values,
        "male": d.patientId.map(pm.SEX).map({"Male": 1.0, "Female": 0.0}).values,
        "npm1": (d.NPM1 == "Positive").values,
        "itd": (d.FLT3_ITD == "Positive").values, "itd_known": d.FLT3_ITD.isin(["Positive", "Negative"]).values,
        "itd_strict": (d.FLT3_ITD == "Positive").values,
        "tkd": tkd(r"^(D835[A-Z]$|I836)"), "tkd_ext": tkd(r"^(D835[A-Z]$|I836|D839|N841|Y842|S840|R834|M837|H838)"),
        "itd_source": "clinical PCR/CE",
        "itd_allelic_ratio": pd.to_numeric(d.ALLELIC_RATIO, errors="coerce").values,
        "wbc": pd.to_numeric(d.sampleId.map(B.WBC_COUNT), errors="coerce").values, "bm_blast": d.sampleId.map(B.PERC_BLASTS_IN_BM).map(parse_bounded).values,
        "bm_blast_bounded": d.sampleId.map(B.PERC_BLASTS_IN_BM).astype(str).str.match(r"^[<>]").values,
        "pb_blast": pd.to_numeric(d.sampleId.map(B.PERC_BLASTS_IN_PB), errors="coerce").values,
        "nk": [normal_karyotype(x) for x in d.KARYOTYPE], "cyto_adverse": np.nan,
        "de_novo": d.IS_DENOVO.map({"TRUE": 1.0, "FALSE": 0.0}).values, "hsct": np.nan,
        "os_months": pd.to_numeric(d.patientId.map(pm.OS_MONTHS), errors="coerce").values,
        "os_event": d.patientId.map(pm.OS_STATUS).map({"1:DECEASED": 1, "0:LIVING": 0}).values,
        "efs_months": np.nan, "efs_event": np.nan,
    })
    for g in ["DNMT3A", "TET2"]:
        a, r, v, r882 = gene_flags(M, out.sample_id, g)
        out[g.lower()], out[g.lower() + "_restr"], out[g.lower() + "_vaf"] = a, r, v
        if g == "DNMT3A":
            out["dnmt3a_r882"] = r882
    # NB: the FLT3-ITD placeholder record's VAF is not interpretable -> no FLT3 VAF exported for BeatAML
    out = out[out.itd_known].drop(columns="itd_known")
    flow.append(("BeatAML", "clinical FLT3-ITD status known", len(out)))
    return out


def build_amlsg(flow):
    import pyreadr
    cl = pyreadr.read_r(f"{RAW}/amlsg/AMLSG_Clinical_Anon.RData")["clinicalData"]
    g = pd.read_csv(f"{RAW}/amlsg/AMLSG_Genetic.txt", sep="\t")
    require(cl, ["PDID", "AOD", "gender", "NPM1", "FLT3_ITD", "FLT3_TKD", "OS", "Status", "efs", "EFSSTAT",
                 "wbc", "BM_Blasts", "PB_Blasts", "NK", "C_Risk", "TypeAML", "TPL_o"], "AMLSG clinical")
    require(g, ["SAMPLE_NAME", "GENE", "AA_CHANGE", "CONSEQUENCE", "Result", "VAF"], "AMLSG genetic")
    flow.append(("AMLSG", "all patients", len(cl)))
    cl = cl[cl.AOD >= 18]; flow.append(("AMLSG", "age >= 18 y", len(cl)))
    key = lambda s: s.astype(str).str.replace(r"[a-z]$", "", regex=True)
    cl = cl.assign(key=key(cl.PDID)); g = g.assign(key=key(g.SAMPLE_NAME))
    if cl.key.duplicated().any():
        raise DataError("AMLSG duplicated patient keys")
    out = pd.DataFrame({
        "cohort": "AMLSG", "patient_id": cl.PDID.astype(str).values, "sample_id": cl.PDID.astype(str).values,
        "age": cl.AOD.values, "male": cl.gender.map({1: 1.0, 2: 0.0}).values,
        "npm1": (cl.NPM1 == 1).values, "itd": (cl.FLT3_ITD == 1).values, "itd_strict": (cl.FLT3_ITD == 1).values,
        "tkd": (cl.FLT3_TKD == 1).values, "tkd_ext": (cl.FLT3_TKD == 1).values, "itd_source": "clinical",
        "itd_allelic_ratio": np.nan,
        "wbc": cl.wbc.values, "bm_blast": cl.BM_Blasts.values, "pb_blast": cl.PB_Blasts.values,
        "nk": cl.NK.astype(float).values,
        "cyto_adverse": np.where(cl.C_Risk.isna(), np.nan, (cl.C_Risk == "Adverse").astype(float)),
        "de_novo": (cl.TypeAML == "AML").astype(float).values, "hsct": cl.TPL_o.astype(float).values,
        "os_months": (cl.OS / 30.4375).values, "os_event": cl.Status.values,
        "efs_months": (cl.efs / 30.4375).values, "efs_event": cl.EFSSTAT.values,
    })
    keys = pd.Series(cl.key.values)
    for gene in ["DNMT3A", "TET2"]:
        x = g[g.GENE == gene].copy()
        out[gene.lower()] = keys.isin(set(x.key)).values
        out[gene.lower() + "_restr"] = keys.isin(set(x.loc[x.Result == "ONCOGENIC", "key"])).values
        x["v"] = pd.to_numeric(x.VAF, errors="coerce") / 100.0
        out[gene.lower() + "_vaf"] = keys.map(x.groupby("key").v.max()).values
        if gene == "DNMT3A":
            out["dnmt3a_r882"] = keys.isin(set(x.loc[x.AA_CHANGE.astype(str).str.startswith("p.R882"), "key"])).values
    return out


def finalize(df, flow):
    for c in ["npm1", "itd", "itd_strict", "tkd", "tkd_ext", "dnmt3a", "tet2", "dnmt3a_restr", "tet2_restr", "dnmt3a_r882"]:
        df[c] = df[c].astype(bool)
    df["epi"] = df.dnmt3a | df.tet2
    df["epi_restr"] = df.dnmt3a_restr | df.tet2_restr
    df["flt3_group"] = np.select([df.itd & df.tkd, df.itd, df.tkd], ["ITD+TKD", "ITD", "TKD"], "none")
    df["os_months"] = df.os_months.clip(lower=0.01)
    df["age10"] = df.age / 10.0
    for c, d in df.groupby("cohort", sort=False):
        n1 = d[d.npm1]; flow.append((c, "NPM1 mutated", len(n1)))
        for k in ["ITD", "TKD", "ITD+TKD", "none"]:
            flow.append((c, f"  FLT3 status: {k}", int((n1.flt3_group == k).sum())))
        e = n1[n1.flt3_group.isin(["ITD", "TKD"])]; flow.append((c, "ITD-only or TKD-only", len(e)))
        e = e[e.os_months.notna() & e.os_event.notna()]; flow.append((c, "OS available (primary analysis set)", len(e)))
    df["primary"] = df.npm1 & df.flt3_group.isin(["ITD", "TKD"]) & df.os_months.notna() & df.os_event.notna()
    df["group4"] = np.where(df.flt3_group.isin(["ITD", "TKD"]),
                            "NPM1/" + df.flt3_group + "/" + np.where(df.epi, "epi-mut", "epi-wt"), "")
    p = df[df.primary]
    if p.duplicated(["cohort", "patient_id"]).any():
        raise DataError("duplicate patient in primary set")
    if (p.age < 18).any():
        raise DataError("pediatric patient in primary set")
    if ((p.wbc < 0) | (p.wbc > 1000)).any() or ((p.bm_blast < 0) | (p.bm_blast > 100)).any():
        raise DataError("implausible WBC/blast value")
    return df


DICTIONARY = [
    ("cohort", "study", "-", "TCGA = laml_tcga_pub; BeatAML = aml_ohsu_2022; AMLSG = AMLSG trials", ["cohort"]),
    ("patient_id", "patientId / PDID", "-", "patient identifier", ["patient_id"]),
    ("sample_id", "sampleId / PDID", "-", "diagnostic sample analysed", ["sample_id"]),
    ("age", "AGE / AGE_AT_DIAGNOSIS / AOD", "years", "age at diagnosis", ["age"]),
    ("male", "SEX / SEX / gender", "0/1", "1 = male", ["male"]),
    ("npm1", "mutation file NPM1 / NPM1 (PCR/CE) / NPM1", "bool", "NPM1 mutated", ["npm1"]),
    ("itd", "FLT3 classifier / FLT3_ITD / FLT3_ITD", "bool", "FLT3-ITD present (TCGA: verified+probable+JM splice)", ["itd"]),
    ("itd_strict", "FLT3 classifier", "bool", "TCGA: duplication-verified ITD only", ["itd_strict"]),
    ("tkd", "mutation file / WES / FLT3_TKD", "bool", "FLT3 D835X or I836 alteration", ["tkd"]),
    ("tkd_ext", "mutation file", "bool", "TKD plus other activation-loop alterations R834/M837/H838/D839/S840/N841/Y842 (cBioPortal cohorts)", ["tkd_ext"]),
    ("flt3_group", "derived", "-", "ITD, TKD, ITD+TKD, none", ["flt3_group"]),
    ("dnmt3a, tet2", "mutation file / AMLSG_Genetic.txt", "bool", "any protein-altering variant (AMLSG: ONCOGENIC or POSSIBLE)", ["dnmt3a", "tet2"]),
    ("dnmt3a_restr, tet2_restr", "derived", "bool", "truncating or catalytic-domain variant (AMLSG: ONCOGENIC)", ["dnmt3a_restr", "tet2_restr"]),
    ("dnmt3a_r882", "proteinChange / AA_CHANGE", "bool", "DNMT3A R882 variant", ["dnmt3a_r882"]),
    ("dnmt3a_vaf, tet2_vaf", "alt/(alt+ref) / VAF", "fraction", "max VAF of the gene's variants (NA if wild-type)", ["dnmt3a_vaf", "tet2_vaf"]),
    ("epi", "derived", "bool", "DNMT3A and/or TET2 mutated", ["epi"]),
    ("wbc", "WBC / BLOOD_MEASUREMENTS:WBC_COUNT / wbc", "x10^9/L", "baseline white-cell count", ["wbc"]),
    ("bm_blast", "BM_BLAST_PERCENTAGE / PERC_BLASTS_IN_BM / BM_Blasts", "%", "bone-marrow blasts", ["bm_blast"]),
    ("pb_blast", "PB_BLAST_PERCENTAGE / PERC_BLASTS_IN_PB / PB_Blasts", "%", "peripheral-blood blasts", ["pb_blast"]),
    ("nk", "CYTOGENETICS / KARYOTYPE / NK", "0/1", "normal karyotype (46,XX or 46,XY only)", ["nk"]),
    ("cyto_adverse", "RISK_CYTO / - / C_Risk", "0/1", "adverse cytogenetic risk (not available for BeatAML)", ["cyto_adverse"]),
    ("de_novo", "source publication / IS_DENOVO / TypeAML", "0/1", "de novo AML", ["de_novo"]),
    ("hsct", "TRANSPLANT_TYPE / - / TPL_o", "0/1", "any HSCT recorded (not time-dependent; not available for BeatAML)", ["hsct"]),
    ("itd_allelic_ratio", "ALLELIC_RATIO", "ratio", "BeatAML only", ["itd_allelic_ratio"]),
    ("os_months, os_event", "OS_MONTHS,OS_STATUS / OS (days),Status", "months, 0/1", "overall survival from diagnosis", ["os_months", "os_event"]),
    ("efs_months, efs_event", "efs (days), EFSSTAT", "months, 0/1", "AMLSG only", ["efs_months", "efs_event"]),
]


def main():
    fetch_amlsg()
    for st in ["laml_tcga_pub", "aml_ohsu_2022"]:
        fetch_cbioportal(st)
    classify = make_flt3_classifier(fetch_flt3_reference())
    flow = []
    df = pd.concat([build_tcga(classify, flow), build_beataml(flow), build_amlsg(flow)], ignore_index=True)
    df = finalize(df, flow)
    df.to_csv(f"{OUT}/analytic_dataset.csv", index=False)
    fl = pd.DataFrame(flow, columns=["cohort", "step", "n"]); fl.to_csv(f"{OUT}/cohort_flow.csv", index=False)
    p = df[df.primary]
    dd = []
    for var, src, unit, defin, cols in DICTIONARY:
        miss = {c: {coh: f"{p.loc[p.cohort == coh, c].isna().mean()*100:.0f}%" for coh in ["TCGA", "BeatAML", "AMLSG"]} for c in cols}
        dd.append(dict(variable=var, source_fields=src, unit=unit, definition=defin,
                       missing_in_primary_set_by_cohort=json.dumps(miss)))
    pd.DataFrame(dd).to_csv(f"{OUT}/data_dictionary.csv", index=False)
    return df, fl


if __name__ == "__main__":
    df, fl = main()
    print(fl.to_string())
    print(df[df.primary].groupby(["cohort", "group4"]).agg(n=("patient_id", "size"), deaths=("os_event", "sum")))
