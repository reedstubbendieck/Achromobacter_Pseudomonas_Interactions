import os
import re
import warnings

import numpy as np
import pandas as pd
from scipy import stats

INFILE = "rawData/proteomics/Combined proteinGroups.xlsx"   # 6 coculture + 6 monoculture
LOCFILE = "references/PAO1/subcellular_localization.csv"  # curated pseudomonas.com
OUTFILE = "derivedData/proteomics/volcano_plot_data.tsv"

CO_LFQ = [f"LFQ intensity C{i}" for i in range(1, 7)]   # C1..C6
MONO_LFQ = [f"LFQ intensity M{i}" for i in range(1, 7)]  # M1..M6
LFQ = CO_LFQ + MONO_LFQ

FC_CUTOFF = 1.0 
FDR_CUTOFF = 0.05
MSMS_MIN = 2


VF_FILE = "references/PAO1/virulenceFactors.csv"
PW_FILE = "references/PAO1/pathways.txt"
CLASSIFIER = "PseudoCAP" 

SECMET_VF = {"Phenazine biosynthesis", "Rhamnolipid biosynthesis",
             "Hydrogen cyanide production"}
SECMET_PW = {"Phenazine biosynthesis"}

IRON_VF = {"Pyoverdine", "Pyochelin", "Pyochelin receptor"}
IRON_PW = {"Pyoverdine synthesis", "Pyochelin synthesis"}
IRON_GENE = (r"^(?:pvd|fpv|pvc|pch|fpt|has[AR]|phu|hxu|feo|fec|fox|fiu|fem|"
             r"hit|pfe|tonB|exb|fur|bfr|ftn)")
IRON_DESC = (r"pyoverdin|pyochelin|ferripyoverdine|ferripyochelin|siderophore|"
             r"heme acquisition|heme/hemoglobin|heme-transport|heme uptake|"
             r"hemoglobin uptake|ferric|ferrous|iron uptake|iron transport|"
             r"iron-regulated|tonb-dependent|bacterioferritin|ferritin|"
             r"ferric uptake regul")


def header_fields(acc, headers):
    parts = str(headers).split(";")
    chosen = next((h for h in parts if acc in h), parts[0]).strip()
    desc = chosen.split(" ", 1)[1].split(" OS=")[0] if " " in chosen else ""
    gene = chosen.split(" GN=")[1].split()[0] if " GN=" in chosen else ""
    return gene, desc


def load_curated():
    """VFDB-only locus/gene -> Virulence Category, and locus -> set(Pathway Names)."""
    vf = pd.read_csv(VF_FILE)
    vf.columns = [c.strip() for c in vf.columns]
    vf = vf[vf["Source"].astype(str).str.strip().str.upper() == "VFDB"]   # VFDB only
    vf["Locus Tag"] = vf["Locus Tag"].astype(str).str.strip()
    vf["Gene Name"] = vf["Gene Name"].fillna("").astype(str).str.strip()
    loc2cat, gene2cat = {}, {}
    for lt, gn, cat in zip(vf["Locus Tag"], vf["Gene Name"], vf["Virulence Category"]):
        loc2cat.setdefault(lt, cat)
        if gn:
            gene2cat.setdefault(gn, cat)
    pw = pd.read_csv(PW_FILE, sep="\t")
    pw = pw[pw["Pathway Database"].astype(str).str.strip() == CLASSIFIER]
    pw["Locus Tag"] = pw["Locus Tag"].astype(str).str.strip()
    loc2pw = {}
    for lt, name in zip(pw["Locus Tag"], pw["Pathway Name"]):
        loc2pw.setdefault(lt, set()).add(str(name))
    return loc2cat, gene2cat, loc2pw


def highlight_curated(gene, desc, keys, loc2cat, gene2cat, loc2pw):
    """Iron: broad regex + VFDB siderophores. Sec Met / Virulence: VFDB-only curated."""
    vfcat = None
    for k in keys:
        if k in loc2cat:
            vfcat = loc2cat[k]; break
        if k in gene2cat:
            vfcat = gene2cat[k]; break
    pws = set().union(*(loc2pw.get(k, set()) for k in keys)) if keys else set()
    iron_regex = bool(re.search(IRON_GENE, str(gene), re.I) or
                      re.search(IRON_DESC, str(desc), re.I))
    if iron_regex or vfcat in IRON_VF or (pws & IRON_PW):
        return "Iron"
    if vfcat in SECMET_VF or (pws & SECMET_PW):
        return "Secondary Metabolism"
    if vfcat is not None:
        return "Virulence"
    return "Not highlighted"


def load_localization():
    """Curated pseudomonas.com subcellular localization -> Accession -> localization string.
    One row per accession, preferring an Extracellular call when a gene has several."""
    lz = pd.read_csv(LOCFILE, dtype=str)
    lz.columns = [c.strip().lstrip("﻿") for c in lz.columns]
    acc_c = next(c for c in lz.columns if re.search(r"uniprot.*accession", c, re.I))
    loc_c = next(c for c in lz.columns if re.search(r"local|subcell", c, re.I))
    lz = lz[[acc_c, loc_c]].rename(columns={acc_c: "Accession", loc_c: "localization"})
    lz["Accession"] = lz["Accession"].astype(str).str.strip()
    lz["localization"] = lz["localization"].astype(str).str.strip().str.lower()
    lz = lz[lz["Accession"].ne("") & lz["Accession"].ne("nan")]
    lz["is_ex"] = lz["localization"].eq("extracellular")
    return lz.sort_values("is_ex", ascending=False).drop_duplicates("Accession", keep="first")


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    adj = np.empty(n); adj[o] = np.minimum(r, 1.0)
    return adj


def main():
    df = pd.read_excel(INFILE)
    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df[df[col].isna() | (df[col] != "+")]
    df["Accession"] = df["Majority protein IDs"].astype(str).str.split(";").str[0].str.strip()
    gd = [header_fields(a, h) for a, h in zip(df["Accession"], df["Fasta headers"])]
    df["Gene"] = [g for g, _ in gd]
    df["Description"] = [d for _, d in gd]

    any_lfq = (df[LFQ] > 0).any(axis=1)
    df = df[any_lfq & (df["MS/MS count"] >= MSMS_MIN)].copy().reset_index(drop=True)

    # detection-limit (MinDet) imputation: zeros -> 1st-percentile of each sample's
    # measured log2 LFQ (left-censoring at the limit of detection).
    l2 = np.log2(df[LFQ].astype(float).replace(0, np.nan))
    floors = l2.quantile(0.01)
    l2 = l2.apply(lambda col: col.fillna(floors[col.name]))
    idx = np.array([int(re.search(r"(\d+)$", c).group(1)) for c in LFQ])
    cond = np.array([1.0 if c in CO_LFQ else 0.0 for c in LFQ])   # 1 = coculture
    batch = (idx >= 4).astype(float)                              # 0 = run1, 1 = run2
    Xd = np.column_stack([np.ones(len(LFQ)), cond, batch])
    Y = l2[LFQ].values                            # proteins x samples
    XtXi = np.linalg.inv(Xd.T @ Xd)
    Beta = XtXi @ Xd.T @ Y.T                       # [intercept, condition, batch] x proteins
    resid = Y.T - Xd @ Beta
    dfree = len(LFQ) - Xd.shape[1]
    sigma2 = (resid ** 2).sum(0) / dfree
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = Beta[1] / np.sqrt(sigma2 * XtXi[1, 1])
    pv = np.nan_to_num(2 * stats.t.sf(np.abs(tstat), dfree), nan=1.0)
    df["log2FC"] = Beta[1]                          # batch-adjusted coculture - monoculture
    df["p_value"] = pv
    df["adj_p"] = bh_fdr(pv)
    df["neg_log10_adj_p"] = -np.log10(np.clip(df["adj_p"], 1e-300, None))

    df["nC"] = (df[CO_LFQ] > 0).sum(axis=1)
    df["nM"] = (df[MONO_LFQ] > 0).sum(axis=1)
    # overall protein intensity (MaxQuant total Intensity), log10, for point size
    df["log10_intensity"] = np.log10(df["Intensity"].astype(float).clip(lower=1))
    df["label"] = np.where(df["Gene"].astype(str).str.len() > 0, df["Gene"], df["Accession"])
    df["label"] = df["label"].replace({"eta": "toxA", "lap": "paaP",   # PaaP = PA2939
                                        "PA3866": "Pyocin S4",
                                        "PA4625": "cdrA"})             # CdrA (pseudomonas.com)
    loc2cat, gene2cat, loc2pw = load_curated()
    df["highlight"] = [highlight_curated(g, d, [str(g), str(l)], loc2cat, gene2cat, loc2pw)
                       for g, d, l in zip(df["Gene"], df["Description"], df["label"])]
    # curator-directed category overrides (collaborator: these are virulence factors)
    HIGHLIGHT_OVERRIDE = {"paaP": "Virulence",   # PaaP / PA2939 (aminopeptidase)
                          "plcB": "Virulence"}   # phospholipase C
    for lb, cat in HIGHLIGHT_OVERRIDE.items():
        df.loc[df["label"] == lb, "highlight"] = cat
    df["significant"] = (df["adj_p"] < FDR_CUTOFF) & (df["log2FC"].abs() >= FC_CUTOFF)
    df["do_label"] = (df["highlight"] != "Not highlighted") & df["significant"]

    lz = load_localization()
    df = df.merge(lz[["Accession", "localization"]], on="Accession", how="left")
    df["curated_extracellular"] = df["localization"].fillna("").eq("extracellular")
    df["secreted"] = df["curated_extracellular"]

    df.loc[df["Accession"] == "Q03023", "secreted"] = True

    keep = ["Accession", "Gene", "Description", "log2FC", "p_value", "adj_p",
            "neg_log10_adj_p", "significant", "highlight", "label", "do_label",
            "nC", "nM", "log10_intensity", "secreted", "localization"]
    out = df[keep].copy()
    for b in ["significant", "do_label", "secreted"]:
        out[b] = out[b].astype(int)
    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
    out.to_csv(OUTFILE, sep="\t", index=False)

    n_up = int((df["significant"] & (df["log2FC"] > 0)).sum())
    n_dn = int((df["significant"] & (df["log2FC"] < 0)).sum())
    print(f"confident proteins tested (MinDet detection-limit): {len(df)}")
    print(f"significant: {int(df['significant'].sum())}  (up {n_up}, down {n_dn})")
    print(f"log2FC range: {df['log2FC'].min():.1f} .. {df['log2FC'].max():.1f}")
    print(f"one-condition proteins (nC==0 or nM==0): {int(((df.nC==0)|(df.nM==0)).sum())}")
    print(f"secreted = curated EXTRACELLULAR only: {int(df['secreted'].sum())}")
    print(f"curated localization coverage: {int(df['localization'].notna().sum())}/{len(df)}")
    print("\nhighlighted & significant:")
    q = df[(df.highlight != 'Not highlighted') & df.significant]
    for _, r in q.sort_values("log2FC", ascending=False).iterrows():
        print(f"  {r['label']:<9} {r['highlight']:<22} log2FC={r['log2FC']:+6.2f} "
              f"adjP={r['adj_p']:.2g}  C{int(r['nC'])}/M{int(r['nM'])}")
    print(f"\nwrote {OUTFILE}")


if __name__ == "__main__":
    main()
