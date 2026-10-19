
import os
import re

import numpy as np
import pandas as pd

INFILE = "rawData/proteomics/proteinGroups.xlsx"
ENRICH = "derivedData/proteomics/coculture_enrichment.xlsx"
OUTDIR = "derivedData/proteomics"

CO_LFQ = ["LFQ intensity C1", "LFQ intensity C5", "LFQ intensity C6"]
MONO_LFQ = ["LFQ intensity M1", "LFQ intensity M6", "LFQ intensity M7"]
CO_MSMS = ["MS/MS count C1", "MS/MS count C5", "MS/MS count C6"]
MONO_MSMS = ["MS/MS count M1", "MS/MS count M6", "MS/MS count M7"]
FC_CUTOFF = np.log2(1)
FDR_CUTOFF = 0.05
MSMS_MIN = 2

GROUPS = {
    "Pyoverdine (synthesis & uptake)": (
        r"^(?:pvd|fpv|pvc)", r"pyoverdin|ferripyoverdine"),
    "Pyochelin (synthesis & uptake)": (
        r"^(?:pch|fpt)", r"pyochelin|ferripyochelin|salicylate|dihydroaeruginoate"),
    "Heme uptake": (
        r"^(?:has|phu|hxu|hem[oO])", r"heme|haem"),
    "Ferrous / other iron transport": (
        r"^(?:feo|fec|fox|fiu|fem|hit|pfe|hasR|tonB|exb|opt|cir|hxu)",
        r"ferric|ferrous|iron[ -]|iron$|siderophore|TonB|iron-regulated|"
        r"iron uptake|iron transport"),
    "Iron storage / regulation": (
        r"^(?:fur|bfr|ftn|dps|sod[BM])", r"bacterioferritin|ferritin|"
        r"ferric uptake regulator|superoxide dismutase"),
    "Secreted proteases": (
        r"^(?:lasA|lasB|aprA|piv|prpL|mucD|ctpA|lap)",
        r"elastase|alkaline protease|protease IV|protease 4|staphylolytic|"
        r"aminopeptidase|serralysin|serine protease|zinc metalloprotease"),
    "Toxins / phospholipases / hemolysins": (
        r"^(?:toxA|eta|exo[STUY]|plc|rhl[AB]|cbpD|hcn)",
        r"exotoxin|phospholipase|hemolysin|leukocidin|cytotoxin|"
        r"chitin-binding|hydrogen cyanide"),
    "Type III / VI secretion": (
        r"^(?:psc|pcr|exs|pop|tse|tss|hcp|clp[VB]|vgr|imp)",
        r"type III secretion|type VI secretion|T3SS|T6SS|translocon"),
}


def header_fields(acc, headers):
    parts = str(headers).split(";")
    chosen = next((h for h in parts if acc in h), parts[0]).strip()
    desc = chosen.split(" ", 1)[1].split(" OS=")[0] if " " in chosen else ""
    gene = chosen.split(" GN=")[1].split()[0] if " GN=" in chosen else ""
    return gene, desc


def main():
    df = pd.read_excel(INFILE)
    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df[df[col].isna() | (df[col] != "+")]

    df["Accession"] = df["Majority protein IDs"].astype(str).str.split(";").str[0].str.strip()
    gd = [header_fields(a, h) for a, h in zip(df["Accession"], df["Fasta headers"])]
    df["Gene"] = [g for g, _ in gd]
    df["Description"] = [d for _, d in gd]

    co = np.log2(df[CO_LFQ].replace(0, np.nan))
    mono = np.log2(df[MONO_LFQ].replace(0, np.nan))
    df["nC"] = co.notna().sum(axis=1).values
    df["nM"] = mono.notna().sum(axis=1).values
    df["log2_C"] = co.mean(axis=1).values
    df["log2_M"] = mono.mean(axis=1).values
    df["log2FC"] = df["log2_C"] - df["log2_M"]
    df["msms_total"] = df["MS/MS count"]
    df["msms_C"] = df[CO_MSMS].sum(axis=1)
    df["msms_M"] = df[MONO_MSMS].sum(axis=1)
    any_lfq = (df[CO_LFQ + MONO_LFQ] > 0).any(axis=1)
    df["confident_id"] = any_lfq & (df["msms_total"] >= MSMS_MIN)

    at = pd.read_excel(ENRICH, "all_tested")[["Accession", "p_value", "adj_p_value_BH"]]
    df = df.merge(at, on="Accession", how="left")

    def category(r):
        if not r["confident_id"]:
            return "not detected" if not (r[CO_LFQ + MONO_LFQ] > 0).any() else "low-confidence"
        if r["nC"] >= 2 and r["nM"] == 0 and r["msms_C"] >= MSMS_MIN:
            return "coculture-ONLY"
        if r["nM"] >= 2 and r["nC"] == 0 and r["msms_M"] >= MSMS_MIN:
            return "monoculture-ONLY"
        if pd.notna(r["adj_p_value_BH"]):
            if r["log2FC"] >= FC_CUTOFF and r["adj_p_value_BH"] < FDR_CUTOFF:
                return "enriched"
            if r["log2FC"] <= -FC_CUTOFF and r["adj_p_value_BH"] < FDR_CUTOFF:
                return "depleted"
            return "ns"
        return "sparse"

    df["category"] = df.apply(lambda r: category(r), axis=1)

    def fmt(r):
        gene = r["Gene"] or "-"
        cat = r["category"]
        fc = f"{r['log2FC']:+5.2f}" if cat in ("enriched", "depleted", "ns") else "  --  "
        adjp = f"{r['adj_p_value_BH']:.2e}" if pd.notna(r["adj_p_value_BH"]) else "   --   "
        flag = "*" if cat in ("enriched", "coculture-ONLY") else " "
        return (f"  {flag}{gene:<9} {r['Accession']:<8} {cat:<16} log2FC={fc}  "
                f"adjP={adjp}  C[{int(r['nC'])}/3] M[{int(r['nM'])}/3] "
                f"MS/MS C{int(r['msms_C'])}/M{int(r['msms_M'])}  {str(r['Description'])[:44]}")

    order = {"enriched": 0, "coculture-ONLY": 1, "sparse": 2, "ns": 3,
             "depleted": 4, "monoculture-ONLY": 5, "low-confidence": 6, "not detected": 7}
    matched = set()
    lines, rows = [], []
    for title, (gre, dre) in GROUPS.items():
        m = df[(df["Gene"].str.contains(gre, case=False, regex=True, na=False) |
                df["Description"].str.contains(dre, case=False, regex=True, na=False))
               & ~df.index.isin(matched)]
        matched |= set(m.index)
        if len(m) == 0:
            continue
        m = m.assign(_o=m["category"].map(order)).sort_values(
            ["_o", "log2FC"], ascending=[True, False], na_position="last")
        lines.append(f"\n### {title}  (n={len(m)})")
        for _, r in m.iterrows():
            lines.append(fmt(r))
            rr = r.copy(); rr["functional_group"] = title; rows.append(rr)

    print(f"Legend: * = enriched in coculture. Cutoffs: 1.5-fold (log2FC>={FC_CUTOFF:.3f}), "
          f"BH p<{FDR_CUTOFF}; confidence LFQ>0 & MS/MS>={MSMS_MIN}.")
    print("\n".join(lines))

    out = pd.DataFrame(rows)
    cols = ["functional_group", "Accession", "Gene", "Description", "category",
            "log2FC", "p_value", "adj_p_value_BH", "log2_C", "log2_M",
            "nC", "nM", "msms_C", "msms_M", "msms_total", "confident_id"]
    outfile = os.path.join(OUTDIR, "coculture_iron_virulence.tsv")
    out[cols].to_csv(outfile, sep="\t", index=False)
    print(f"\nwrote {outfile}  ({len(out)} proteins)")

    print("\nNamed targets:")
    for label, pat in [("Protease IV (PrpL/piv)", r"protease IV|prpL|\bpiv\b|lysyl endopeptidase"),
                       ("Exotoxin A (toxA/eta)", r"exotoxin|toxA|\beta\b")]:
        hit = df[df["Description"].str.contains(pat, case=False, regex=True, na=False) |
                 df["Gene"].str.contains(pat, case=False, regex=True, na=False)]
        for _, r in hit.iterrows():
            print("  " + fmt(r).strip())


if __name__ == "__main__":
    main()
