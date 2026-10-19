import os

import numpy as np
import pandas as pd
from scipy import stats

INFILE = "rawData/proteomics/proteinGroups.xlsx"
OUTDIR = "derivedData/proteomics"
OUTFILE = os.path.join(OUTDIR, "coculture_enrichment.xlsx")

CO_LFQ = ["LFQ intensity C1", "LFQ intensity C5", "LFQ intensity C6"]
MONO_LFQ = ["LFQ intensity M1", "LFQ intensity M6", "LFQ intensity M7"]
CO_MSMS = ["MS/MS count C1", "MS/MS count C5", "MS/MS count C6"]
MONO_MSMS = ["MS/MS count M1", "MS/MS count M6", "MS/MS count M7"]

FC_CUTOFF = np.log2(1) 
FDR_CUTOFF = 0.05
MSMS_MIN = 2               # confident ID: total MS/MS count >= 2


def bh_fdr(pvals):
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    return adj


def header_fields(acc, headers):
    """Return (gene, description) from the fasta header matching accession `acc`."""
    parts = str(headers).split(";")
    chosen = next((h for h in parts if acc in h), parts[0])
    chosen = chosen.strip()
    desc = chosen.split(" ", 1)[1].split(" OS=")[0] if " " in chosen else ""
    gene = chosen.split(" GN=")[1].split()[0] if " GN=" in chosen else ""
    return gene, desc


def main():
    df = pd.read_excel(INFILE)

    # MaxQuant QC flags
    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df[df[col].isna() | (df[col] != "+")]

    # representative protein = leading Majority protein ID
    df["Accession"] = df["Majority protein IDs"].astype(str).str.split(";").str[0].str.strip()
    gd = [header_fields(a, h) for a, h in zip(df["Accession"], df["Fasta headers"])]
    df["Gene"] = [g for g, _ in gd]
    df["Description"] = [d for _, d in gd]

    # MS/MS support
    df["msms_total"] = df["MS/MS count"]
    df["msms_C"] = df[CO_MSMS].sum(axis=1)
    df["msms_M"] = df[MONO_MSMS].sum(axis=1)

    # log2 LFQ, zeros -> NaN
    co = np.log2(df[CO_LFQ].replace(0, np.nan))
    mono = np.log2(df[MONO_LFQ].replace(0, np.nan))
    df["n_valid_coculture"] = co.notna().sum(axis=1).values
    df["n_valid_monoculture"] = mono.notna().sum(axis=1).values
    df["log2_mean_coculture"] = co.mean(axis=1).values
    df["log2_mean_monoculture"] = mono.mean(axis=1).values

    any_lfq = (df[CO_LFQ + MONO_LFQ] > 0).any(axis=1)
    df["confident_id"] = any_lfq & (df["msms_total"] >= MSMS_MIN)

    n_before = len(df)
    df = df[df["confident_id"]].copy()
    co = np.log2(df[CO_LFQ].replace(0, np.nan))
    mono = np.log2(df[MONO_LFQ].replace(0, np.nan))
    n_co = df["n_valid_coculture"]
    n_mono = df["n_valid_monoculture"]

    # quantitative comparison: >=2 valid in both groups
    testable = (n_co >= 2) & (n_mono >= 2)
    log2fc = df["log2_mean_coculture"] - df["log2_mean_monoculture"]
    pvals = np.full(len(df), np.nan)
    idx = np.where(testable.values)[0]
    pv = [
        stats.ttest_ind(c[~np.isnan(c)], m[~np.isnan(m)], equal_var=False).pvalue
        for c, m in zip(co.values[idx], mono.values[idx])
    ]
    pvals[idx] = pv
    adj = np.full(len(df), np.nan)
    adj[idx] = bh_fdr(pv)

    df["log2FC_co_vs_mono"] = log2fc.values
    df["p_value"] = pvals
    df["adj_p_value_BH"] = adj

    def categorize(r):
        if r["n_valid_coculture"] >= 2 and r["n_valid_monoculture"] == 0 and r["msms_C"] >= MSMS_MIN:
            return "coculture-only"
        if r["n_valid_monoculture"] >= 2 and r["n_valid_coculture"] == 0 and r["msms_M"] >= MSMS_MIN:
            return "monoculture-only"
        if pd.notna(r["adj_p_value_BH"]):
            if r["log2FC_co_vs_mono"] >= FC_CUTOFF and r["adj_p_value_BH"] < FDR_CUTOFF:
                return "enriched"
            if r["log2FC_co_vs_mono"] <= -FC_CUTOFF and r["adj_p_value_BH"] < FDR_CUTOFF:
                return "depleted"
            return "ns"
        return "sparse"

    df["category"] = df.apply(categorize, axis=1)

    keep = [
        "Accession", "Gene", "Description", "category",
        "log2FC_co_vs_mono", "p_value", "adj_p_value_BH",
        "log2_mean_coculture", "log2_mean_monoculture",
        "n_valid_coculture", "n_valid_monoculture",
        "msms_total", "msms_C", "msms_M",
        "Razor + unique peptides", "Unique peptides", "Score",
    ]
    keep_pa = [c for c in keep if c not in ("log2FC_co_vs_mono", "p_value", "adj_p_value_BH")]

    enriched = df[df["category"] == "enriched"][keep].sort_values(
        "log2FC_co_vs_mono", ascending=False)
    co_only = df[df["category"] == "coculture-only"][keep_pa].sort_values(
        "log2_mean_coculture", ascending=False)
    depleted = df[df["category"] == "depleted"][keep].sort_values("log2FC_co_vs_mono")
    mono_only = df[df["category"] == "monoculture-only"][keep_pa].sort_values(
        "log2_mean_monoculture", ascending=False)
    all_tested = df[testable.values][keep].sort_values("p_value")
    all_proteins = df[keep].sort_values(
        ["category", "log2FC_co_vs_mono"], ascending=[True, False], na_position="last")

    os.makedirs(OUTDIR, exist_ok=True)
    with pd.ExcelWriter(OUTFILE, engine="openpyxl") as xl:
        enriched.to_excel(xl, sheet_name="enriched_in_coculture", index=False)
        co_only.to_excel(xl, sheet_name="coculture_only", index=False)
        depleted.to_excel(xl, sheet_name="depleted_in_coculture", index=False)
        mono_only.to_excel(xl, sheet_name="monoculture_only", index=False)
        all_tested.to_excel(xl, sheet_name="all_tested", index=False)
        all_proteins.to_excel(xl, sheet_name="all_proteins", index=False)

    print(f"proteins after MaxQuant QC flags:       {n_before}")
    print(f"passed confidence filter (LFQ>0 & MS/MS>={MSMS_MIN}): {len(df)}  "
          f"(removed {n_before - len(df)})")
    print(f"testable (>=2 valid LFQ per group):     {int(testable.sum())}")
    print(f"ENRICHED in coculture (>=1.5x, BH p<{FDR_CUTOFF}): {len(enriched)}")
    print(f"coculture-only (>=2 valid C, 0 in M):   {len(co_only)}")
    print(f"depleted in coculture:                  {len(depleted)}")
    print(f"monoculture-only:                       {len(mono_only)}")
    print(f"\nwrote {OUTFILE}")

    print("\nTop enriched (by log2FC):")
    for _, r in enriched.head(30).iterrows():
        g = r["Gene"] if isinstance(r["Gene"], str) and r["Gene"] else "-"
        print(f"  {r['Accession']:<8} {g:<9} log2FC={r['log2FC_co_vs_mono']:+.2f} "
              f"({2**r['log2FC_co_vs_mono']:.1f}x) adjP={r['adj_p_value_BH']:.3g}  "
              f"{str(r['Description'])[:52]}")


if __name__ == "__main__":
    main()
