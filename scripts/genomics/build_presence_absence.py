#!/usr/bin/env python3

import csv, sys, re, glob, os
from collections import defaultdict

ALC_REF = "BGC0002464"          # alcaligin BGC MiBIG ref
MET_STRAIN = "DSM26089"         # strain with NRP-metallophore
MET_PRODUCT = "NRP-metallophore"

def strain_of(gbk):
    b = os.path.basename(gbk)
    if b.startswith("BGC") or "MIBIG" in b.upper():
        return None
    return b.split("_")[0].split(".")[0]

def main(mix, as_dir, out_dir):
    if mix == "auto":
        cands = glob.glob(os.path.join(out_dir, "output_files", "*", "mix",
                                       "mix_clustering_c0.3.tsv"))
        mix = max(cands, key=os.path.getmtime)
    sys.stderr.write(f"mix clustering: {mix}\n")

    fam_strains = defaultdict(set)      
    fam_has_alc = set()                 
    rec_fam = {}                        
    with open(mix) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            fam = row["Family"]
            gbk = row["GBK"]
            rec_fam[os.path.basename(gbk)] = fam
            if ALC_REF in gbk:
                fam_has_alc.add(fam)
            s = strain_of(gbk)
            if s:
                fam_strains[fam].add(s)

    # find antiSMASH regions whose product is NRP-metallophore
    met_fams = set()
    sdir = os.path.join(as_dir, MET_STRAIN)
    for gbk in glob.glob(os.path.join(sdir, "*.region*.gbk")):
        txt = open(gbk, encoding="utf-8", errors="ignore").read()
        m = re.search(r'/product="([^"]+)"', txt)
        if m and m.group(1) == MET_PRODUCT:
            base = os.path.basename(gbk)[:-4]            
            fam = rec_fam.get(base) or rec_fam.get(base + ".gbk_region_1")
            if not fam:
                for k, v in rec_fam.items():
                    if k.startswith(base):
                        fam = v; break
            if fam:
                met_fams.add(fam)

    fams = [f for f in fam_strains if fam_strains[f]]
    strains = sorted({s for f in fams for s in fam_strains[f]})
    fams.sort(key=lambda f: (-len(fam_strains[f]), f))
    new_id = {f: f"FAM_{i+1:05d}" for i, f in enumerate(fams)}

    os.makedirs(out_dir, exist_ok=True)
    mat_p = os.path.join(out_dir, "absence_presence_clean.tsv")
    with open(mat_p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([""] + [new_id[f] for f in fams])
        for s in strains:
            w.writerow([s] + [1 if s in fam_strains[f] else 0 for f in fams])

    roles_p = os.path.join(out_dir, "family_roles.tsv")
    with open(roles_p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["role", "fam", "n_strains"])
        for f in sorted(fam_has_alc, key=lambda f: new_id[f]):
            w.writerow(["alc", new_id[f], len(fam_strains[f])])
        for f in sorted(met_fams, key=lambda f: new_id[f]):
            w.writerow(["nrp_metallophore", new_id[f], len(fam_strains[f])])

    sys.stderr.write(f"wrote {mat_p} ({len(strains)} strains x {len(fams)} families)\n")
    sys.stderr.write(f"alc -> {[new_id[f] for f in fam_has_alc]}; "
                     f"nrp_metallophore({MET_STRAIN}) -> {[new_id[f] for f in met_fams]}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
