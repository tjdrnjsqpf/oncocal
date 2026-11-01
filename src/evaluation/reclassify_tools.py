"""2-way supervised/unsupervised reclassification of tools + recompute driver under-scoring rate (R3).

Principle (user-defined):
  supervised        = discriminative model trained on clinical/pathogenicity labels (ClinVar/HGMD/disease variants) or proxy discriminative labels
  unsupervised/self = scores derived without clinical labels, from conservation, evolution, generative language models, or population genetics
Borderline cases are decided on the actual training method (see the RATIONALE log below). No guessing — grounded in the training method per the literature.
Output: add type2 to driver_miss_rate.tsv, recompute the two group means.
Usage: python -m src.evaluation.reclassify_tools
"""
from __future__ import annotations
import os
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
from src.paths import table_out_dirs
FILES = [os.path.join(d, "driver_miss_rate.tsv") for d in table_out_dirs()]

# Scores derived without clinical labels → unsupervised/self-supervised (keyed by the short tool names in the file)
UNSUP = {
    # Conservation/evolution (no labels)
    "GERPpp_RS", "GERP_92_mammals", "phyloP17way_primate", "phyloP100way_vertebrate",
    "phyloP470way_mammalian", "phastCons17way_primate", "phastCons100way_vertebrate",
    "phastCons470way_mammalian", "SIFT", "SIFT4G", "PROVEAN", "MutationAssessor",
    # Generative/self-supervised language & evolution models
    "ESM1b", "popEVE", "AlphaMissense",
    # Population genetics (frequency/selection coefficient) & unsupervised integration
    "bStatistic", "MisFit_D", "MisFit_S", "Eigen-raw_coding", "Eigen-PC-raw_coding",
    # Based on population (primate common) labels (not clinical labels)
    "PrimateAI",
}
# Log of rationale for borderline decisions
RATIONALE = {
    "AlphaMissense": "unsup — protein LM + AlphaFold structure, weak supervision from population frequency (common=benign). Does not use clinical/ClinVar labels (Cheng 2023, Science).",
    "PrimateAI": "unsup — discriminative CNN using primate/population common variants as benign. Does not use clinical pathogenicity labels (ClinVar/HGMD) (Sundaram 2018, Nat Genet). [discriminative model but no clinical labels → unsup by principle]",
    "Eigen-raw_coding": "unsup — unsupervised meta-score integrating annotations via spectral/eigen-decomposition without labels (Ionita-Laza 2016, Nat Genet).",
    "Eigen-PC-raw_coding": "unsup — PC variant of Eigen, likewise unsupervised.",
    "MisFit_D": "unsup — population-genetics model estimating selection coefficient/deleteriousness from population allele frequencies. Does not use clinical labels.",
    "MisFit_S": "unsup — same as above (selection coefficient S).",
    "DANN": "supervised — neural network trained on the same discriminative training set as CADD (observed proxy-benign vs simulated proxy-deleterious) (Quang 2015). Consistent with classifying CADD as supervised. (If strictly 'clinical labels' only, CADD/DANN move to unsup — recomputation needed in that case.)",
}


def classify():
    df = pd.read_csv(FILES[0] if os.path.exists(FILES[0]) else FILES[1], sep="\t")
    df["type2"] = df["tool"].map(lambda t: "unsupervised" if t in UNSUP else "supervised")
    return df


def main():
    df = classify()
    for f in FILES:
        if os.path.exists(f):
            df.to_csv(f, sep="\t", index=False)

    print("=== Rationale for borderline decisions ===")
    for k, v in RATIONALE.items():
        print("  [%s] %s" % (k, v))

    print("\n=== Tools by type2 group ===")
    for t2 in ["unsupervised", "supervised"]:
        g = df[df.type2 == t2].sort_values("miss_rate(<0.5)", ascending=False)
        print("\n-- %s (n=%d) --" % (t2, len(g)))
        print("  tools:", ", ".join(g.tool.tolist()))
        print("  mean under-scoring(rec≥50 driver, rankscore<0.5): %.3f (%.1f%%)"
              % (g["miss_rate(<0.5)"].mean(), 100 * g["miss_rate(<0.5)"].mean()))

    u = df[df.type2 == "unsupervised"]["miss_rate(<0.5)"].mean()
    s = df[df.type2 == "supervised"]["miss_rate(<0.5)"].mean()
    print("\n=== R3 new numbers (replacing current 14.3/9.1) ===")
    print("  unsupervised %.1f%% vs supervised %.1f%% (diff %.1f%%p)" % (100 * u, 100 * s, 100 * (u - s)))
    print("\n[done] type2 added -> driver_miss_rate.tsv")


if __name__ == "__main__":
    main()
