"""check_reproduction.py — compare a fresh run (results/results.json) with the published reference values.

    python tests/check_reproduction.py [path/to/results.json]

Counts must match exactly; Cox coefficients within 1e-6; bootstrap summaries within 1e-6
(same seed, 20261026). A mismatch in counts almost always means the upstream source data
changed — compare results/data_provenance.json checksums first.
"""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ref = json.load(open(os.path.join(ROOT, "tests", "reference_results.json")))
new = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "results.json")))
pc = {d["covariate"]: d["coef"] for d in new["primary_coefficients"]}
bs = new["S7_bootstrap_interaction"]
fails = []
for k in ["n_primary", "deaths", "n_by_cohort", "cell_events"]:
    if new[k] != ref[k]:
        fails.append(f"{k}: {new[k]} != {ref[k]}")
for k, v in ref["coef"].items():
    if abs(pc.get(k, float("nan")) - v) > 1e-6 or pc.get(k) is None:
        fails.append(f"coef {k}: {pc.get(k)} vs {v}")
for k, nk in [("bootstrap_ratio_median", "ratio_median"), ("bootstrap_lo", "lo"), ("bootstrap_hi", "hi")]:
    if abs(bs[nk] - ref[k]) > 1e-6:
        fails.append(f"{k}: {bs[nk]} vs {ref[k]}")
print("REPRODUCED: all reference values match" if not fails else "MISMATCH:\n  " + "\n  ".join(fails))
sys.exit(1 if fails else 0)
