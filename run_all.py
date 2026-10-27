"""
run_all.py — reproduce the full analysis from public sources.

    python run_all.py               # data download/harmonisation -> statistics -> figures

Environment variables (optional): AML_RAW (default data/raw), AML_OUT (default results).
Runtime ~6-8 min on a laptop (dominated by the 2,000-replicate bootstrap).
"""
import os, sys, runpy, time, argparse
ROOT = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("AML_RAW", os.path.join(ROOT, "data", "raw"))
os.environ.setdefault("AML_OUT", os.path.join(ROOT, "results"))
sys.path[:0] = [os.path.join(ROOT, "src")]


def step(name, path):
    t = time.time(); print(f"[{name}] running {os.path.relpath(path, ROOT)}")
    runpy.run_path(path, run_name="__main__")
    print(f"[{name}] done in {time.time()-t:.0f} s")


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__).parse_args()
    step("1/3 data", os.path.join(ROOT, "src", "build_data.py"))
    step("1b provenance", os.path.join(ROOT, "src", "provenance.py"))
    step("2/3 statistics", os.path.join(ROOT, "src", "analyze.py"))
    step("3/3 figures", os.path.join(ROOT, "src", "figures.py"))
