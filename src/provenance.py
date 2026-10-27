"""provenance.py — record what was downloaded, from where, and when.

Writes <AML_OUT>/data_provenance.json with: UTC timestamp, cBioPortal API build info,
the AMLSG GitHub commit that the requested ref resolved to, and SHA-256 checksums of
every raw file. Compare checksums against the published data_provenance.json to
confirm that you analysed byte-identical inputs.
"""
import os, json, hashlib, datetime, urllib.request

RAW = os.environ.get("AML_RAW", "raw"); OUT = os.environ.get("AML_OUT", "out")
REF = os.environ.get("AMLSG_REF", "a3f96ba727bfa3f30e47d9419a9faf1caab44d51")


def _get(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=60) as r:
            return json.loads(r.read().decode())
    except Exception as e:  # provenance must never break the pipeline
        return {"error": repr(e)[:200]}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if __name__ == "__main__":
    files = {}
    for d, _, fs in os.walk(RAW):
        for fn in sorted(fs):
            p = os.path.join(d, fn)
            files[os.path.relpath(p, RAW).replace(os.sep, "/")] = {"bytes": os.path.getsize(p), "sha256": sha256(p)}
    c = _get(f"https://api.github.com/repos/gerstung-lab/AML-multistage/commits/{REF}")
    rec = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "cbioportal_api_info": _get("https://www.cbioportal.org/api/info"),
        "cbioportal_studies": ["laml_tcga_pub", "aml_ohsu_2022", "aml_ohsu_2018 (classifier validation only)"],
        "amlsg_repository": "https://github.com/gerstung-lab/AML-multistage",
        "amlsg_ref_requested": REF,
        "amlsg_commit_resolved": c.get("sha", c),
        "amlsg_commit_date": (c.get("commit") or {}).get("committer", {}).get("date") if isinstance(c, dict) else None,
        "raw_files": files,
    }
    os.makedirs(OUT, exist_ok=True)
    json.dump(rec, open(os.path.join(OUT, "data_provenance.json"), "w"), indent=2)
    print(f"provenance: {len(files)} raw files hashed; AMLSG commit {str(rec['amlsg_commit_resolved'])[:12]}")
