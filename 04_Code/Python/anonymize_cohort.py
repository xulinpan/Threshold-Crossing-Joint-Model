#!/usr/bin/env python3
"""
anonymize_cohort.py — irreversible re-anonymisation of the CML cohort.

WHY THIS IS NEEDED
------------------
The existing coded identifiers P0001...P0096 are assigned in *source-file
order*. Anyone holding the original spreadsheet can therefore map row N back to
patient P000N without ever seeing `private/patient_key.csv`. That makes the
current files pseudonymised, not anonymised.

Additionally, in a single-centre cohort of 87 patients, the combination
(exact age, sex, treatment duration) is a quasi-identifier.

WHAT THIS DOES
--------------
1. Re-maps every patient to an opaque, non-sequential ID (e.g. PT-3F9A) drawn
   from a random permutation, and re-assigns the Stan index `patient_num`
   accordingly, so neither identifier preserves source order.
2. Re-sorts every output file by the new index, destroying original row order.
3. Coarsens quasi-identifiers: age -> 10-year bands, duration -> integer years.
4. Writes anonymised analysis files to 03_Data/Anonymized/.
5. Writes the crosswalk to 03_Data/Processed/private/ ONLY. Delete or secure
   that file to make the anonymisation irreversible.

Original files are never modified.
"""
import csv, os, random, hashlib, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PROC = os.path.join(ROOT, "03_Data", "Processed")
OUT  = os.path.join(ROOT, "03_Data", "Anonymized")
PRIV = os.path.join(PROC, "private")
os.makedirs(OUT, exist_ok=True); os.makedirs(PRIV, exist_ok=True)

# Set to None for a non-reproducible (more private) mapping; an integer makes
# the mapping reproducible for the authors but is still unguessable externally.
SEED = None

FILES = [
    "real_longitudinal_analysis.csv",
    "real_interval_survival_analysis.csv",
    "real_patient_level_analysis.csv",
    "real_patient_level_complete_covariates.csv",
    "clean_patient_covariates.csv",
]

def read(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f); return list(r), r.fieldnames

# ---- 1. build the mapping from the patient-level file ------------------
base, _ = read(os.path.join(PROC, "real_patient_level_analysis.csv"))
ids = sorted({r["patient_id"] for r in base})
rng = random.Random(SEED)
order = ids[:]; rng.shuffle(order)

used, id_map, num_map = set(), {}, {}
for new_num, old in enumerate(order, start=1):
    while True:
        tok = "PT-" + hashlib.sha256(f"{old}{rng.random()}".encode()).hexdigest()[:4].upper()
        if tok not in used: used.add(tok); break
    id_map[old] = tok
    num_map[old] = new_num
print(f"mapped {len(id_map)} patients to opaque IDs")

# ---- 2. quasi-identifier coarsening -----------------------------------
def band_age(v):
    try: a = float(v)
    except (TypeError, ValueError): return "NA"
    lo = int(a // 10) * 10
    return f"{lo}-{lo+9}"

def coarsen(row):
    if "age" in row:
        row["age_band"] = band_age(row.pop("age"))
    if "duration_years" in row and row["duration_years"] not in ("", "NA"):
        try: row["duration_years"] = str(int(round(float(row["duration_years"]))))
        except ValueError: pass
    return row

# ---- 3. apply to every file -------------------------------------------
for fn in FILES:
    src = os.path.join(PROC, fn)
    if not os.path.exists(src):
        print(f"  skip (absent): {fn}"); continue
    rows, cols = read(src)
    out = []
    for r in rows:
        pid = r.get("patient_id")
        if pid not in id_map:      # patient not in the analysis cohort
            continue
        r["patient_id"] = id_map[pid]
        if "patient_num" in r: r["patient_num"] = str(num_map[pid])
        out.append(coarsen(r))
    # destroy original ordering
    out.sort(key=lambda r: (int(r.get("patient_num", 0)),
                            float(r.get("t_months", 0) or 0),
                            int(r.get("visit_index", 0) or 0)))
    cols = [c for c in cols if c != "age"]
    if any("age_band" in r for r in out) and "age_band" not in cols:
        cols.append("age_band")
    dst = os.path.join(OUT, fn.replace("real_", "anon_").replace("clean_", "anon_"))
    with open(dst, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(out)
    print(f"  wrote {os.path.basename(dst)}  ({len(out)} rows)")

# ---- 4. crosswalk -> private only --------------------------------------
xw = os.path.join(PRIV, "anonymization_crosswalk.csv")
with open(xw, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["original_patient_id", "anon_patient_id", "anon_patient_num"])
    for old in ids: w.writerow([old, id_map[old], num_map[old]])
print(f"\ncrosswalk -> {xw}")
print("DELETE OR SECURE THAT FILE to make the anonymisation irreversible.")

# ---- 5. verification ---------------------------------------------------
bad = 0
for fn in os.listdir(OUT):
    with open(os.path.join(OUT, fn), encoding="utf-8") as f:
        txt = f.read()
    if any(ord(c) > 0x4e00 and ord(c) < 0x9fff for c in txt):
        print(f"  FAIL: CJK characters in {fn}"); bad += 1
    if "P00" in txt:
        print(f"  FAIL: original-style ID in {fn}"); bad += 1
print("verification:", "PASS - no names, no original IDs" if bad == 0 else f"{bad} problem(s)")
