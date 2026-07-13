"""
step07_tractability — Tractability + druggability + novelty

Open Targets / ChEMBL / DGIdb / UniProt -> graded tractability + separate novelty axis.

Reconstructed from the Claude Science artifact lineage of the GBM small-molecule
target-discovery project. See docs/METHODS.md for the scientific description of
this step and README.md for how paths and run modes work.

Inputs / outputs use the repository-relative helpers in common.py
(raw/ cache/ result). Run cached (default) or live (GBM_LIVE=1); see common.py.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    ROOT, DATA, DATA_RAW, DATA_CACHE, RESULTS, STEP_DIRS,
    LIVE, use_cache, raw, cache, result, require, ensure_dirs,
)
from adapters import host  # portable shim for host.llm / host.mcp / host.current_model

ensure_dirs()

# Tractability + druggability + novelty


import json, pandas as pd, numpy as np

cand = pd.read_csv(result("step01_dependency", "glioma_candidates.csv")).set_index("gene")
ot = json.load(open(cache("ot_tractability.json")))
dg = json.load(open(cache("dgidb.json")))
ch = json.load(open(cache("chembl.json")))
up = json.load(open(cache("uniprot.json")))
clean = json.load(open(cache("chembl_clean_potent.json")))

rows = []
for g in cand.index:
    o = ot.get(g, {}); d = dg.get(g, {}); c = ch.get(g, {}); u = up.get(g, {})
    n_potent = clean.get(g, c.get("chembl_n_potent", 0))
    tract = 0.0
    if o.get("ot_sm_hq_pocket") or o.get("ot_sm_hq_ligand"): tract += 0.35
    if o.get("ot_sm_struct_ligand"): tract += 0.25
    if o.get("ot_sm_druggable_family"): tract += 0.20
    if n_potent >= 10: tract += 0.20
    tract = min(tract, 1.0)
    has_approved = bool(o.get("ot_sm_approved_drug")) or (d.get("dgidb_n_approved_drugs", 0) > 0)
    has_clinical = bool(o.get("ot_sm_advanced_clinical")) or bool(o.get("ot_sm_phase1")) or (o.get("ot_drug_count", 0) > 0)
    novelty = 0.0 if has_approved else (0.5 if has_clinical else 1.0)
    acts = set(c.get("chembl_action_types", []))
    has_inhibitor_moa = bool(acts & {"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR", "BLOCKER"})
    rows.append({"gene": g, "tractability_score": round(tract, 3), "novelty_score": novelty,
        "ot_sm_approved_drug": bool(o.get("ot_sm_approved_drug")),
        "ot_sm_druggable_family": bool(o.get("ot_sm_druggable_family")),
        "ot_sm_hq_pocket": bool(o.get("ot_sm_hq_pocket")),
        "ot_sm_struct_ligand": bool(o.get("ot_sm_struct_ligand")),
        "ot_drug_count": o.get("ot_drug_count", 0),
        "dgidb_n_drugs": d.get("dgidb_n_drugs", 0), "dgidb_n_approved_drugs": d.get("dgidb_n_approved_drugs", 0),
        "chembl_target_id": c.get("chembl_target_id"), "chembl_n_potent": n_potent,
        "chembl_action_types": ";".join(sorted(acts)) if acts else "", "has_inhibitor_moa": has_inhibitor_moa,
        "uniprot_acc": u.get("uniprot_acc"), "uniprot_kinase": bool(u.get("uniprot_kinase")),
        "uniprot_enzyme": bool(u.get("uniprot_enzyme")), "uniprot_membrane": bool(u.get("uniprot_membrane")),
        "uniprot_loc": ";".join(u.get("uniprot_loc", [])[:3])})

tr = pd.DataFrame(rows).set_index("gene")
tr.to_csv(result("step07_tractability", "tractability_features.csv"))
