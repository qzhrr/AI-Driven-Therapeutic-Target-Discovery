"""
step03_synthetic_lethality — Synthetic lethality

Genotype-stratified SL against 9 recurrent GBM lesions + expression-context paralog scan.

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

# Synthetic-lethality (genotype + expression-context paralog scan)


import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import urllib.request
import ssl

# Load candidates
cand = pd.read_csv(result("step01_dependency", "glioma_candidates.csv"), index_col="gene")
candidate_genes = list(cand.index)

# Download DepMap 24Q4 files
ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE

import os
os.makedirs(DATA_RAW, exist_ok=True)

for fid, name in [(51065297,"Model.csv"), (51064916,"CRISPRInferredCommonEssentials.csv")]:
    url = f"https://ndownloader.figshare.com/files/{fid}"
    req = urllib.request.Request(url, headers={"User-Agent":"python-urllib"})
    data = urllib.request.urlopen(req, timeout=120, context=ctx).read()
    with open(str(raw(name)),"wb") as f: f.write(data)

# Download genomics files
for fid, name in [(51065747,"OmicsSomaticMutationsMatrixDamaging.csv"),
                  (51065750,"OmicsSomaticMutationsMatrixHotspot.csv"),
                  (51065303,"OmicsAbsoluteCNGene.csv")]:
    url = f"https://ndownloader.figshare.com/files/{fid}"
    req = urllib.request.Request(url, headers={"User-Agent":"python-urllib"})
    data = urllib.request.urlopen(req, timeout=300, context=ctx).read()
    with open(str(raw(name)),"wb") as f: f.write(data)

# Download CRISPRGeneEffect (only needed columns)
model = pd.read_csv(str(raw("Model.csv")))
glioma_ids = set(model[(model.OncotreeLineage=="CNS/Brain") &
                       (model.OncotreePrimaryDisease=="Diffuse Glioma")].ModelID)

# GBM lesion gene definitions
lesion_genes = ["CDKN2A","CDKN2B","MTAP","EGFR","PDGFRA","PTEN","NF1","RB1","TP53"]

def cols_for(path, genes):
    hdr = pd.read_csv(path, nrows=0)
    mp = {}
    for c in hdr.columns:
        sym = c.split(" (")[0]
        if sym in genes: mp[sym] = c
    return mp

cn_map = cols_for(str(raw("OmicsAbsoluteCNGene.csv")), set(lesion_genes))
dm_map = cols_for(str(raw("OmicsSomaticMutationsMatrixDamaging.csv")), set(lesion_genes))
ht_map = cols_for(str(raw("OmicsSomaticMutationsMatrixHotspot.csv")), set(lesion_genes))

def load_cols(path, cmap):
    idxcol = pd.read_csv(path, nrows=0).columns[0]
    df = pd.read_csv(path, index_col=0, usecols=[idxcol]+[cmap[g] for g in cmap])
    df.columns = [c.split(" (")[0] for c in df.columns]
    return df

cn = load_cols(str(raw("OmicsAbsoluteCNGene.csv")), cn_map)
dm = load_cols(str(raw("OmicsSomaticMutationsMatrixDamaging.csv")), dm_map)
ht = load_cols(str(raw("OmicsSomaticMutationsMatrixHotspot.csv")), ht_map)

# Build binary lesion status
all_models = cn.index.union(dm.index).union(ht.index)
L = pd.DataFrame(index=all_models)

def del_call(g):   return (cn[g] < 0.5)
def amp_call(g):   return (cn[g] >= 6.0)
def dmg_call(g):   return (dm[g] > 0) if g in dm.columns else pd.Series(False, index=dm.index)
def hot_call(g):   return (ht[g] > 0) if g in ht.columns else pd.Series(False, index=ht.index)

L["CDKN2A_del"] = del_call("CDKN2A").reindex(all_models)
L["CDKN2B_del"] = del_call("CDKN2B").reindex(all_models)
L["MTAP_del"]   = del_call("MTAP").reindex(all_models)
L["EGFR_amp"]   = amp_call("EGFR").reindex(all_models)
L["PDGFRA_amp"] = amp_call("PDGFRA").reindex(all_models)
L["PTEN_loss"]  = (del_call("PTEN")  | dmg_call("PTEN")).reindex(all_models)
L["NF1_loss"]   = (del_call("NF1")   | dmg_call("NF1")).reindex(all_models)
L["RB1_loss"]   = (del_call("RB1")   | dmg_call("RB1")).reindex(all_models)
L["TP53_mut"]   = (dmg_call("TP53")  | hot_call("TP53")).reindex(all_models)
L = L.fillna(False)

L_glioma = L[L.index.isin(glioma_ids)]

# Glioma prevalence per lesion
glioma_prev = {c: L_glioma[c].mean() for c in L.columns}

# Download CRISPRGeneEffect - only candidate columns
ge_hdr = pd.read_csv("https://ndownloader.figshare.com/files/51064667", nrows=0)
sym2col = {c.split(" (")[0]: c for c in ge_hdr.columns[1:]}
cand_cols = [sym2col[g] for g in candidate_genes if g in sym2col]
idxname = ge_hdr.columns[0]

# Stream only needed columns
req = urllib.request.Request("https://ndownloader.figshare.com/files/51064667",
                             headers={"User-Agent":"python-urllib"})
import io
raw = urllib.request.urlopen(req, timeout=600, context=ctx).read()
GE = pd.read_csv(io.BytesIO(raw), index_col=0, usecols=[idxname]+cand_cols)
GE.columns = [c.split(" (")[0] for c in GE.columns]

# Align
common = GE.index.intersection(L.index)
GEc = GE.loc[common]; Lc = L.loc[common]
glioma_mask_c = GEc.index.isin(glioma_ids)

def cohend(a, b):
    a = a.dropna(); b = b.dropna()
    if len(a) < 3 or len(b) < 3: return np.nan
    na, nb = len(a), len(b)
    sp = np.sqrt(((na-1)*a.std()**2 + (nb-1)*b.std()**2) / (na+nb-2))
    return (a.mean() - b.mean()) / sp if sp > 0 else np.nan

rows = []
for gene in GEc.columns:
    eff = GEc[gene]
    for lesion in Lc.columns:
        alt = eff[Lc[lesion].values]
        wt  = eff[~Lc[lesion].values]
        alt_d, wt_d = alt.dropna(), wt.dropna()
        if len(alt_d) < 5 or len(wt_d) < 5: continue
        d = cohend(alt_d, wt_d)
        try: u, p = stats.mannwhitneyu(alt_d, wt_d, alternative="two-sided")
        except: p = np.nan
        ge_g = eff[glioma_mask_c]; lg = Lc[lesion][glioma_mask_c]
        ga, gw = ge_g[lg.values].dropna(), ge_g[~lg.values].dropna()
        gdir = np.nan
        if len(ga) >= 3 and len(gw) >= 3: gdir = ga.mean() - gw.mean()
        rows.append([gene, lesion, alt_d.mean(), wt_d.mean(), d, p, len(alt_d), gdir])

sl = pd.DataFrame(rows, columns=["gene","lesion","alt_mean","wt_mean","cohens_d","p","n_alt","glioma_dir"])
sl["fdr"] = multipletests(sl["p"].fillna(1), method="fdr_bh")[1]
sl["glioma_lesion_prev"] = sl["lesion"].map(glioma_prev)

FLOOR = 0.08
sl["is_SL_pan"]         = (sl.cohens_d < -0.3) & (sl.fdr < 0.05)
sl["glioma_consistent"] = sl.glioma_dir < 0
sl["glioma_relevant"]   = sl.glioma_lesion_prev >= FLOOR
sl["is_SL_anchored"]    = sl.is_SL_pan & sl.glioma_consistent & sl.glioma_relevant

anch = sl[sl.is_SL_anchored].copy().sort_values("cohens_d")

TSG_LESIONS = {"CDKN2A_del","CDKN2B_del","PTEN_loss","NF1_loss","RB1_loss","TP53_mut"}
ONCOGENE_LESIONS = {"EGFR_amp","PDGFRA_amp"}
PASSENGER_LESIONS = {"MTAP_del"}
lesion_tsg = {"CDKN2A_del":"CDKN2A","CDKN2B_del":"CDKN2B","PTEN_loss":"PTEN",
              "NF1_loss":"NF1","RB1_loss":"RB1","TP53_mut":"TP53"}

feat = []
for g in candidate_genes:
    h = anch[anch.gene==g]
    if len(h):
        b = h.iloc[0]
        feat.append([g, -b.cohens_d, b.lesion, b.fdr, len(h), ";".join(h.lesion)])
    else:
        feat.append([g, 0.0, "", 1.0, 0, ""])
slf = pd.DataFrame(feat, columns=["gene","sl_genotype_max_effect_size","sl_best_context",
                                   "sl_fdr","sl_n_contexts","sl_all_contexts"]).set_index("gene")

def ctx_type(l):
    return "TSG-loss" if l in TSG_LESIONS else "oncogene-amp" if l in ONCOGENE_LESIONS else "passenger-del(PRMT5-axis)" if l in PASSENGER_LESIONS else ""
slf["sl_context_type"] = slf.sl_best_context.map(lambda l: ctx_type(l) if l else "")
slf["is_SL_partner_of_lost_TSG"] = False
slf["lost_TSG_partner"] = ""
for g in slf.index:
    th = anch[(anch.gene==g) & (anch.lesion.isin(TSG_LESIONS))]
    if len(th):
        slf.loc[g,"is_SL_partner_of_lost_TSG"] = True
        slf.loc[g,"lost_TSG_partner"] = ";".join(sorted(set(lesion_tsg[l] for l in th.lesion)))

slf.to_csv(result("step03_synthetic_lethality", "sl_candidate_features.csv"))


import pandas as pd
import numpy as np
import re
from scipy import stats

# Load data
rank = pd.read_csv('/Users/qzhang11/Desktop/hackathon/Hackathon 2/GBM_target_discovery/final_candidate_ranking.csv')
coh = pd.read_csv('{{artifact:MISSING:e57fba62-93dc-446d-9602-8ac7c9f2cfa7}}')
glioma_ids = set(coh.ModelID)

candidates = rank['gene'].tolist()

def load_by_symbol(path, symbols):
    h = pd.read_csv(path, nrows=0).columns.tolist()
    idc = h[0]
    name_map = {}
    for s in symbols:
        for c in h:
            if c.split(' (')[0]==s:
                name_map[c]=s; break
    d = pd.read_csv(path, usecols=[idc]+list(name_map))
    d = d.rename(columns={idc:'ModelID', **name_map})
    return d

# Load candidate gene-effect columns
dep = load_by_symbol(raw("CRISPRGeneEffect.csv"), candidates)
dep_genes = [c for c in dep.columns if c!='ModelID']
print(f"candidate gene-effect cols found: {len(dep_genes)}/{len(candidates)}")

# Load FULL expression matrix
eh = pd.read_csv(raw("OmicsExpression.csv"), nrows=0).columns.tolist()
expr_full = pd.read_csv(raw("OmicsExpression.csv"))
expr_full = expr_full.rename(columns={eh[0]:'ModelID'})
expr_genes = [c.split(' (')[0] for c in expr_full.columns if c!='ModelID']
expr_full.columns = ['ModelID']+expr_genes
print(f"expression matrix: {expr_full.shape}")

# Intersect on common cell lines
common = sorted(set(dep['ModelID']) & set(expr_full['ModelID']))
depm = dep.set_index('ModelID').loc[common, dep_genes].astype('float32')
exprm = expr_full.set_index('ModelID').loc[common, :].astype('float32')
del expr_full
print(f"common cell lines: {len(common)} | dep {depm.shape} | expr {exprm.shape}")

# Vectorized Pearson correlation
Dz = (depm - depm.mean()) / depm.std()
Ez = (exprm - exprm.mean()) / exprm.std()
Dz = Dz.fillna(0).values.astype('float32')
Ez = Ez.fillna(0).values.astype('float32')
n = Dz.shape[0]
R = (Dz.T @ Ez) / n
print(f"correlation matrix {R.shape}")

dep_idx = {g:i for i,g in enumerate(dep_genes)}
expr_idx = {g:i for i,g in enumerate(exprm.columns)}

def expr_sl_delta(dep_gene, partner):
    di, pj = dep_idx[dep_gene], expr_idx[partner]
    d = depm.iloc[:, di].values; e = exprm.iloc[:, pj].values
    m = ~(np.isnan(d)|np.isnan(e))
    d,e = d[m], e[m]
    lo = d[e <= np.percentile(e,33)]; hi = d[e >= np.percentile(e,67)]
    sp = np.sqrt((lo.var(ddof=1)+hi.var(ddof=1))/2)
    dd = (lo.mean()-hi.mean())/sp
    _,p = stats.mannwhitneyu(lo,hi,alternative='less')
    return dd, p, lo.mean(), hi.mean()

# Scan every candidate: best expression-SL partner (max positive r), with delta+p
rows=[]
for g in dep_genes:
    row = R[dep_idx[g]].copy()
    row[expr_idx.get(g,-1)] = -9  # exclude self
    j = int(np.argmax(row)); partner = exprm.columns[j]; r = float(row[j])
    dd,p,lom,him = expr_sl_delta(g, partner)
    rows.append(dict(gene=g, top_expr_partner=partner, r=round(r,3),
                     dep_lowExpr=round(lom,3), dep_highExpr=round(him,3),
                     cohens_d=round(dd,3), p=p))
scan = pd.DataFrame(rows)
scan['sig'] = (scan['r']>=0.3) & (scan['cohens_d']<=-0.4) & (scan['p']<0.01)
scan = scan.sort_values('r', ascending=False)

print(f"candidates with a strong expression-context SL partner (r>=0.3, d<=-0.4, p<0.01): {scan['sig'].sum()}/283")

# Flag likely paralogs by shared symbol root
def root(sym):
    base = re.sub(r'\d+[A-Z]?$','',sym)
    return base

scan['gene_root'] = scan['gene'].map(root)
scan['partner_root'] = scan['top_expr_partner'].map(root)
scan['likely_paralog'] = (scan['gene_root']==scan['partner_root']) & (scan['gene']!=scan['top_expr_partner'])

# Literature verdicts for the 9 name-matched paralog SL pairs
verdicts = {
 'VRK1':   ('yes','validated','VRK1←VRK2: Cancer Res 2022 (Tango/TANDEM); VRK2-methylated GBM; isogenic + in-vivo validated'),
 'FAM50A': ('yes','validated','FAM50A←FAM50B: Thompson/Parrish 2021 Nat Commun; PaCT Cell Rep 2022; dual-KO validated'),
 'RPP25L': ('yes','validated','RPP25L←RPP25: PaCT Cell Rep 2022; GBM-specific Neuro-Oncology 2022 EXTH-61 (RPP25 promoter methylation)'),
 'CDS2':   ('yes','validated','CDS2←CDS1: Nature Genetics 2025 (uveal melanoma + pan-cancer); in-vivo + rescue validated'),
 'EAF1':   ('yes','screen','EAF1←EAF2: reported experimentally-verified paralog SL (SL-RFM 2023; grouped with CDK4/6, COPG1/2)'),
 'ELMO2':  ('yes','validated','ELMO2←ELMO3: Nature Communications 2026 (mesenchymal-like/EGFR-TKI-resistant NSCLC; ZEB1 represses ELMO3)'),
 'PTK2':   ('partial','redundancy','PTK2←PTK2B: FAK/PYK2 functional redundancy documented; FAK is SL with GNAQ (uveal mel); not a crisp PTK2B-low→FAK SL landmark'),
 'GOLT1B': ('no','none','GOLT1B←GOLT1A: no SL report found; only prognostic/oncogene studies — candidate NOVEL from this scan'),
 'CHMP3':  ('not_true_paralog','none','CHMP3←CHMP2A: name-root false match — CHMP2A paralog is CHMP2B; CHMP3=Vps24 is a distinct ESCRT-III subunit. Correlation = obligate complex co-membership, not paralog SL; not reported'),
}

def status(row):
    g = row['gene']
    if g in verdicts:
        return pd.Series(verdicts[g], index=['pair_reported_before','report_evidence_level','reported_reference'])
    if row['sig'] and not row['likely_paralog']:
        return pd.Series(['not_assessed','', 'Non-paralog top partner — likely lineage/co-expression confounder; SL not literature-assessed'],
                         index=['pair_reported_before','report_evidence_level','reported_reference'])
    return pd.Series(['','',''], index=['pair_reported_before','report_evidence_level','reported_reference'])

scan[['pair_reported_before','report_evidence_level','reported_reference']] = scan.apply(status, axis=1)

scan = scan.merge(rank[['gene','rank_all']], on='gene', how='left')
# reorder: put rank_all near the front
cols = ['gene','top_expr_partner','r','cohens_d','p','sig','likely_paralog','rank_all',
        'pair_reported_before','report_evidence_level','reported_reference',
        'dep_lowExpr','dep_highExpr','gene_root','partner_root']
cols = [c for c in cols if c in scan.columns]
scan = scan[cols]
scan.to_csv(result("step03_synthetic_lethality", "expression_SL_scan_candidates.csv"), index=False)
print("saved expression_SL_scan_candidates.csv", scan.shape)
