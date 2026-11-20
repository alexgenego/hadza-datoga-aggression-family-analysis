"""Assemble manuscript tables, figures and Supplementary registries.

This module is intentionally downstream-only: it does not recompute inferential
statistics. It reads outputs produced by scripts 01-09 and converts them into
manuscript-ready tables, central registries and selected figures.

Primary products:
- Table 1: sample flow / model-specific N;
- Table 4: LASSO summary;
- Table 5: full-family-graph + grouped RF integrated summary;
- Table 6: kinship-scale profiles for DAT1, HTR1A, HTR2A and SERTPr;
- Figure 6: adult genetic hierarchy / forest-style overview;
- Supplementary Table S_Kinship;
- Supplementary Table S_RF;
- Supplementary provenance manifest.

Existing manuscript Figures 1-5 are not re-analysed here; this script may copy
or regenerate only those panels whose source data are available in public
outputs.
"""

from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TARGET_KINSHIP_LOCI=["DAT1","HTR1A","HTR2A","SERTPr"]

def read_csv_if(path):
    p=Path(path)
    return pd.read_csv(p) if p.exists() else None

def ensure_dir(path):
    p=Path(path); p.mkdir(parents=True,exist_ok=True); return p

def concat_glob(root,pattern):
    files=sorted(Path(root).glob(pattern))
    dfs=[]
    for f in files:
        try:
            x=pd.read_csv(f); x.insert(0,"source_file",str(f)); dfs.append(x)
        except Exception:
            pass
    return pd.concat(dfs,ignore_index=True) if dfs else pd.DataFrame()

def table1_sample_flow(outputs_root):
    root=Path(outputs_root)
    rows=[]
    for f in sorted(root.glob("harmonized/*_sample_flow.csv")):
        x=pd.read_csv(f)
        dataset=f.name.replace("_sample_flow.csv","")
        for _,r in x.iterrows():
            rows.append({"dataset":dataset,"variable":r["variable"],
                         "N_nonmissing":r["N_nonmissing"],"N_missing":r["N_missing"]})
    # add analysis Ns from major downstream summaries where available
    extra=[]
    candidates=[
        ("lasso",root/"lasso"),
        ("random_forest",root/"random_forest"),
        ("family_graph",root/"family_graph"),
    ]
    for kind,p in candidates:
        if not p.exists(): continue
        for f in sorted(p.glob("*summary.csv")):
            x=pd.read_csv(f)
            ncols=[c for c in x.columns if c in {"N","n_nodes","n_families","n_components","n_parent_clusters"}]
            for i,r in x.iterrows():
                rec={"dataset":kind,"variable":f.stem}
                for c in ncols: rec[c]=r[c]
                extra.append(rec)
    base=pd.DataFrame(rows)
    if extra:
        extra=pd.DataFrame(extra)
        return pd.concat([base,extra],ignore_index=True,sort=False)
    return base

def table4_lasso(outputs_root):
    root=Path(outputs_root)/"lasso"
    summaries=concat_glob(root,"*_summary.csv")
    stabilities=concat_glob(root,"*_stability.csv")
    if summaries.empty and stabilities.empty: return pd.DataFrame()
    if not stabilities.empty:
        keep=stabilities.loc[
            (stabilities.get("selection_frequency_lambda_min",0)>=0.10) |
            (stabilities.get("selection_frequency_lambda_1se",0)>0)
        ].copy()
    else:
        keep=pd.DataFrame()
    # compact manuscript table uses summary rows plus a textual selected-marker digest
    rows=[]
    if not summaries.empty:
        for _,r in summaries.iterrows():
            panel=r.get("panel","")
            ph=r.get("phenotype","")
            sol=r.get("solution","")
            freq=""
            if not keep.empty:
                hit=keep[
                    keep["source_file"].astype(str).str.contains(str(panel),regex=False) &
                    keep["source_file"].astype(str).str.contains(str(ph),regex=False)
                ]
                if len(hit):
                    hit=hit.sort_values("selection_frequency_lambda_min",ascending=False)
                    freq="; ".join(
                        f"{z.locus} {z.selection_frequency_lambda_min:.3f}"
                        for z in hit.itertuples()
                    )
            rows.append({
                "panel":panel,"phenotype":ph,"N":r.get("N",np.nan),
                "solution":sol,"OOF_delta_R2":r.get("OOF_delta_R2",np.nan),
                "mean_selected":r.get("mean_selected",np.nan),
                "lambda_min_selection_ge_0.10":freq
            })
    return pd.DataFrame(rows)

def table5_family_rf(outputs_root):
    root=Path(outputs_root)
    rows=[]
    fg=read_csv_if(root/"family_graph_permutations.csv")
    if fg is not None and not fg.empty:
        for _,r in fg.iterrows():
            rows.append({
                "analysis":"Family graph",
                "phenotype":r.get("trait",""),
                "scope":r.get("scope",""),
                "N":r.get("n_nodes",np.nan),
                "families":r.get("n_components",np.nan),
                "effect":r.get("dispersion_ratio",np.nan),
                "metric":"dispersion_ratio",
                "p":r.get("permutation_p",np.nan),
                "q":r.get("BH_q",np.nan)
            })
    rf=read_csv_if(root/"random_forest"/"rf_summary.csv")
    if rf is not None and not rf.empty:
        for _,r in rf.iterrows():
            if r.get("model","")!="genetic": continue
            rows.append({
                "analysis":"Grouped Random Forest",
                "phenotype":r.get("phenotype",""),
                "scope":f"{r.get('cohort','')}; leaf={r.get('min_samples_leaf',np.nan)}",
                "N":r.get("N",np.nan),
                "families":r.get("n_families",np.nan),
                "effect":r.get("delta_R2_genetic_minus_base",np.nan),
                "metric":"delta_R2",
                "p":r.get("paired_signflip_p",np.nan),
                "q":np.nan
            })
    return pd.DataFrame(rows)

def table6_kinship(outputs_root):
    x=read_csv_if(Path(outputs_root)/"kinship_scale_decomposition.csv")
    if x is None or x.empty: return pd.DataFrame()
    y=x.loc[x["locus"].isin(TARGET_KINSHIP_LOCI)].copy()
    order=["sibling","mother_defined_offspring","father_defined_offspring",
           "direct_mother","direct_father","mixed_generation","full_graph"]
    y["scale_order"]=pd.Categorical(y["scale"],categories=order,ordered=True)
    return y.sort_values(["locus","phenotype","scale_order"]).drop(columns=["scale_order"])

def supplementary_kinship(outputs_root):
    x=read_csv_if(Path(outputs_root)/"kinship_scale_decomposition.csv")
    return pd.DataFrame() if x is None else x

def supplementary_rf(outputs_root):
    root=Path(outputs_root)/"random_forest"
    parts=[]
    summary=read_csv_if(root/"rf_summary.csv")
    if summary is not None:
        z=summary.copy(); z.insert(0,"record_type","summary"); parts.append(z)
    for f in sorted(root.glob("*_family_errors.csv")):
        z=pd.read_csv(f); z.insert(0,"record_type","family_error"); z.insert(1,"source_file",f.name); parts.append(z)
    return pd.concat(parts,ignore_index=True,sort=False) if parts else pd.DataFrame()

def adult_genetic_hierarchy(outputs_root):
    rows=[]
    # locus-specific/cross-age results
    ag=read_csv_if(Path(outputs_root)/"age_genotype_models.csv")
    if ag is not None:
        for _,r in ag.iterrows():
            rows.append({"label":f"{r.get('locus','')} {r.get('phenotype','')} {r.get('analysis','')}",
                         "effect":r.get("beta",np.nan),"p":r.get("p",np.nan),
                         "q":r.get("q",np.nan),"source":"age_genotype"})
    # multilocus block
    mb=read_csv_if(Path(outputs_root)/"multilocus_dispersion"/"multilocus_block.csv")
    if mb is not None:
        for _,r in mb.iterrows():
            rows.append({"label":f"Core block {r.get('phenotype','')}",
                         "effect":r.get("partial_delta_R2",np.nan),"p":r.get("permutation_p",np.nan),
                         "q":np.nan,"source":"multilocus"})
    # kinship
    ks=read_csv_if(Path(outputs_root)/"kinship_scale_decomposition.csv")
    if ks is not None:
        for _,r in ks.loc[ks["locus"].isin(TARGET_KINSHIP_LOCI)].iterrows():
            rows.append({"label":f"{r.get('locus','')} {r.get('phenotype','')} {r.get('scale','')}",
                         "effect":r.get("rho",np.nan),"p":r.get("permutation_p",np.nan),
                         "q":r.get("q9",np.nan),"source":"kinship"})
    return pd.DataFrame(rows)

def make_figure6(registry,outfile):
    if registry.empty: return
    x=registry.dropna(subset=["effect"]).copy()
    if x.empty:return
    # keep a readable hierarchy: strongest supported + top nominal entries
    x["support_rank"]=np.where(x["q"].fillna(1)<=0.05,0,
                       np.where(x["p"].fillna(1)<=0.05,1,2))
    x=x.sort_values(["support_rank","p"],na_position="last").head(20).iloc[::-1]
    fig,ax=plt.subplots(figsize=(9,max(5,0.32*len(x))))
    yy=np.arange(len(x))
    ax.scatter(x["effect"],yy)
    ax.axvline(0,linewidth=1)
    ax.set_yticks(yy)
    ax.set_yticklabels(x["label"])
    ax.set_xlabel("Effect estimate (scale depends on analysis)")
    ax.set_title("Adult genetic hierarchy: confirmed and nominal signals")
    fig.tight_layout()
    fig.savefig(outfile,dpi=300,bbox_inches="tight")
    plt.close(fig)

def make_family_component_figure(outputs_root,outfile):
    nodes=read_csv_if(Path(outputs_root)/"family_graph"/"social_graph_nodes.csv")
    if nodes is None or nodes.empty:return
    sizes=(nodes.groupby(["population","component_id"]).size()
           .reset_index(name="n").sort_values(["population","n"],ascending=[True,False]))
    fig,ax=plt.subplots(figsize=(9,5))
    labels=sizes["population"].astype(str)+":"+sizes["component_id"].astype(str)
    ax.bar(np.arange(len(sizes)),sizes["n"])
    ax.set_xticks(np.arange(len(sizes))); ax.set_xticklabels(labels,rotation=90)
    ax.set_ylabel("Adults per connected component")
    ax.set_title("Strict adult family graph: connected-component sizes")
    fig.tight_layout(); fig.savefig(outfile,dpi=300,bbox_inches="tight"); plt.close(fig)

def manifest(outputs_root):
    root=Path(outputs_root)
    rows=[]
    for f in sorted(root.rglob("*.csv")):
        try:
            x=pd.read_csv(f)
            rows.append({"file":str(f.relative_to(root)),"rows":len(x),"columns":len(x.columns),
                         "column_names":"|".join(map(str,x.columns))})
        except Exception as e:
            rows.append({"file":str(f.relative_to(root)),"rows":np.nan,"columns":np.nan,
                         "column_names":"","read_error":str(e)})
    return pd.DataFrame(rows)

def write_xlsx(path,sheets):
    with pd.ExcelWriter(path,engine="openpyxl") as w:
        wrote=False
        for name,df in sheets.items():
            if df is not None and not df.empty:
                df.to_excel(w,sheet_name=name[:31],index=False); wrote=True
        if not wrote:
            pd.DataFrame({"status":["No source outputs available yet"]}).to_excel(w,sheet_name="README",index=False)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--outputs-root",default="outputs")
    p.add_argument("--manuscript-dir",default="outputs/manuscript")
    a=p.parse_args()
    out=ensure_dir(a.manuscript_dir)
    tables={
        "Table1_sample_flow":table1_sample_flow(a.outputs_root),
        "Table4_LASSO":table4_lasso(a.outputs_root),
        "Table5_family_RF":table5_family_rf(a.outputs_root),
        "Table6_kinship_scale":table6_kinship(a.outputs_root),
    }
    write_xlsx(out/"manuscript_tables.xlsx",tables)
    supp={
        "S_Kinship":supplementary_kinship(a.outputs_root),
        "S_RF":supplementary_rf(a.outputs_root),
        "Provenance":manifest(a.outputs_root),
    }
    write_xlsx(out/"supplementary_registries.xlsx",supp)
    reg=adult_genetic_hierarchy(a.outputs_root)
    reg.to_csv(out/"adult_genetic_hierarchy_registry.csv",index=False)
    make_figure6(reg,out/"Figure6_adult_genetic_hierarchy.png")
    make_family_component_figure(a.outputs_root,out/"Supplementary_Figure_S2_component_sizes.png")
    with open(out/"build_manifest.json","w",encoding="utf-8") as fh:
        json.dump({
            "tables":[k for k,v in tables.items() if v is not None and not v.empty],
            "supplementary":[k for k,v in supp.items() if v is not None and not v.empty],
            "figures":[p.name for p in out.glob("*.png")]
        },fh,ensure_ascii=False,indent=2)

if __name__=="__main__":
    main()
