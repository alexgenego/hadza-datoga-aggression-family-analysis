"""Kinship-scale decomposition of locus-specific family associations.

Scales implemented:
- sibling components;
- mother-defined and father-defined offspring lines;
- direct parent genotype -> mean offspring phenotype;
- mixed-generation components;
- full connected pedigree components.

Adult phenotypes are residualized within population as
phenotype ~ sex + centered age + centered age^2 using behavior-wave age >=18.

For family-mean scales, a family unit enters the screen when at least one
member has both the locus and adjusted phenotype available. This reproduces
the manuscript scale definitions (for example, DAT1 mother-defined verbal:
10 lines / 19 offspring). Locus dosage is shuffled among eligible adults
within population while family membership is fixed. For direct-parent tests, parent
genotypes are shuffled among eligible parents within population. Two-sided
permutation p-values are calculated from |Spearman rho|. BH-FDR q9 is then
calculated across loci within each phenotype x kinship-scale screen.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

CORE9=["DRD4P","SERTPr","SERT_In2","DRD4E3","DAT1","DRD2Taq1","COMT","HTR1A","HTR2A"]

def read_table(path):
    path=Path(path)
    if path.suffix.lower() in {".xlsx",".xls"}: return pd.read_excel(path)
    return pd.read_csv(path,sep="\t" if path.suffix.lower() in {".tsv",".txt"} else ",")

def bh_fdr(p):
    p=np.asarray(list(p),float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    vals=p[ok]; order=np.argsort(vals); ranked=vals[order]
    q=ranked*len(ranked)/np.arange(1,len(ranked)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    inv=np.empty_like(order); inv[order]=np.arange(len(order))
    out[np.where(ok)[0]]=q[inv]; return out

def encode_sex(s):
    if not pd.api.types.is_numeric_dtype(s): s=s.astype("string")
    vals=sorted(s.dropna().unique())
    if len(vals)>2: raise ValueError(f"sex has >2 categories: {vals}")
    return s.map({v:i for i,v in enumerate(vals)}).astype(float)

def age_adjusted_residuals(nodes,phenotype):
    out=pd.Series(np.nan,index=nodes.index,dtype=float)
    for _,idx in nodes.groupby("population",dropna=False).groups.items():
        x=nodes.loc[idx].copy()
        x["age"]=pd.to_numeric(x["age"],errors="coerce")
        x=x.loc[x["age"].ge(18)].dropna(subset=[phenotype,"sex","age"])
        if len(x)<4: continue
        age=x["age"].astype(float); ac=age-age.mean()
        X=pd.DataFrame({"sex":encode_sex(x["sex"]),"age_c":ac,"age_c2":ac**2},index=x.index)
        X=sm.add_constant(X,has_constant="add")
        r=sm.OLS(pd.to_numeric(x[phenotype],errors="coerce"),X).fit().resid
        sd=r.std(ddof=1)
        if np.isfinite(sd) and sd>0: out.loc[r.index]=(r-r.mean())/sd
    return out

def rho(x,y):
    z=pd.DataFrame({"x":x,"y":y}).dropna()
    if len(z)<3 or z["x"].nunique()<2 or z["y"].nunique()<2: return np.nan
    return float(spearmanr(z["x"],z["y"]).statistic)

def component_labels_from_edges(nodes,edges):
    g=nx.Graph()
    for r in edges.itertuples(index=False):
        g.add_edge(getattr(r,"row_1"),getattr(r,"row_2"))
    cmap={}
    for i,c in enumerate(nx.connected_components(g),1):
        for n in c: cmap[n]=f"SIB{i:04d}"
    return nodes["node_id"].map(cmap)

def family_means(nodes,groups,geno,pheno,min_members=1):
    z=pd.DataFrame({"group":groups,"population":nodes["population"],
                    "genotype":geno,"phenotype":pheno},index=nodes.index).dropna(subset=["group"])
    rows=[]
    for gid,g in z.groupby("group"):
        both=g.dropna(subset=["genotype","phenotype"])
        if len(both)<min_members: continue
        pops=both["population"].dropna().unique()
        if len(pops)!=1: raise ValueError(f"Family {gid} spans populations.")
        rows.append({"group":gid,"population":pops[0],"n_members":len(both),
                     "genotype_mean":both["genotype"].astype(float).mean(),
                     "phenotype_mean":both["phenotype"].astype(float).mean()})
    return pd.DataFrame(rows)

def perm_family(nodes,groups,locus,pheno,n_perm,seed):
    geno=pd.to_numeric(nodes[locus],errors="coerce")
    tab=family_means(nodes,groups,geno,pheno)
    obs=rho(tab.get("genotype_mean",pd.Series(dtype=float)),tab.get("phenotype_mean",pd.Series(dtype=float)))
    if not np.isfinite(obs): return {"rho":np.nan,"permutation_p":np.nan,"n_units":len(tab)}
    rng=np.random.default_rng(seed)
    valid=groups.notna() & geno.notna()
    popidx=[np.asarray(list(i),dtype=object) for i in nodes.loc[valid].groupby("population").groups.values()]
    null=np.empty(n_perm)
    for b in range(n_perm):
        gp=geno.copy()
        for labels in popidx: gp.loc[labels]=rng.permutation(geno.loc[labels].to_numpy())
        t=family_means(nodes,groups,gp,pheno)
        null[b]=rho(t["genotype_mean"],t["phenotype_mean"]) if len(t) else np.nan
    null=null[np.isfinite(null)]
    p=(1+np.sum(np.abs(null)>=abs(obs)))/(len(null)+1)
    return {"rho":obs,"permutation_p":float(p),"n_units":len(tab)}

def line_means(nodes,lines,geno,pheno):
    values=pd.DataFrame({"node_id":nodes["node_id"],"genotype":geno,"phenotype":pheno})
    z=lines.merge(values,left_on="offspring_id",right_on="node_id",how="left")
    rows=[]
    for lid,g in z.groupby("line_id"):
        both=g.dropna(subset=["genotype","phenotype"])
        if len(both)<1: continue
        rows.append({"group":lid,"population":both["population"].iloc[0],"n_members":len(both),
                     "genotype_mean":both["genotype"].astype(float).mean(),
                     "phenotype_mean":both["phenotype"].astype(float).mean()})
    return pd.DataFrame(rows)

def perm_lines(nodes,lines,locus,pheno,n_perm,seed):
    geno=pd.to_numeric(nodes[locus],errors="coerce")
    tab=line_means(nodes,lines,geno,pheno)
    obs=rho(tab.get("genotype_mean",pd.Series(dtype=float)),tab.get("phenotype_mean",pd.Series(dtype=float)))
    if not np.isfinite(obs): return {"rho":np.nan,"permutation_p":np.nan,"n_units":len(tab)}
    rng=np.random.default_rng(seed); eligible=set(lines["offspring_id"])
    valid=nodes["node_id"].isin(eligible) & geno.notna()
    popidx=[np.asarray(list(i),dtype=object) for i in nodes.loc[valid].groupby("population").groups.values()]
    null=np.empty(n_perm)
    for b in range(n_perm):
        gp=geno.copy()
        for labels in popidx: gp.loc[labels]=rng.permutation(geno.loc[labels].to_numpy())
        t=line_means(nodes,lines,gp,pheno)
        null[b]=rho(t["genotype_mean"],t["phenotype_mean"]) if len(t) else np.nan
    null=null[np.isfinite(null)]
    p=(1+np.sum(np.abs(null)>=abs(obs)))/(len(null)+1)
    return {"rho":obs,"permutation_p":float(p),"n_units":len(tab)}

def direct_parent(nodes,lines,locus,pheno,n_perm,seed):
    gmap=pd.Series(pd.to_numeric(nodes[locus],errors="coerce").to_numpy(),index=nodes["node_id"])
    pmap=pd.Series(pheno.to_numpy(),index=nodes["node_id"])
    rows=[]
    for lid,g in lines.groupby("line_id"):
        parent=g["parent_id"].iloc[0]; pop=g["population"].iloc[0]
        pg=gmap.get(parent,np.nan); ov=pmap.reindex(g["offspring_id"]).dropna()
        if pd.isna(pg) or not len(ov): continue
        rows.append({"line_id":lid,"parent_id":parent,"population":pop,
                     "parent_genotype":float(pg),"offspring_phenotype":float(ov.mean())})
    tab=pd.DataFrame(rows)
    obs=rho(tab.get("parent_genotype",pd.Series(dtype=float)),tab.get("offspring_phenotype",pd.Series(dtype=float)))
    if not np.isfinite(obs): return {"rho":np.nan,"permutation_p":np.nan,"n_units":len(tab)}
    rng=np.random.default_rng(seed); null=np.empty(n_perm)
    for b in range(n_perm):
        t=tab.copy()
        for _,idx in t.groupby("population").groups.items():
            labels=list(idx)
            t.loc[labels,"parent_genotype"]=rng.permutation(t.loc[labels,"parent_genotype"].to_numpy())
        null[b]=rho(t["parent_genotype"],t["offspring_phenotype"])
    p=(1+np.sum(np.abs(null)>=abs(obs)))/(n_perm+1)
    return {"rho":obs,"permutation_p":float(p),"n_units":len(tab)}

def run(nodes,phenotypes,loci,sibling_edges=None,mother_lines=None,father_lines=None,
        mixed_component_col=None,full_component_col="component_id",n_perm=10000,seed=20261115):
    rows=[]; counter=0
    sib=component_labels_from_edges(nodes,sibling_edges) if sibling_edges is not None else None
    scales=[]
    if sib is not None: scales.append(("sibling",sib))
    if mixed_component_col and mixed_component_col in nodes: scales.append(("mixed_generation",nodes[mixed_component_col]))
    if full_component_col in nodes: scales.append(("full_graph",nodes[full_component_col]))
    for phenotype in phenotypes:
        pheno=age_adjusted_residuals(nodes,phenotype)
        for scale,groups in scales:
            for locus in loci:
                s=perm_family(nodes,groups,locus,pheno,n_perm,seed+counter); counter+=1
                rows.append({"phenotype":phenotype,"scale":scale,"locus":locus,**s})
        for role,lines in (("mother",mother_lines),("father",father_lines)):
            if lines is None: continue
            for locus in loci:
                s=perm_lines(nodes,lines,locus,pheno,n_perm,seed+counter); counter+=1
                rows.append({"phenotype":phenotype,"scale":f"{role}_defined_offspring","locus":locus,**s})
                s=direct_parent(nodes,lines,locus,pheno,n_perm,seed+counter); counter+=1
                rows.append({"phenotype":phenotype,"scale":f"direct_{role}","locus":locus,**s})
    res=pd.DataFrame(rows); res["q9"]=np.nan
    for _,idx in res.groupby(["phenotype","scale"]).groups.items():
        res.loc[idx,"q9"]=bh_fdr(res.loc[idx,"permutation_p"])
    return res

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--nodes",required=True); p.add_argument("--sibling-edges",default=None)
    p.add_argument("--mother-lines",default=None); p.add_argument("--father-lines",default=None)
    p.add_argument("--mixed-component-col",default=None); p.add_argument("--full-component-col",default="component_id")
    p.add_argument("--phenotypes",nargs="+",default=["physical_aggression","verbal_aggression","anger","hostility"])
    p.add_argument("--loci",nargs="+",default=CORE9); p.add_argument("--n-perm",type=int,default=10000)
    p.add_argument("--seed",type=int,default=20261115)
    p.add_argument("--output",default="outputs/kinship_scale_decomposition.csv")
    a=p.parse_args()
    nodes=read_table(a.nodes)
    sibling=read_table(a.sibling_edges) if a.sibling_edges else None
    mother=read_table(a.mother_lines) if a.mother_lines else None
    father=read_table(a.father_lines) if a.father_lines else None
    res=run(nodes,a.phenotypes,a.loci,sibling,mother,father,a.mixed_component_col,
            a.full_component_col,a.n_perm,a.seed)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); res.to_csv(out,index=False)

if __name__=="__main__":
    main()
