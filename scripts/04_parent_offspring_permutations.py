"""Vertical parent-offspring phenotype correlations with permutation tests.

Each parent is one analytical unit. Offspring phenotype is averaged within the
directed parent cluster. The phenotype adjustment model is fit on the full
adult archive within population and then applied to parents/offspring.

Models:
- max_sample: phenotype ~ sex within population;
- age_adjusted: phenotype ~ sex + centered age + centered age^2 within
  population, using behavior-wave age >=18 only.

Permutation null: parent residual values/ranks are permuted among eligible
parents within population while parent identity, parent role, offspring
membership and offspring count remain fixed. The frozen analysis used 50,000
permutations.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

def read_table(path):
    p=Path(path)
    if p.suffix.lower() in {".xlsx",".xls"}: return pd.read_excel(p)
    return pd.read_csv(p,sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ",")

def bh_fdr(p):
    p=np.asarray(list(p),float); out=np.full(len(p),np.nan)
    ok=np.isfinite(p)
    if not ok.any(): return out
    vals=p[ok]; order=np.argsort(vals); ranked=vals[order]
    q=ranked*len(ranked)/np.arange(1,len(ranked)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    inv=np.empty_like(order); inv[order]=np.arange(len(order))
    out[np.where(ok)[0]]=q[inv]
    return out

def residualize_full_archive(df,trait,model):
    out=pd.Series(np.nan,index=df.index,dtype=float)
    for _,idx in df.groupby("population",dropna=False).groups.items():
        x=df.loc[idx].copy()
        x["age"]=pd.to_numeric(x["age"],errors="coerce")
        cols=[trait,"sex"]
        if model=="age_adjusted":
            x=x.loc[x["age"].ge(18)]
            cols.append("age")
        elif model!="max_sample":
            raise ValueError(model)
        x=x.dropna(subset=cols)
        if len(x)<3: continue
        sex=(x["sex"].astype(str).str.lower()=="male").astype(float)
        X=pd.DataFrame({"sex":sex},index=x.index)
        if model=="age_adjusted":
            ac=x["age"].astype(float)-x["age"].astype(float).mean()
            X["age_c"]=ac; X["age_c2"]=ac**2
        X=sm.add_constant(X,has_constant="add")
        r=sm.OLS(pd.to_numeric(x[trait],errors="coerce"),X).fit().resid
        sd=r.std(ddof=1)
        if np.isfinite(sd) and sd>0:
            out.loc[r.index]=(r-r.mean())/sd
    return out

def cluster_table(links,nodes,residuals,role):
    x=nodes.copy(); x["__resid"]=residuals
    nm=x.set_index("adult_id" if "adult_id" in x.columns else "node_id")
    rows=[]
    sub=links.loc[links["parent_role"].astype(str).str.lower().eq(role)]
    for pid,g in sub.groupby("parent_id"):
        if pid not in nm.index: continue
        pv=nm.at[pid,"__resid"]
        ov=nm.reindex(g["offspring_id"])["__resid"].dropna()
        if pd.isna(pv) or ov.empty: continue
        rows.append({
            "parent_id":pid,
            "population":g["population"].iloc[0],
            "parent_role":role,
            "parent_trait":float(pv),
            "offspring_trait":float(ov.mean()),
            "n_offspring":int(len(ov))
        })
    return pd.DataFrame(rows)

def rho(tab):
    if len(tab)<3: return np.nan
    return float(spearmanr(tab["parent_trait"],tab["offspring_trait"]).statistic)

def permutation_test(tab,n_perm=50000,seed=20261115):
    obs=rho(tab)
    if not np.isfinite(obs):
        return obs,np.nan
    rng=np.random.default_rng(seed); null=np.empty(n_perm)
    groups=[np.asarray(list(idx),dtype=int)
            for idx in tab.groupby("population").groups.values()]
    parent=tab["parent_trait"].to_numpy(float)
    offspring=tab["offspring_trait"].to_numpy(float)
    for b in range(n_perm):
        pp=parent.copy()
        for idx in groups:
            pp[idx]=rng.permutation(parent[idx])
        null[b]=spearmanr(pp,offspring).statistic
    p=float((1+np.sum(np.abs(null)>=abs(obs)))/(n_perm+1))
    return obs,p

def run(nodes,links,traits,model,n_perm=50000,seed=20261115):
    rows=[]; k=0
    for role in ("mother","father"):
        for trait in traits:
            resid=residualize_full_archive(nodes,trait,model)
            tab=cluster_table(links,nodes,resid,role)
            r,p=permutation_test(tab,n_perm,seed+k); k+=1
            rows.append({
                "model":model,"parent_role":role,"trait":trait,
                "n_parent_clusters":len(tab),"Spearman_rho":r,
                "permutation_p":p
            })
    res=pd.DataFrame(rows)
    res["BH_q"]=bh_fdr(res["permutation_p"])
    return res

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--links",required=True)
    p.add_argument("--nodes",required=True,
                   help="Full adult archive; IDs must match parent_id/offspring_id.")
    p.add_argument("--traits",nargs="+",
                   default=["physical_aggression","verbal_aggression","anger","hostility"])
    p.add_argument("--model",choices=["both","max_sample","age_adjusted"],default="both")
    p.add_argument("--n-perm",type=int,default=50000)
    p.add_argument("--seed",type=int,default=20261115)
    p.add_argument("--output",default="outputs/parent_offspring_permutations.csv")
    a=p.parse_args()
    nodes=read_table(a.nodes); links=read_table(a.links)
    models=["max_sample","age_adjusted"] if a.model=="both" else [a.model]
    frames=[run(nodes,links,a.traits,m,a.n_perm,a.seed+10000*i) for i,m in enumerate(models)]
    res=pd.concat(frames,ignore_index=True)
    # Frozen workbook applies BH-FDR within each model-specific 8-test family.
    res["BH_q"]=np.nan
    for _,idx in res.groupby("model").groups.items():
        res.loc[idx,"BH_q"]=bh_fdr(res.loc[idx,"permutation_p"])
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    res.to_csv(out,index=False)

if __name__=="__main__":
    main()
