"""Permutation tests for family-graph phenotype clustering.

IMPORTANT: residual models are fit on the full adult analysis-ready archive,
not only on adults that belong to the reconstructed graph. Graph membership is
then used only when calculating the observed family statistic. The permutation
null draws pseudo-family residuals from the full eligible population pool while
preserving population and observed family/component sizes.

Statistic: for every connected component with >=2 usable graph adults,
calculate the mean pairwise squared distance between standardized residuals;
components receive equal weight.

Max-sample residuals: phenotype ~ sex within population.
Age-adjusted residuals: phenotype ~ sex + centered age + centered age^2 within
population, behavior-wave age >=18 only.

The historical frozen workbook used 30,000 permutations. A fixed public seed
is supplied for deterministic reconstruction; historical RNG state was not
preserved.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

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
    vals=sorted(s.dropna().astype(str).unique())
    if len(vals)>2: raise ValueError(f"sex has >2 categories: {vals}")
    return s.astype("string").map({v:i for i,v in enumerate(vals)}).astype(float)

def residualize_trait(df,trait,model,population_col="population",sex_col="sex",age_col="age"):
    """Fit adjustment models on every eligible archive record in each population."""
    out=pd.Series(np.nan,index=df.index,dtype=float)
    for _,idx in df.groupby(population_col,dropna=False).groups.items():
        x=df.loc[idx].copy()
        cols=[trait,sex_col]
        if model=="age_adjusted":
            x[age_col]=pd.to_numeric(x[age_col],errors="coerce")
            x=x.loc[x[age_col].ge(18)]; cols.append(age_col)
        elif model!="max_sample":
            raise ValueError(model)
        x=x.dropna(subset=cols)
        if len(x)<3: continue
        X=pd.DataFrame({"sex":encode_sex(x[sex_col])},index=x.index)
        if model=="age_adjusted":
            age=x[age_col].astype(float); ac=age-age.mean()
            X["age_c"]=ac; X["age_c2"]=ac**2
        X=sm.add_constant(X,has_constant="add")
        y=pd.to_numeric(x[trait],errors="coerce").astype(float)
        r=sm.OLS(y,X,missing="drop").fit().resid
        sd=r.std(ddof=1)
        if np.isfinite(sd) and sd>0:
            out.loc[r.index]=(r-r.mean())/sd
    return out

def infer_components_from_edges(nodes,edges,node_col="node_id"):
    import networkx as nx
    g=nx.Graph()
    for row in edges.itertuples(index=False):
        g.add_edge(getattr(row,"row_1"),getattr(row,"row_2"))
    cmap={}
    for i,comp in enumerate(nx.connected_components(g),1):
        for n in comp: cmap[n]=f"R{i:04d}"
    return nodes[node_col].map(cmap)

def _component_arrays(df,residuals,component,scope,population_col):
    graph_mask=residuals.notna() & component.notna()
    if scope!="Pooled":
        graph_mask &= df[population_col].astype(str).eq(scope)
    arrays=[]
    for _,idx in component.loc[graph_mask].groupby(component.loc[graph_mask]).groups.items():
        a=np.asarray(list(idx),dtype=int)
        if len(a)>=2: arrays.append(a)
    return arrays,graph_mask

def _stat(values,arrays):
    if not arrays: return np.nan
    # Mean pairwise squared distance = 2 * sample variance.
    return float(np.mean([2.0*np.var(values[a],ddof=1) for a in arrays]))

def permutation_test(df,residuals,component,scope="Pooled",population_col="population",
                     n_perm=30000,seed=20261115):
    vals=residuals.to_numpy(float)
    arrays,graph_mask=_component_arrays(df,residuals,component,scope,population_col)
    obs=_stat(vals,arrays)
    if not np.isfinite(obs):
        return dict(observed_dispersion=np.nan,null_mean=np.nan,dispersion_ratio=np.nan,
                    permutation_p=np.nan,n_components=0,n_nodes=0)

    # Eligible reference pool = all archive adults with an adjusted residual.
    ref_mask=residuals.notna()
    if scope!="Pooled":
        ref_mask &= df[population_col].astype(str).eq(scope)

    pools={}
    graph_positions={}
    for pop,idx in df.loc[ref_mask].groupby(population_col).groups.items():
        pools[str(pop)]=np.asarray(list(idx),dtype=int)
    for pop,idx in df.loc[graph_mask].groupby(population_col).groups.items():
        graph_positions[str(pop)]=np.asarray(list(idx),dtype=int)

    rng=np.random.default_rng(seed)
    null=np.empty(n_perm,dtype=float)
    work=vals.copy()
    for b in range(n_perm):
        rp=work.copy()
        for pop,gidx in graph_positions.items():
            pool=pools[pop]
            chosen=rng.choice(pool,size=len(gidx),replace=False)
            rp[gidx]=vals[chosen]
        null[b]=_stat(rp,arrays)

    null=null[np.isfinite(null)]
    p=float((1+np.sum(null<=obs))/(len(null)+1))
    return dict(
        observed_dispersion=obs,
        null_mean=float(np.mean(null)),
        dispersion_ratio=float(obs/np.mean(null)),
        permutation_p=p,
        n_components=len(arrays),
        n_nodes=int(graph_mask.sum())
    )

def run_screen(nodes,traits,model,component_col="component_id",relation_edges=None,
               scopes=("Pooled","Hadza","Datoga"),n_perm=30000,seed=20261115):
    comp=infer_components_from_edges(nodes,relation_edges) if relation_edges is not None else nodes[component_col]
    rows=[]
    # Residualize once on the full archive for each trait.
    residuals={trait:residualize_trait(nodes,trait,model) for trait in traits}
    for si,scope in enumerate(scopes):
        for j,trait in enumerate(traits):
            s=permutation_test(nodes,residuals[trait],comp,scope=scope,n_perm=n_perm,
                               seed=seed+si*1000+j)
            rows.append({"scope":scope,"model":model,"trait":trait,**s})
    res=pd.DataFrame(rows)
    res["BH_q"]=np.nan
    pooled=res["scope"].eq("Pooled")
    res.loc[pooled,"BH_q"]=bh_fdr(res.loc[pooled,"permutation_p"])
    res.loc[~pooled,"BH_q"]=bh_fdr(res.loc[~pooled,"permutation_p"])
    return res

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--nodes",required=True,
        help="Full adult analysis-ready archive; component_id is blank for non-graph adults.")
    p.add_argument("--traits",nargs="+",default=["physical_aggression","verbal_aggression","anger","hostility"])
    p.add_argument("--model",choices=["both","max_sample","age_adjusted"],default="both")
    p.add_argument("--component-col",default="component_id")
    p.add_argument("--relation-edges",default=None)
    p.add_argument("--n-perm",type=int,default=30000)
    p.add_argument("--seed",type=int,default=20261115)
    p.add_argument("--output",default="outputs/family_graph_permutations.csv")
    a=p.parse_args()
    nodes=read_table(a.nodes)
    edges=read_table(a.relation_edges) if a.relation_edges else None
    models=["max_sample","age_adjusted"] if a.model=="both" else [a.model]
    frames=[run_screen(nodes,a.traits,m,a.component_col,edges,n_perm=a.n_perm,
                       seed=a.seed+10000*i) for i,m in enumerate(models)]
    res=pd.concat(frames,ignore_index=True)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    res.to_csv(out,index=False)

if __name__=="__main__":
    main()
