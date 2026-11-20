"""Multilocus block, phenotype-dispersion and genetic-distance QC analyses.

The manuscript uses fixed adult complete-case panels before phenotype-specific
testing:
- core9: age >=18, both physical and verbal aggression available, all 9 loci
  complete (historical N=300; 293 Hadza, 7 Datoga);
- sensitivity7: same rule after dropping SERT_In2 and DRD4E3 (historical N=325;
  307 Hadza, 18 Datoga).

For the joint additive block, age is standardized within population and
phenotypes are standardized within population on the fixed panel. Covariates
include sex, age_z, age_z^2 and population. Genotype rows are shuffled within
population for the permutation null. The historical frozen workbook used
3,000 permutations; its original RNG state was not preserved.

Dispersion tests operate on covariate-adjusted standardized adult phenotypes.
Genetic-distance QC shuffles complete multilocus profiles within population
while keeping observed family pairs fixed.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import levene, fligner
import statsmodels.api as sm

CORE9=["DRD4P","SERTPr","SERT_In2","DRD4E3","DAT1","DRD2Taq1","COMT","HTR1A","HTR2A"]
SENS7=["DRD4P","SERTPr","DAT1","DRD2Taq1","COMT","HTR1A","HTR2A"]

def read_table(path):
    p=Path(path)
    if p.suffix.lower() in {".xlsx",".xls"}: return pd.read_excel(p)
    return pd.read_csv(p,sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ",")

def fixed_panel(df,loci):
    x=df.copy()
    x["age"]=pd.to_numeric(x["age"],errors="coerce")
    need=["age","sex","population","physical_aggression","verbal_aggression",*loci]
    x=x.loc[x["age"].ge(18)].dropna(subset=need).copy()
    for c in ["physical_aggression","verbal_aggression",*loci]:
        x[c]=pd.to_numeric(x[c],errors="coerce")
    return x.dropna(subset=["physical_aggression","verbal_aggression",*loci])

def prepare_panel(df,phenotype):
    x=df.copy()
    x["sex_code"]=(x["sex"].astype(str).str.lower()=="male").astype(float)
    x["population_code"]=(x["population"].astype(str)=="Datoga").astype(float)
    x["age_z"]=x.groupby("population")["age"].transform(
        lambda s:(s-s.mean())/s.std(ddof=1))
    x["age_z2"]=x["age_z"]**2
    x["outcome_z"]=x.groupby("population")[phenotype].transform(
        lambda s:(s-s.mean())/s.std(ddof=1))
    return x

def base_matrix(x):
    return sm.add_constant(
        x[["sex_code","age_z","age_z2","population_code"]].astype(float),
        has_constant="add")

def full_matrix(x,loci):
    return sm.add_constant(
        x[["sex_code","age_z","age_z2","population_code",*loci]].astype(float),
        has_constant="add")

def partial_delta_r2(x,loci):
    y=x["outcome_z"].astype(float)
    m0=sm.OLS(y,base_matrix(x)).fit()
    m1=sm.OLS(y,full_matrix(x,loci)).fit()
    f,p,df_num=m1.compare_f_test(m0)
    return {
        "partial_delta_R2":float(m1.rsquared-m0.rsquared),
        "F":float(f),"parametric_p":float(p),"effective_df":int(round(df_num))
    }

def permute_genotype_rows(x,loci,rng):
    z=x.copy()
    for _,idx in z.groupby("population").groups.items():
        labels=np.asarray(list(idx))
        perm=rng.permutation(labels)
        z.loc[labels,loci]=x.loc[perm,loci].to_numpy()
    return z

def block_main(df,phenotype,loci,n_perm=3000,seed=20261115):
    panel=fixed_panel(df,loci)
    x=prepare_panel(panel,phenotype)
    obs=partial_delta_r2(x,loci)
    rng=np.random.default_rng(seed)
    null=np.empty(n_perm)
    for b in range(n_perm):
        z=permute_genotype_rows(x,loci,rng)
        null[b]=partial_delta_r2(z,loci)["partial_delta_R2"]
    obs["permutation_p"]=float((1+np.sum(null>=obs["partial_delta_R2"]))/(n_perm+1))
    obs.update({
        "phenotype":phenotype,"N":len(x),
        "Hadza_N":int((x.population=="Hadza").sum()),
        "Datoga_N":int((x.population=="Datoga").sum())
    })
    return obs

def interaction_test(df,phenotype,loci,moderator):
    panel=fixed_panel(df,loci)
    x=prepare_panel(panel,phenotype)
    y=x["outcome_z"].astype(float)
    X0=full_matrix(x,loci)
    X1=X0.copy()
    if moderator=="age":
        for loc in loci:
            X1[f"{loc}_X_age"]=x[loc].astype(float)*x["age_z"]
    elif moderator=="population":
        for loc in loci:
            X1[f"{loc}_X_population"]=x[loc].astype(float)*x["population_code"]
    else:
        raise ValueError(moderator)
    # Remove rank-redundant interaction columns explicitly.
    rank0=np.linalg.matrix_rank(X0.to_numpy(float))
    keep=list(X0.columns)
    current=X0.to_numpy(float)
    for col in [c for c in X1.columns if c not in X0.columns]:
        candidate=np.column_stack([current,X1[col].to_numpy(float)])
        if np.linalg.matrix_rank(candidate)>np.linalg.matrix_rank(current):
            keep.append(col); current=candidate
    X1=X1[keep]
    m0=sm.OLS(y,X0).fit(); m1=sm.OLS(y,X1).fit()
    f,p,df_num=m1.compare_f_test(m0)
    return {
        "phenotype":phenotype,"moderator":moderator,"N":len(x),
        "Hadza_N":int((x.population=="Hadza").sum()),
        "Datoga_N":int((x.population=="Datoga").sum()),
        "delta_R2":float(m1.rsquared-m0.rsquared),
        "F":float(f),"effective_df":int(round(df_num)),"parametric_p":float(p)
    }

def adjusted_outcome(df,phenotype):
    """Standardized residual used for adult genotype-dispersion summaries."""
    x=df.dropna(subset=[phenotype,"population","sex"]).copy()
    x["sex_code"]=(x["sex"].astype(str).str.lower()=="male").astype(float)
    # Full adult-protocol dispersion retains age-missing adults; primary
    # adjustment is population + sex. Age sensitivity is reported separately.
    x["y_z"]=x.groupby("population")[phenotype].transform(
        lambda s:(s-s.mean())/s.std(ddof=1))
    X=sm.add_constant(pd.get_dummies(
        x[["population","sex_code"]],columns=["population"],drop_first=True,dtype=float),
        has_constant="add")
    r=sm.OLS(x["y_z"].astype(float),X.astype(float)).fit().resid
    out=pd.Series(np.nan,index=df.index,dtype=float)
    out.loc[x.index]=r/r.std(ddof=1)
    return out

def dispersion_test(df,locus,phenotype,scope=None):
    x=df.copy()
    if scope is not None:
        x=x.loc[x["population"].astype(str).eq(scope)].copy()
    adj=adjusted_outcome(x,phenotype)
    z=pd.DataFrame({"genotype":x[locus+"_genotype"] if locus+"_genotype" in x else x[locus],
                    "residual":adj}).dropna()
    groups=[g["residual"].to_numpy(float) for _,g in z.groupby("genotype") if len(g)>=2]
    if len(groups)<2:
        return {"locus":locus,"phenotype":phenotype,"scope":scope or "Pooled",
                "N":len(z),"BF_F":np.nan,"BF_p":np.nan,"Fligner_p":np.nan}
    bf,bp=levene(*groups,center="median")
    fl,fp=fligner(*groups,center="median")
    return {"locus":locus,"phenotype":phenotype,"scope":scope or "Pooled",
            "N":len(z),"BF_F":float(bf),"BF_p":float(bp),"Fligner_p":float(fp)}

def genotype_distance(a,b):
    aa=np.asarray(a,dtype=float); bb=np.asarray(b,dtype=float)
    ok=np.isfinite(aa)&np.isfinite(bb)
    if not ok.any(): return np.nan
    return float(np.mean(np.abs(aa[ok]-bb[ok])))

def pair_distance_qc(nodes,pairs,loci,n_perm=10000,seed=20261115):
    complete=nodes.dropna(subset=["node_id","population",*loci]).copy()
    nm=complete.set_index("node_id")
    valid=[]
    for r in pairs.itertuples(index=False):
        relation=str(getattr(r,"relation",""))
        compatible=str(getattr(r,"genetic_parent_compatible","")).lower()
        if "parent_offspring" in relation and compatible=="false":
            continue
        if r.row_1 in nm.index and r.row_2 in nm.index:
            if nm.at[r.row_1,"population"]==nm.at[r.row_2,"population"]:
                valid.append((r.row_1,r.row_2,nm.at[r.row_1,"population"]))
    if not valid:
        return {"n_pairs":0,"distance_ratio":np.nan,"permutation_p":np.nan}
    obs=np.array([genotype_distance(nm.loc[a,loci],nm.loc[b,loci]) for a,b,_ in valid])
    obsmean=float(np.mean(obs)); rng=np.random.default_rng(seed); null=[]
    pops={p:complete.loc[complete.population.eq(p),"node_id"].tolist() for _,_,p in valid}
    profiles={p:nm.loc[ids,loci].to_numpy(float) for p,ids in pops.items()}
    for _ in range(n_perm):
        maps={}
        for p,ids in pops.items():
            maps[p]=dict(zip(ids,profiles[p][rng.permutation(len(ids))]))
        null.append(float(np.mean([
            genotype_distance(maps[p][a],maps[p][b]) for a,b,p in valid
        ])))
    null=np.asarray(null)
    return {"n_pairs":len(valid),"distance_ratio":float(obsmean/null.mean()),
            "permutation_p":float((1+np.sum(null<=obsmean))/(n_perm+1))}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data",required=True)
    p.add_argument("--panel",choices=["core9","sensitivity7"],default="core9")
    p.add_argument("--phenotypes",nargs="+",default=["physical_aggression","verbal_aggression"])
    p.add_argument("--dispersion-locus",default="HTR2A")
    p.add_argument("--dispersion-phenotype",default="verbal_aggression")
    p.add_argument("--pairs",default=None)
    p.add_argument("--n-perm",type=int,default=3000)
    p.add_argument("--outdir",default="outputs/multilocus_dispersion")
    a=p.parse_args()
    df=read_table(a.data); loci=CORE9 if a.panel=="core9" else SENS7
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    main_rows=[]; interaction_rows=[]
    for i,ph in enumerate(a.phenotypes):
        main_rows.append(block_main(df,ph,loci,a.n_perm,20261115+i))
        interaction_rows.append(interaction_test(df,ph,loci,"age"))
        interaction_rows.append(interaction_test(df,ph,loci,"population"))
    pd.DataFrame(main_rows).to_csv(out/f"{a.panel}_multilocus_block.csv",index=False)
    pd.DataFrame(interaction_rows).to_csv(out/f"{a.panel}_multilocus_interactions.csv",index=False)

    disp=[dispersion_test(df,a.dispersion_locus,a.dispersion_phenotype,s)
          for s in [None,"Hadza","Datoga"]]
    pd.DataFrame(disp).to_csv(out/"dispersion_tests.csv",index=False)

    if a.pairs:
        qc=pair_distance_qc(df,read_table(a.pairs),CORE9,10000)
        pd.DataFrame([qc]).to_csv(out/"genetic_distance_qc.csv",index=False)

if __name__=="__main__":
    main()
