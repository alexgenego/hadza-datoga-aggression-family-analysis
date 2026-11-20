"""Nested-CV LASSO and stability selection for adult multilocus panels.

Design follows the manuscript:
- covariates (age, sex, population as requested) are never penalized;
- phenotype and each genetic predictor are residualized using coefficients fit
  only on the relevant training fold;
- outer 5-fold CV is repeated 10 times;
- lambda_min and lambda_1se are chosen by an inner 5-fold CV;
- stability selection uses 500 population x sex stratified 80% subsamples;
- selection frequency >= 0.70 is the pre-specified stability threshold.

The complete-case cohort is fixed before phenotype-specific fitting: adults
age >=18 with both physical and verbal aggression available and complete panel
genotypes. The script does not silently impute missing genotypes for LASSO.
"""

from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

CORE9=["DRD4P","SERTPr","SERT_In2","DRD4E3","DAT1","DRD2Taq1","COMT","HTR1A","HTR2A"]
SENS7=["DRD4P","SERTPr","DAT1","DRD2Taq1","COMT","HTR1A","HTR2A"]
HADZA11=CORE9+["HTR1B","OXTR"]

def read_table(path):
    p=Path(path)
    if p.suffix.lower() in {".xlsx",".xls"}: return pd.read_excel(p)
    return pd.read_csv(p,sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ",")

def strata(df):
    return df["population"].astype(str)+"|"+df["sex"].astype(str)

def covariate_matrix(df,covariates,columns=None):
    x=df[list(covariates)].copy()
    cats=[c for c in covariates if c in {"population","sex"} or not pd.api.types.is_numeric_dtype(x[c])]
    x=pd.get_dummies(x,columns=cats,drop_first=True,dtype=float)
    x=x.astype(float)
    if columns is None: columns=list(x.columns)
    x=x.reindex(columns=columns,fill_value=0.0)
    return x,columns

@dataclass
class Residualizer:
    cov_columns:list[str]
    y_model:LinearRegression
    x_models:list[LinearRegression]
    x_scaler:StandardScaler

def fit_residualizer(train,phenotype,loci,covariates):
    C,cols=covariate_matrix(train,covariates)
    y=train[phenotype].astype(float).to_numpy()
    G=train[loci].astype(float).to_numpy()
    ym=LinearRegression().fit(C,y)
    yr=y-ym.predict(C)
    xmods=[]; gr=np.empty_like(G,dtype=float)
    for j in range(G.shape[1]):
        m=LinearRegression().fit(C,G[:,j]); xmods.append(m)
        gr[:,j]=G[:,j]-m.predict(C)
    scaler=StandardScaler().fit(gr)
    return Residualizer(cols,ym,xmods,scaler),yr,scaler.transform(gr)

def apply_residualizer(obj,test,phenotype,loci,covariates):
    C,_=covariate_matrix(test,covariates,obj.cov_columns)
    y=test[phenotype].astype(float).to_numpy()
    G=test[loci].astype(float).to_numpy()
    yr=y-obj.y_model.predict(C)
    gr=np.empty_like(G,dtype=float)
    for j,m in enumerate(obj.x_models): gr[:,j]=G[:,j]-m.predict(C)
    return yr,obj.x_scaler.transform(gr)

def alpha_grid(X,y,n=100):
    # sklearn's Lasso objective: ||y-Xb||^2/(2n)+alpha*|b|_1
    amax=max(np.max(np.abs(X.T@y))/len(y),1e-6)
    return np.geomspace(amax,amax*1e-4,n)

def inner_select(train,phenotype,loci,covariates,n_splits=5,seed=1):
    skf=KFold(n_splits=n_splits,shuffle=True,random_state=seed)
    # one common grid based on full training residuals; fold scores still use fold-only residualization
    full_obj,yf,Xf=fit_residualizer(train,phenotype,loci,covariates)
    alphas=alpha_grid(Xf,yf)
    fold_mse=np.full((n_splits,len(alphas)),np.nan)
    for k,(tr,va) in enumerate(skf.split(train)):
        a=train.iloc[tr]; b=train.iloc[va]
        obj,yt,Xt=fit_residualizer(a,phenotype,loci,covariates)
        yv,Xv=apply_residualizer(obj,b,phenotype,loci,covariates)
        for j,alpha in enumerate(alphas):
            m=Lasso(alpha=float(alpha),max_iter=100000).fit(Xt,yt)
            fold_mse[k,j]=np.mean((yv-m.predict(Xv))**2)
    means=fold_mse.mean(0); se=fold_mse.std(0,ddof=1)/np.sqrt(n_splits)
    imin=int(np.argmin(means))
    threshold=means[imin]+se[imin]
    eligible=np.where(means<=threshold)[0]
    # alphas descend from strongest to weakest; choose strongest eligible
    i1=int(eligible[0])
    return float(alphas[imin]),float(alphas[i1])

def valid_complete_case(df,phenotype,loci,covariates,population=None,
                        required_phenotypes=("physical_aggression","verbal_aggression")):
    """Return the manuscript-defined strict adult panel.

    The panel is fixed before phenotype-specific fitting: age >=18, both
    physical and verbal aggression available, requested covariates available,
    and all loci in the panel complete. This reproduces N=300 for core9,
    N=325 for sensitivity7, and N=290 for the Hadza-only 11-locus panel.
    """
    x=df.copy()
    if population is not None:
        x=x.loc[x["population"].astype(str).eq(population)]
    if "age" in x:
        x["age"]=pd.to_numeric(x["age"],errors="coerce")
        x=x.loc[x["age"].ge(18)]
    need=[*required_phenotypes,*loci,*covariates]
    x=x.dropna(subset=list(dict.fromkeys(need))).copy()
    for col in [*required_phenotypes,*loci]:
        x[col]=pd.to_numeric(x[col],errors="coerce")
    return x.dropna(subset=[*required_phenotypes,*loci])

def outer_cv(df,phenotype,loci,covariates,repeats=10,n_splits=5,seed=20261115):
    records=[]
    for rep in range(repeats):
        skf=KFold(n_splits=n_splits,shuffle=True,random_state=seed+rep)
        for fold,(tr,te) in enumerate(skf.split(df),1):
            train=df.iloc[tr]; test=df.iloc[te]
            amin,a1=inner_select(train,phenotype,loci,covariates,n_splits,seed+1000*rep+fold)
            obj,yt,Xt=fit_residualizer(train,phenotype,loci,covariates)
            yv,Xv=apply_residualizer(obj,test,phenotype,loci,covariates)
            baseline_sse=float(np.sum(yv**2))
            for label,alpha in [("lambda_min",amin),("lambda_1se",a1)]:
                m=Lasso(alpha=alpha,max_iter=100000).fit(Xt,yt)
                pred=m.predict(Xv)
                sse=float(np.sum((yv-pred)**2))
                records.append({"repeat":rep+1,"fold":fold,"solution":label,
                    "alpha":alpha,"n_test":len(test),"baseline_sse":baseline_sse,
                    "model_sse":sse,"selected":int(np.sum(np.abs(m.coef_)>1e-12))})
    r=pd.DataFrame(records)
    summary=[]
    for sol,g in r.groupby("solution"):
        b=g["baseline_sse"].sum(); s=g["model_sse"].sum()
        summary.append({"solution":sol,"OOF_delta_R2":float(1-s/b),
                        "n_outer_predictions":int(g["n_test"].sum()),
                        "mean_selected":float(g["selected"].mean())})
    return r,pd.DataFrame(summary)

def stratified_subsample(df,fraction,rng):
    chosen=[]
    for _,idx in df.groupby(["population","sex"],dropna=False).groups.items():
        idx=np.asarray(list(idx))
        n=max(1,int(round(len(idx)*fraction)))
        chosen.extend(rng.choice(idx,size=min(n,len(idx)),replace=False).tolist())
    return df.loc[chosen].copy()

def stability_selection(df,phenotype,loci,covariates,n_iter=500,fraction=.8,seed=20261115):
    rng=np.random.default_rng(seed)
    counts_min=np.zeros(len(loci),int); counts_1se=np.zeros(len(loci),int)
    coef_min=[]; coef_1se=[]
    for b in range(n_iter):
        sub=stratified_subsample(df,fraction,rng)
        amin,a1=inner_select(sub,phenotype,loci,covariates,5,seed+b)
        obj,y,X=fit_residualizer(sub,phenotype,loci,covariates)
        row_min=[]; row_1se=[]
        for alpha,row,counts in [(amin,row_min,counts_min),(a1,row_1se,counts_1se)]:
            m=Lasso(alpha=alpha,max_iter=100000).fit(X,y)
            co=m.coef_
            counts += (np.abs(co)>1e-12)
            row.extend(co.tolist())
        coef_min.append(row_min); coef_1se.append(row_1se)
    return pd.DataFrame({
        "locus":loci,
        "selection_frequency_lambda_min":counts_min/n_iter,
        "selection_frequency_lambda_1se":counts_1se/n_iter,
        "mean_coefficient_lambda_min":np.mean(np.asarray(coef_min),axis=0),
        "mean_coefficient_lambda_1se":np.mean(np.asarray(coef_1se),axis=0),
        "stable_ge_0.70":counts_min/n_iter>=.70
    })

def full_fit(df,phenotype,loci,covariates,seed):
    amin,a1=inner_select(df,phenotype,loci,covariates,5,seed)
    obj,y,X=fit_residualizer(df,phenotype,loci,covariates)
    rows=[]
    for label,alpha in [("lambda_min",amin),("lambda_1se",a1)]:
        m=Lasso(alpha=alpha,max_iter=100000).fit(X,y)
        for locus,coef in zip(loci,m.coef_):
            rows.append({"solution":label,"alpha":alpha,"locus":locus,"coefficient":coef,
                         "selected":abs(coef)>1e-12})
    return pd.DataFrame(rows)

def run_panel(df,panel_name,loci,phenotypes,covariates,population=None,outdir="outputs/lasso",seed=20261115):
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    summaries=[]
    for i,pheno in enumerate(phenotypes):
        x=valid_complete_case(df,pheno,loci,covariates,population)
        folds,summ=outer_cv(x,pheno,loci,covariates,10,5,seed+i*10000)
        stab=stability_selection(x,pheno,loci,covariates,500,.8,seed+i*20000)
        fit=full_fit(x,pheno,loci,covariates,seed+i*30000)
        stem=f"{panel_name}_{pheno}"
        folds.to_csv(out/f"{stem}_outer_folds.csv",index=False)
        stab.to_csv(out/f"{stem}_stability.csv",index=False)
        fit.to_csv(out/f"{stem}_full_fit.csv",index=False)
        summ.insert(0,"phenotype",pheno); summ.insert(0,"N",len(x)); summ.insert(0,"panel",panel_name)
        summaries.append(summ)
    pd.concat(summaries,ignore_index=True).to_csv(out/f"{panel_name}_summary.csv",index=False)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data",required=True)
    p.add_argument("--panel",choices=["core9","sensitivity7","hadza11"],default="core9")
    p.add_argument("--phenotypes",nargs="+",default=["physical_aggression","verbal_aggression"])
    p.add_argument("--covariates",nargs="+",default=["age","sex","population"])
    p.add_argument("--outdir",default="outputs/lasso"); p.add_argument("--seed",type=int,default=20261115)
    a=p.parse_args()
    panels={"core9":(CORE9,None),"sensitivity7":(SENS7,None),"hadza11":(HADZA11,"Hadza")}
    loci,pop=panels[a.panel]
    run_panel(read_table(a.data),a.panel,loci,a.phenotypes,a.covariates,pop,a.outdir,a.seed)

if __name__=="__main__":
    main()
