"""Grouped Random Forest with leave-one-family-component-out validation.

Manuscript design:
- family component is the validation unit, not a predictor;
- records without a family-component assignment are excluded before fitting;
- the max-sample cohort requires all four adult aggression traits (historical
  N=108/30 families); age-conditioned models additionally require age >=18
  (historical N=86/28 families);
- relatives never appear in training and test simultaneously;
- baseline and genetic models use identical family-blocked folds;
- each training family receives equal total weight (individual weight 1/n_family);
- genotype missing values are imputed with training-fold medians only;
- no missingness indicators are used;
- primary RF: 100 trees, min_samples_leaf=4, max_features=1.0;
- sensitivity checks use min_samples_leaf=2 and 6;
- performance: family-balanced MSE/RMSE and weighted R²;
- base vs genetic models are compared by a paired Monte Carlo sign-flip test
  on family-specific MSE differences.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

CORE9=["DRD4P","SERTPr","SERT_In2","DRD4E3","DAT1","DRD2Taq1","COMT","HTR1A","HTR2A"]

def read_table(path):
    p=Path(path)
    if p.suffix.lower() in {".xlsx",".xls"}: return pd.read_excel(p)
    return pd.read_csv(p,sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ",")

def family_weights(df,component_col):
    sizes=df.groupby(component_col)[component_col].transform("size").astype(float)
    return 1.0/sizes

def family_balanced_metrics(df,y,pred,component_col):
    tmp=pd.DataFrame({"family":df[component_col].to_numpy(),"y":y,"pred":pred})
    fam=[]
    for f,g in tmp.groupby("family"):
        mse=float(np.mean((g["y"]-g["pred"])**2))
        fam.append((f,mse,len(g),float(g["y"].mean()),float(g["pred"].mean())))
    ft=pd.DataFrame(fam,columns=["family","mse","n","y_mean","pred_mean"])
    fb_mse=float(ft["mse"].mean()); fb_rmse=float(np.sqrt(fb_mse))
    w=np.asarray(family_weights(df,component_col),float)
    y=np.asarray(y,float); pred=np.asarray(pred,float)
    ybar=np.average(y,weights=w)
    sse=np.sum(w*(y-pred)**2); sst=np.sum(w*(y-ybar)**2)
    r2=float(1-sse/sst) if sst>0 else np.nan
    return fb_mse,fb_rmse,r2,ft

def build_model(base_features,genetic_features,n_estimators,min_leaf,random_state):
    cats=[c for c in base_features if c in {"population","sex"}]
    nums=[c for c in base_features if c not in cats]
    transformers=[]
    if cats:
        transformers.append(("cat",OneHotEncoder(handle_unknown="ignore"),cats))
    if nums:
        transformers.append(("num","passthrough",nums))
    if genetic_features:
        transformers.append(("gen",Pipeline([
            ("impute",SimpleImputer(strategy="median"))
        ]),genetic_features))
    prep=ColumnTransformer(transformers,remainder="drop")
    rf=RandomForestRegressor(
        n_estimators=n_estimators,
        min_samples_leaf=min_leaf,
        max_features=1.0,
        random_state=random_state,
        n_jobs=-1
    )
    return Pipeline([("prep",prep),("rf",rf)])

def leave_one_family_out(df,phenotype,component_col,base_features,genetic_features,
                         n_estimators=100,min_leaf=4,seed=20261115):
    rows=[]; family_errors=[]
    families=list(pd.unique(df[component_col].dropna()))
    for k,fam in enumerate(families):
        test=df.loc[df[component_col].eq(fam)].copy()
        train=df.loc[~df[component_col].eq(fam)].copy()
        if train.empty or test.empty: continue
        ytr=train[phenotype].astype(float).to_numpy()
        yte=test[phenotype].astype(float).to_numpy()
        sw=np.asarray(family_weights(train,component_col),float)
        models={
            "baseline":build_model(base_features,[],n_estimators,min_leaf,seed+k),
            "genetic":build_model(base_features,genetic_features,n_estimators,min_leaf,seed+k)
        }
        preds={}
        for label,m in models.items():
            feats=base_features+(genetic_features if label=="genetic" else [])
            m.fit(train[feats],ytr,rf__sample_weight=sw)
            preds[label]=m.predict(test[feats])
            for node,y,p in zip(test.index,yte,preds[label]):
                rows.append({"family":fam,"row_index":node,"model":label,
                             "observed":float(y),"predicted":float(p),"min_samples_leaf":min_leaf})
        b=float(np.mean((yte-preds["baseline"])**2))
        g=float(np.mean((yte-preds["genetic"])**2))
        family_errors.append({"family":fam,"n_test":len(test),
                              "baseline_mse":b,"genetic_mse":g,"delta_mse_genetic_minus_base":g-b,
                              "min_samples_leaf":min_leaf})
    pred=pd.DataFrame(rows); ferr=pd.DataFrame(family_errors)
    return pred,ferr

def metrics_from_predictions(df,pred,component_col):
    rows=[]
    for model,g in pred.groupby("model"):
        idx=g["row_index"].to_numpy()
        sub=df.loc[idx]
        mse,rmse,r2,_=family_balanced_metrics(sub,g["observed"].to_numpy(),
                                              g["predicted"].to_numpy(),component_col)
        rows.append({"model":model,"family_balanced_MSE":mse,
                     "family_balanced_RMSE":rmse,"weighted_R2":r2})
    out=pd.DataFrame(rows)
    if set(out["model"])=={"baseline","genetic"}:
        b=float(out.loc[out["model"].eq("baseline"),"weighted_R2"].iloc[0])
        g=float(out.loc[out["model"].eq("genetic"),"weighted_R2"].iloc[0])
        out["delta_R2_genetic_minus_base"]=np.where(out["model"].eq("genetic"),g-b,np.nan)
    return out

def sign_flip_test(differences,n_iter=100000,seed=20261115):
    d=np.asarray(differences,float); d=d[np.isfinite(d)]
    if len(d)==0: return np.nan
    obs=float(np.mean(d))
    rng=np.random.default_rng(seed)
    extreme=0
    # two-sided paired Monte Carlo sign-flip test
    for _ in range(n_iter):
        stat=float(np.mean(d*rng.choice([-1.0,1.0],size=len(d))))
        if abs(stat)>=abs(obs): extreme+=1
    return float((extreme+1)/(n_iter+1))

def run(df,phenotypes,component_col,base_features,genetic_features,
        min_leaves=(4,2,6),n_estimators=100,n_signflip=100000,seed=20261115,
        outdir="outputs/random_forest",cohort_name="max_sample"):
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    summaries=[]
    # Family-blocked validation is defined only for adults assigned to the
    # reconstructed graph. Non-graph archive records must never enter training.
    graph=df.loc[df[component_col].notna()].copy()
    # Manuscript cohorts are fixed before outcome-specific fitting.
    traits=["physical_aggression","verbal_aggression","anger","hostility"]
    graph=graph.dropna(subset=traits+[component_col,"population","sex"]).copy()
    if "age" in base_features:
        graph["age"]=pd.to_numeric(graph["age"],errors="coerce")
        graph=graph.loc[graph["age"].ge(18)].dropna(subset=["age"]).copy()
    for pi,pheno in enumerate(phenotypes):
        need=[pheno,component_col,*base_features]
        x=graph.dropna(subset=need).copy()
        x[pheno]=pd.to_numeric(x[pheno],errors="coerce")
        x=x.dropna(subset=[pheno])
        for li,leaf in enumerate(min_leaves):
            pred,ferr=leave_one_family_out(x,pheno,component_col,base_features,genetic_features,
                                           n_estimators,leaf,seed+pi*10000+li*1000)
            met=metrics_from_predictions(x,pred,component_col)
            p=sign_flip_test(ferr["delta_mse_genetic_minus_base"],n_signflip,
                             seed+pi*20000+li*2000)
            met.insert(0,"cohort",cohort_name)
            met.insert(1,"phenotype",pheno); met.insert(2,"N",len(x))
            met.insert(3,"n_families",x[component_col].nunique()); met.insert(4,"min_samples_leaf",leaf)
            met["paired_signflip_p"]=p
            summaries.append(met)
            stem=f"{cohort_name}_{pheno}_leaf{leaf}"
            pred.to_csv(out/f"{stem}_predictions.csv",index=False)
            ferr.to_csv(out/f"{stem}_family_errors.csv",index=False)
    return pd.concat(summaries,ignore_index=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data",required=True)
    p.add_argument("--phenotypes",nargs="+",default=["physical_aggression","verbal_aggression","anger","hostility"])
    p.add_argument("--component-col",default="component_id")
    p.add_argument("--mode",choices=["both","max_sample","age_known"],default="both")
    p.add_argument("--genetic-features",nargs="+",default=CORE9)
    p.add_argument("--min-leaves",nargs="+",type=int,default=[4,2,6])
    p.add_argument("--n-estimators",type=int,default=100)
    p.add_argument("--n-signflip",type=int,default=100000)
    p.add_argument("--seed",type=int,default=20261115)
    p.add_argument("--outdir",default="outputs/random_forest")
    a=p.parse_args()
    df=read_table(a.data); frames=[]
    if a.mode in {"both","max_sample"}:
        frames.append(run(df,a.phenotypes,a.component_col,["population","sex"],a.genetic_features,
                          tuple(a.min_leaves),a.n_estimators,a.n_signflip,a.seed,
                          a.outdir,"max_sample"))
    if a.mode in {"both","age_known"}:
        frames.append(run(df,a.phenotypes,a.component_col,["population","sex","age"],a.genetic_features,
                          tuple(a.min_leaves),a.n_estimators,a.n_signflip,a.seed+50000,
                          a.outdir,"age_known"))
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    pd.concat(frames,ignore_index=True).to_csv(out/"rf_summary.csv",index=False)

if __name__=="__main__":
    main()
