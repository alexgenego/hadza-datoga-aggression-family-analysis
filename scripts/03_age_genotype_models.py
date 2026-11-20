"""Cross-age HTR1A/HTR2A models used in the manuscript.

Input is a de-identified harmonized child/adult table containing:
record_id, protocol, stage, population, sex, age,
physical_aggression_rank_normal_z, verbal_aggression_rank_normal_z,
HTR1A_G_dosage, HTR2A_harmonized_dosage.

Stage definitions for the primary strict comparison:
- child: child protocol;
- adult_gt18: adult protocol with behavioral age >18.
Archive adolescents and age-missing adults are excluded from the strict
child-vs-adult interaction but can be retained in the input for sensitivity
analyses.

Interaction models use age standardized within stage. Significance is assessed
with a Freedman-Lane residual permutation test (10,000 permutations), blocked
by stage × population. BH-FDR is applied within the predeclared four-test
family (2 loci × 2 outcomes).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

LOCI={
    "HTR1A":"HTR1A_G_dosage",
    "HTR2A":"HTR2A_harmonized_dosage"
}
OUTCOMES={
    "physical_aggression":"physical_aggression_rank_normal_z",
    "verbal_aggression":"verbal_aggression_rank_normal_z"
}

def read_table(path):
    p=Path(path)
    if p.suffix.lower() in {".xlsx",".xls"}: return pd.read_excel(p)
    return pd.read_csv(p,sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ",")

def bh_fdr(p):
    p=np.asarray(list(p),float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    vals=p[ok]; order=np.argsort(vals); ranked=vals[order]
    q=ranked*len(ranked)/np.arange(1,len(ranked)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    inv=np.empty_like(order); inv[order]=np.arange(len(order))
    out[np.where(ok)[0]]=q[inv]
    return out

def add_stage_age_z(x):
    x=x.copy()
    x["age"]=pd.to_numeric(x["age"],errors="coerce")
    x["age_within_stage_z"]=np.nan
    for _,idx in x.groupby("stage",dropna=False).groups.items():
        a=x.loc[idx,"age"].astype(float)
        sd=a.std(ddof=1)
        if np.isfinite(sd) and sd>0:
            x.loc[idx,"age_within_stage_z"]=(a-a.mean())/sd
    return x

def encode_model_columns(x,dosage,outcome):
    z=x.copy()
    z["sex_code"]=(z["sex"].astype(str).str.lower()=="female").astype(float)+1.0
    z["group"]=(z["stage"].astype(str)=="adult_gt18").astype(float)
    z["pop_datoga"]=(z["population"].astype(str)=="Datoga").astype(float)
    z["interaction"]=pd.to_numeric(z[dosage],errors="coerce")*z["group"]
    cols=[outcome,dosage,"group","pop_datoga","sex_code","age_within_stage_z","interaction"]
    z=z.dropna(subset=cols).copy()
    return z

def matrices(z,dosage,outcome,with_interaction=True):
    y=z[outcome].astype(float).to_numpy()
    base=pd.DataFrame({
        "const":1.0,
        "dosage":z[dosage].astype(float),
        "group":z["group"].astype(float),
        "population":z["pop_datoga"].astype(float),
        "sex":z["sex_code"].astype(float),
        "age_within_stage_z":z["age_within_stage_z"].astype(float),
    },index=z.index)
    if with_interaction:
        base["dosage_x_group"]=z["interaction"].astype(float)
    return y,base

def interaction_model(df,dosage,outcome,n_perm=10000,seed=20261115):
    x=df.loc[df["stage"].isin(["child","adult_gt18"])].copy()
    # The historical age-within-stage z score was calculated on the actual
    # locus × outcome analysis subset, after all required values were present.
    x=x.dropna(subset=[outcome,dosage,"population","sex","age","stage"]).copy()
    x=add_stage_age_z(x)
    z=encode_model_columns(x,dosage,outcome)
    y,Xfull=matrices(z,dosage,outcome,True)
    _,Xred=matrices(z,dosage,outcome,False)
    full=sm.OLS(y,Xfull).fit(); red=sm.OLS(y,Xred).fit()
    beta=float(full.params["dosage_x_group"])
    ci=full.conf_int().loc["dosage_x_group"].to_numpy(float)
    obs_f=float(full.compare_f_test(red)[0])

    fitted=np.asarray(red.fittedvalues,float)
    resid=np.asarray(red.resid,float)
    block_keys=(z["stage"].astype(str)+"|"+z["population"].astype(str))
    blocks=[np.asarray(list(idx),dtype=int)
            for idx in block_keys.reset_index(drop=True).groupby(block_keys.reset_index(drop=True)).groups.values()]
    rng=np.random.default_rng(seed); extreme=0
    # Matrix projection form avoids refitting formula encoders 10,000 times.
    xf=Xfull.to_numpy(float); xr=Xred.to_numpy(float)
    pinvf=np.linalg.pinv(xf); pinvr=np.linalg.pinv(xr)
    rankf=np.linalg.matrix_rank(xf); rankr=np.linalg.matrix_rank(xr)
    df_num=rankf-rankr; df_den=len(y)-rankf
    for _ in range(n_perm):
        rp=resid.copy()
        for idx in blocks:
            rp[idx]=rng.permutation(resid[idx])
        yp=fitted+rp
        br=pinvr@yp; bf=pinvf@yp
        ssr_r=float(np.sum((yp-xr@br)**2))
        ssr_f=float(np.sum((yp-xf@bf)**2))
        f=((ssr_r-ssr_f)/df_num)/(ssr_f/df_den)
        if f>=obs_f: extreme+=1
    p=float((extreme+1)/(n_perm+1))
    return {
        "N":len(z),"interaction_beta":beta,
        "CI_low":float(ci[0]),"CI_high":float(ci[1]),
        "partial_F":obs_f,"permutation_p":p,
        "child_N":int((z.stage=="child").sum()),
        "adult_gt18_N":int((z.stage=="adult_gt18").sum())
    }

def adult_age_model(df,dosage,outcome):
    x=df.loc[(df["protocol"].astype(str)=="adult") &
             pd.to_numeric(df["age"],errors="coerce").gt(18)].copy()
    x["age"]=pd.to_numeric(x["age"],errors="coerce")
    x["sex_code"]=(x["sex"].astype(str).str.lower()=="female").astype(float)+1.0
    x["pop_datoga"]=(x["population"].astype(str)=="Datoga").astype(float)
    need=[outcome,dosage,"age","sex_code","pop_datoga"]
    x=x.dropna(subset=need)
    age=(x["age"]-x["age"].mean())/x["age"].std(ddof=1)
    g=pd.to_numeric(x[dosage],errors="coerce")
    X=pd.DataFrame({
        "const":1.0,"dosage":g,"age_z":age,
        "population":x["pop_datoga"],"sex":x["sex_code"],
        "dosage_x_age":g*age
    },index=x.index)
    m=sm.OLS(x[outcome].astype(float),X.astype(float)).fit()
    return {
        "N":len(x),"interaction_beta":float(m.params["dosage_x_age"]),
        "parametric_p":float(m.pvalues["dosage_x_age"])
    }

def run(df,loci=("HTR1A","HTR2A"),outcomes=("physical_aggression","verbal_aggression"),
        n_perm=10000,seed=20261115):
    rows=[]; k=0
    for locus in loci:
        dosage=LOCI[locus]
        for outcome_name in outcomes:
            outcome=OUTCOMES[outcome_name]
            r=interaction_model(df,dosage,outcome,n_perm,seed+k); k+=1
            rows.append({
                "analysis":"strict_child_vs_adult_gt18",
                "locus":locus,"phenotype":outcome_name,**r
            })
    res=pd.DataFrame(rows)
    res["BH_q"]=bh_fdr(res["permutation_p"])

    adult=[]
    for locus in loci:
        for outcome_name in outcomes:
            adult.append({
                "analysis":"adult_genotype_x_age",
                "locus":locus,"phenotype":outcome_name,
                **adult_age_model(df,LOCI[locus],OUTCOMES[outcome_name])
            })
    return res,pd.DataFrame(adult)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data",required=True)
    p.add_argument("--n-perm",type=int,default=10000)
    p.add_argument("--seed",type=int,default=20261115)
    p.add_argument("--output",default="outputs/age_genotype_models.csv")
    p.add_argument("--adult-output",default="outputs/adult_genotype_age_models.csv")
    a=p.parse_args()
    primary,adult=run(read_table(a.data),n_perm=a.n_perm,seed=a.seed)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    primary.to_csv(out,index=False)
    adult.to_csv(Path(a.adult_output),index=False)

if __name__=="__main__":
    main()
