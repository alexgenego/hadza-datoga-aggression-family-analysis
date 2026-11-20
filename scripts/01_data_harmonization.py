"""Harmonize de-identified child/adult phenotype and demographic inputs.

Rows and original missingness are preserved. In addition to standardized column
names, this module creates the cross-age phenotype transformation used in the
manuscript:

1. physical and verbal aggression are ranked within protocol × population;
2. tied observations receive average ranks;
3. inverse-normal scores use q=(rank-0.5)/N;
4. the resulting scores are standardized to mean 0 and sample SD 1 within the
   same protocol × population block.

Raw child and adult scores are therefore never pooled directly.
Model-specific filtering is deliberately deferred to downstream scripts.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

PHENOTYPES=["physical_aggression","verbal_aggression","anger","hostility"]

def read_table(path,sheet=None):
    p=Path(path)
    if p.suffix.lower() in {".xlsx",".xls"}:
        return pd.read_excel(p,sheet_name=sheet or 0)
    return pd.read_csv(p,sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ",")

def clean_name(x):
    s=re.sub(r"[^0-9A-Za-z_]+","_",str(x).strip()).strip("_").lower()
    aliases={
        "pha":"physical_aggression","verba":"verbal_aggression","ang":"anger",
        "host":"hostility","age":"age","sex":"sex","population":"population",
        "physical_aggression":"physical_aggression","verbal_aggression":"verbal_aggression"
    }
    return aliases.get(s,s)

def normalize_sex(s):
    m={"1":"male","2":"female","boy":"male","girl":"female",
       "m":"male","male":"male","f":"female","female":"female"}
    raw=s.astype("string").str.strip()
    return raw.str.lower().map(m).fillna(raw)

def normalize_population(s):
    x=s.astype("string").str.strip(); low=x.str.lower()
    return pd.Series(
        np.where(low.eq("hadza"),"Hadza",np.where(low.eq("datoga"),"Datoga",x)),
        index=s.index,dtype="string")

def inverse_normal_block(v):
    """Average-rank inverse normal transform, standardized to sample SD=1."""
    v=pd.to_numeric(v,errors="coerce")
    out=pd.Series(np.nan,index=v.index,dtype=float)
    ok=v.notna()
    if ok.sum()<2: return out
    r=v.loc[ok].rank(method="average")
    q=(r-0.5)/len(r)
    z=pd.Series(norm.ppf(q),index=r.index,dtype=float)
    sd=z.std(ddof=1)
    if np.isfinite(sd) and sd>0:
        out.loc[z.index]=(z-z.mean())/sd
    return out

def add_cross_age_transforms(df,protocol_col="protocol"):
    out=df.copy()
    if protocol_col not in out:
        out[protocol_col]="unspecified"
    if "population" not in out:
        return out
    for ph in ["physical_aggression","verbal_aggression"]:
        if ph not in out: continue
        z=pd.Series(np.nan,index=out.index,dtype=float)
        for _,idx in out.groupby([protocol_col,"population"],dropna=False).groups.items():
            z.loc[idx]=inverse_normal_block(out.loc[idx,ph])
        out[f"{ph}_rank_normal_z"]=z
    return out

def harmonize(df,dataset_label,protocol_col="protocol"):
    x=df.copy(); x.columns=[clean_name(c) for c in x.columns]
    x.insert(0,"dataset",dataset_label)
    if "population" in x: x["population"]=normalize_population(x["population"])
    if "sex" in x: x["sex"]=normalize_sex(x["sex"])
    for c in ["age",*PHENOTYPES]:
        if c in x: x[c]=pd.to_numeric(x[c],errors="coerce")
    x=add_cross_age_transforms(x,protocol_col)
    return x

def sample_flow(df):
    return pd.DataFrame([
        {"variable":c,"N_nonmissing":int(df[c].notna().sum()),
         "N_missing":int(df[c].isna().sum())}
        for c in df.columns
    ])

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",required=True); p.add_argument("--sheet",default=None)
    p.add_argument("--dataset-label",required=True)
    p.add_argument("--protocol-col",default="protocol")
    p.add_argument("--outdir",default="outputs/harmonized")
    a=p.parse_args()
    x=harmonize(read_table(a.input,a.sheet),a.dataset_label,a.protocol_col)
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    x.to_csv(out/f"{a.dataset_label}_analysis_ready.csv",index=False)
    sample_flow(x).to_csv(out/f"{a.dataset_label}_sample_flow.csv",index=False)

if __name__=="__main__":
    main()
