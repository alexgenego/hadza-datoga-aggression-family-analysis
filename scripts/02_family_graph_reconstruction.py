"""Reconstruct de-identified family graphs and relation-specific subgraphs.

The public pipeline uses a de-identified edge registry. It preserves a strict
social/field graph for phenotype analyses and a nested Mendelian-compatible
genetic graph in which explicitly incompatible parent-offspring edges are
removed from the biological-parentage layer without deleting the social link.

Required columns: row_1, row_2, population, relation.
Optional: parent_role, sibling_role, genetic_parent_compatible, parent_status,
sibling_groups, parent_id, offspring_id.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import networkx as nx
import pandas as pd

PARENT = "parent_offspring"
FULL = "full_sib"
MAT = "maternal_sib"
PAT = "paternal_sib"

def read_table(path, sheet=None):
    path=Path(path)
    if path.suffix.lower() in {".xlsx",".xls"}:
        return pd.read_excel(path, sheet_name=sheet or 0)
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv",".txt"} else ",")

def truthy(v):
    if pd.isna(v): return None
    if isinstance(v,bool): return v
    s=str(v).strip().lower()
    if s in {"true","1","yes","y"}: return True
    if s in {"false","0","no","n"}: return False
    return None

def relations(v):
    if pd.isna(v): return set()
    return {x.strip() for x in str(v).split("+") if x.strip()}

def validate(edges):
    req={"row_1","row_2","population","relation"}
    miss=req-set(edges.columns)
    if miss: raise ValueError(f"Missing required columns: {sorted(miss)}")
    if edges[["row_1","row_2"]].isna().any().any():
        raise ValueError("row_1/row_2 may not contain missing IDs.")

def build_graph(edges, genetic=False):
    validate(edges)
    g=nx.Graph()
    for row in edges.itertuples(index=False):
        rel=relations(getattr(row,"relation"))
        if genetic and PARENT in rel and truthy(getattr(row,"genetic_parent_compatible",None)) is False:
            rel.remove(PARENT)
        if not rel: continue
        a,b=getattr(row,"row_1"),getattr(row,"row_2")
        pop=getattr(row,"population")
        g.add_node(a,population=pop); g.add_node(b,population=pop)
        g.add_edge(a,b,
            relation="+".join(sorted(rel)),
            parent_role=getattr(row,"parent_role",None),
            sibling_role=getattr(row,"sibling_role",None),
            genetic_parent_compatible=getattr(row,"genetic_parent_compatible",None),
            parent_status=getattr(row,"parent_status",None),
            sibling_groups=getattr(row,"sibling_groups",None))
    return g

def assign_components(g):
    rows=[]; by_pop={}
    for comp in nx.connected_components(g):
        pops={g.nodes[n].get("population") for n in comp}
        if len(pops)!=1: raise ValueError(f"Component spans populations: {pops}")
        by_pop.setdefault(next(iter(pops)),[]).append(set(comp))
    for pop, comps in sorted(by_pop.items()):
        prefix=str(pop)[0].upper() if str(pop) else "C"
        for i,comp in enumerate(sorted(comps,key=lambda c:(-len(c),min(map(str,c)))),1):
            cid=f"{prefix}{i:02d}"
            for node in sorted(comp,key=str):
                rows.append({"node_id":node,"population":pop,"component_id":cid})
    return pd.DataFrame(rows)

def edge_table(g, components):
    cmap=components.set_index("node_id")["component_id"].to_dict()
    return pd.DataFrame([
        {"component_id":cmap[a],"population":g.nodes[a].get("population"),
         "row_1":a,"row_2":b,**d}
        for a,b,d in g.edges(data=True)
    ])

def relation_subgraph(edges, kind):
    rs=edges["relation"].map(relations)
    if kind=="parent_offspring": mask=rs.map(lambda s:PARENT in s)
    elif kind=="full_sib": mask=rs.map(lambda s:FULL in s)
    elif kind=="shared_mother_only": mask=rs.map(lambda s:MAT in s and FULL not in s)
    elif kind=="shared_father_only": mask=rs.map(lambda s:PAT in s and FULL not in s)
    elif kind=="sibling": mask=rs.map(lambda s:bool(s & {FULL,MAT,PAT}))
    else: raise ValueError(kind)
    return edges.loc[mask].copy()

def parent_defined_offspring_lines(edges, role):
    if role not in {"mother","father"}: raise ValueError(role)
    req={"parent_id","offspring_id","parent_role","population"}
    miss=req-set(edges.columns)
    if miss:
        raise ValueError("Parent-defined lines require explicit parent_id/offspring_id; "
                         f"missing {sorted(miss)}. Direction is never guessed from legacy unordered pairs.")
    x=edges.copy()
    x=x.loc[x["relation"].map(lambda v:PARENT in relations(v))]
    x=x.loc[x["parent_role"].astype(str).str.lower().eq(role)]
    x=x.dropna(subset=["parent_id","offspring_id"])
    x=x[["population","parent_id","offspring_id","parent_role"]].drop_duplicates()
    x["line_id"]=x["population"].astype(str)+":"+role+":"+x["parent_id"].astype(str)
    return x[["line_id","population","parent_role","parent_id","offspring_id"]]

def write_outputs(edges, outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    social=build_graph(edges,False); genetic=build_graph(edges,True)
    sn=assign_components(social); gn=assign_components(genetic)
    se=edge_table(social,sn); ge=edge_table(genetic,gn)
    sn.to_csv(outdir/"social_graph_nodes.csv",index=False)
    se.to_csv(outdir/"social_graph_edges.csv",index=False)
    gn.to_csv(outdir/"genetic_graph_nodes.csv",index=False)
    ge.to_csv(outdir/"genetic_graph_edges.csv",index=False)
    for kind in ["parent_offspring","full_sib","shared_mother_only","shared_father_only","sibling"]:
        relation_subgraph(se,kind).to_csv(outdir/f"{kind}_edges.csv",index=False)
    if {"parent_id","offspring_id","parent_role"}.issubset(edges.columns):
        for role in ("mother","father"):
            parent_defined_offspring_lines(edges,role).to_csv(
                outdir/f"{role}_defined_offspring_lines.csv",index=False)
    pd.DataFrame([
        {"graph":"social","n_nodes":social.number_of_nodes(),"n_edges":social.number_of_edges(),
         "n_components":nx.number_connected_components(social)},
        {"graph":"genetic","n_nodes":genetic.number_of_nodes(),"n_edges":genetic.number_of_edges(),
         "n_components":nx.number_connected_components(genetic)}
    ]).to_csv(outdir/"graph_summary.csv",index=False)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--edges",required=True); p.add_argument("--sheet",default=None)
    p.add_argument("--outdir",default="outputs/family_graph")
    a=p.parse_args()
    write_outputs(read_table(a.edges,a.sheet),Path(a.outdir))

if __name__=="__main__":
    main()
