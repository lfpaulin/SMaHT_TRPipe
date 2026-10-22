#!/usr/bin/env python
import sys
import gzip
import polars as pl
import argparse


def load_gene_meta(filepath):
    data = {}
    if filepath.endswith("gz"):
        f = gzip.open(filepath, 'rt')
    else:
        f = open(filepath)
    for line in f:
        chrom, pos, end, gid, gname, gtype, *_ = line.rstrip("\n").split("\t")
        gid_raw = gid.split(".").pop(0)
        if gid_raw not in data:
            data[gid_raw] = (chrom, pos, end, gname, gtype, gid)
        else:
            print(2, gid_raw, gid)
    f.close()
    return data


def load_gene_exp(filepath, g_meta):
    known_colnames = ['chrom', 'pos', 'end', 'gene_id', 'gene_name']
    # known_colnames = ['chrom', 'pos', 'end', 'gene_id', 'gene_id_raw', 'gene_version', 'gene_name']
    # gene_exp_schema = {
    #     "chrom": pl.String,
    #     "gos": pl.UInt32,
    #     "gnd": pl.UInt32,
    #     "gene_id": pl.String,
    #     "gene_id_raw": pl.String, 
    #     "gene_version": pl.UInt8,
    #     "gene_name": pl.String
    # }
    gene_exp = pl.read_csv(filepath, separator='\t', has_header=True)

    gene_id_raw, gene_id_version = [], []
    for gid in gene_exp["gene_id"].to_list():
        g, v = gid.split(".")
        gene_id_raw.append(g) 
        gene_id_version.append(v)
    chroms, poss, ends, gnames = [], [], [], []
    # update with raw Gene ID
    for gid in gene_id_raw:
        if gid in g_meta:
            chrom, pos, end, gname, *_ = g_meta[gid]
            chroms.append(chrom)
            poss.append(int(pos))
            ends.append(int(end))
            gnames.append(gname)
        else:
            print(f'WARN: missing {gid}')

        # gene_id_raw = pl.Series(gene_id_raw),
        # gene_version = pl.Series(gene_id_version),
    gene_exp = gene_exp.with_columns( 
        chrom = pl.Series(chroms),
        pos = pl.Series(poss),
        end = pl.Series(ends),
        gene_name = pl.Series(gnames)
    )

    sample_cols = [coln for coln in gene_exp.columns if coln not in known_colnames]
    print(sample_cols)
    
    gene_exp_final: pl.DataFrame = pl.DataFrame({})
    do_once = True
    for each_sample in sample_cols:
        tmp = gene_exp[tuple(known_colnames)]
        donor, tissue, core, *_ = each_sample.split("_")
        tmp = tmp.with_columns(
            TPM = gene_exp[each_sample],
            Donor = pl.Series([donor] * gene_exp.height),
            Tissue = pl.Series([tissue] * gene_exp.height),
            Core = pl.Series([core] * gene_exp.height),
            SampleID = pl.Series([each_sample] * gene_exp.height)
        )
        if do_once:
            gene_exp_final = tmp
            do_once = False
        else:
            gene_exp_final = pl.concat([gene_exp_final, tmp])

    return gene_exp_final


def get_config():
    parser = argparse.ArgumentParser(description="Process transcript catalog and gene expression files.")
    parser.add_argument("--gene-exp", required=True, help="Path to the gene expression file (tab-separated)")
    parser.add_argument("--gene-meta", required=True, help="Path to the gene metadata fom gencode (parsed bed file)")
    parser.add_argument("--output", required=True, help="Path to the output file")
    return parser.parse_args()


def main():
    config = get_config()
    gene_meta = load_gene_meta(config.gene_meta)
    # print(len(gene_meta))
    gene_exp = load_gene_exp(config.gene_exp, gene_meta)
    gene_exp.write_csv(f'{config.output}.bed', separator="\t", include_header=False)
    gene_exp.write_parquet(f'{config.output}.pq')
    
    #trc = tr_catalog.join(
    #    other=gene_exp,
    #    on=("chrom", "pos", "end"),
    #    how="right",  # strategy: inner, left, right, full, semi, anti, cross
    #)
    #print(trc)
    

if __name__ == "__main__":
    main()
