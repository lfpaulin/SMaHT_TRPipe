#!/usr/bin/env python
import os
import sys
from typing import Any
import numpy as np
import pandas as pd
import polars as pl
import argparse
from pysam import TabixFile as BED
from pysam import VariantFile as VCF
from pysam import VariantRecord as Var
from pysam import VariantRecordSamples as VCFSample


tissue_name_converter = {
    "3A": "Blood_Whole",
    "3B": "Buccal_Swab",
    "3C": "Esophagus",
    "3D": "Esophagus",
    "3E": "Colon_Asc",
    "3F": "Colon_Asc",
    "3G": "Colon_Desc",
    "3H": "Colon_Desc",
    "3I": "Liver_Sample",
    "3J": "Liver_Sample",
    "3K": "Adrena_GL",
    "3L": "Adrena_GL",
    "3M": "Adrena_GR",
    "3N": "Adrena_GR",
    "3O": "Aorta_Abd",
    "3P": "Aorta_Abd",
    "3Q": "Lung",
    "3R": "Lung",
    "3S": "Heart_LV",
    "3T": "Heart_LV",
    "3U": "Testis_L",
    "3V": "Testis_L",
    "3W": "Testis_R",
    "3X": "Testis_R",
    "3Y": "Ovary_L",
    "3Z": "Ovary_L",
    "3AA": "Ovary_R",
    "3AB": "Ovary_R",
    "3AC": "Fibroblast",
    "3AD": "Skin_Calf",
    "3AE": "Skin_Calf",
    "3AF": "Skin_Abdomen",
    "3AG": "Skin_Abdomen",
    "3AH": "Muscle",
    "3AI": "Muscle",
    "3AJ": "Brain",
    "3AK": "Front_Lobe_BrLH",
    "3AL": "Temp_Lobe_BrLH",
    "3AM": "Cereb_BrLH",
    "3AN": "Hippoc_BrLH",
    "3AO": "Hippoc_BrRH",
    "3AP": "Front_Lobe_BrLH",
    "3AQ": "Temp_Lobe_BrLH",
    "3AR": "Cereb_BrLH",
    "3AS": "Hippoc_BrLH",
    "3AT": "Hippoc_BrRH"
}


tissue_derms = {
    "Blood_Whole": "Mesoderm",
    "Buccal_Swab": "Other",
    "Esophagus": "Endoderm",
    "Esophagus": "Endoderm",
    "Colon_Asc": "Endoderm",
    "Colon_Asc": "Endoderm",
    "Colon_Desc": "Endoderm",
    "Colon_Desc": "Endoderm",
    "Liver_Sample": "Endoderm",
    "Liver_Sample": "Endoderm",
    "Adrena_GL": "Ectoderm",
    "Adrena_GL": "Ectoderm",
    "Adrena_GR": "Ectoderm",
    "Adrena_GR": "Ectoderm",
    "Aorta_Abd": "Mesoderm",
    "Aorta_Abd": "Mesoderm",
    "Lung": "Endoderm",
    "Lung": "Endoderm",
    "Heart_LV": "Mesoderm",
    "Heart_LV": "Mesoderm",
    "Testis_L": "Other",
    "Testis_L": "Other",
    "Testis_R": "Other",
    "Testis_R": "Other",
    "Ovary_L": "Other",
    "Ovary_L": "Other",
    "Ovary_R": "Other",
    "Ovary_R": "Other",
    "Fibroblast": "Other",
    "Skin_Calf": "Ectoderm",
    "Skin_Calf": "Ectoderm",
    "Skin_Abdomen": "Ectoderm",
    "Skin_Abdomen": "Ectoderm",
    "Muscle": "Mesoderm",
    "Muscle": "Mesoderm",
    "Brain": "Ectoderm",
    "Front_Lobe_BrLH": "Ectoderm",
    "Temp_Lobe_BrLH": "Ectoderm",
    "Cereb_BrLH": "Ectoderm",
    "Hippoc_BrLH": "Ectoderm",
    "Hippoc_BrRH": "Ectoderm",
    "Front_Lobe_BrLH": "Ectoderm",
    "Temp_Lobe_BrLH": "Ectoderm",
    "Cereb_BrLH": "Ectoderm",
    "Hippoc_BrLH": "Ectoderm",
    "Hippoc_BrRH": "Ectoderm"
}


smallify_gcc = {
    "broad": "BR",
    "uwsc": "UW",
    "bcm": "BC",
    "nysc": "NY",
    "wu": "WU"
}


seq_tech = {
    "B001": "PacBio Revio HiFi WG",
    "D001": "ONT PromethION 24 WG",
    "D003": "ONT PromethION 24 Ultra-Long WGS",
    "B004": "PacBio Revio HiFi Fiber-seq"
}


tissue_tech_converter = {
    "B001": "PacBio",
    "D001": "ONT",
    "D003": "ONT",
    "B004": "PacBio"
}


class TRConfig:
    def __init__(self, args: argparse.Namespace ):
        self.cmd = args.cmd
        self.donor = args.donor
        self.infile = args.infile
        self.tr_bed = args.tr_bed
        self.genes = args.genes

    def __repr__(self) -> str:
        return f'{self.cmd}  {self.infile}  <{self.donor}> <{self.tr_bed}> <{self.genes}>'


class VCFInfo:
    def __init__(self):
        self.info = "INFO"
        self.supp_vec =  [("ID", "SUPP_VEC"),    ("Number", 1), ("Type", "String"), ("Description", "Support vector")]
        self.supp_germ = [("ID", "SUPP_GERM"),   ("Number", 1), ("Type", "String"), ("Description", "Support vector for gemrline alleles")]
        self.supp_mosc = [("ID", "SUPP_MOSAIC"), ("Number", 1), ("Type", "String"), ("Description", "Support vector for mosaic alleles")]
        self.tr_annot =  [("ID", "TR_ANNOT"),    ("Number", 1), ("Type", "String"), ("Description", "TR instability/annotation")]
        self.tr_annot_add =  [("ID", "TR_ANNOT_INFO"),    ("Number", 1), ("Type", "String"), ("Description", "TR instability/annotation additional information")]
        self.tr_annot_lid =  [("ID", "TR_ANNOT_LOCUS"),    ("Number", 1), ("Type", "String"), ("Description", "TR annotation in other alleles of the same locus")]
        self.tr_annot_gene =  [("ID", "GENE"),    ("Number", 1), ("Type", "String"), ("Description", "TR annotation of gene overlaps")]


class MyLog():
    def __init__(self) -> None:
        pass

    @staticmethod
    def info(text: str):
        sys.stdout.write(f'INFO: {text}\n')

    @staticmethod
    def warn(text: str):
        sys.stdout.write(f'WARNING: {text}\n')

    @staticmethod
    def err(text: str):
        sys.stderr.write(f'ERROR: {text}\n')

    @staticmethod
    def dev(text: str, debug = False):
        if debug:
            sys.stderr.write(f'ERROR: {text}\n')

log = MyLog()


def parse_table(args: TRConfig):
    assert args.donor != "", f'Donor is needed, none provided: <{args.donor}>'
    locus_reads, donor = args.infile, args.donor
    need_close = False
    current_tissue = ""
    current_tissue_name = ""
    donor_alen = []
    donor_tissues = []
    counter = 0
    tissue = ""
    lid = ""
    if "-" == locus_reads:
        reader = sys.stdin
    else:
        assert os.path.isfile(locus_reads), f'File {locus_reads} does not exixst'
        reader = open(locus_reads)
        need_close = True
    print(f'donor <- "{donor}"\nalens <- list()')
    for line in reader:
        # file is 2 columns, col1=tissue code, col2=allele len
        # save all lengths in an R list
        line = line.rstrip("\n")
        tcode, lid, alen = line.split("\t")
        tissue = tissue_name_converter[tcode]
        if "" == current_tissue:
            #init
            current_tissue = tcode
            current_tissue_name = tissue
            donor_alen = [alen]
            donor_tissues = [f'"{current_tissue_name}"']
            counter += 1
        else:
            if tcode != current_tissue:
                alen_str = ", ".join(donor_alen)
                print(f'alens[["{current_tissue_name}"]] <- c({alen_str}) #{counter}')
                counter += 1
                # restart
                current_tissue = tcode
                current_tissue_name = tissue
                donor_tissues.append(f'"{tissue}"')
                donor_alen = [alen]
            else:
                donor_alen.append(alen)
    alen_str = ", ".join(donor_alen)
    print(f'alens[["{tissue}"]] <- c({alen_str})')
    donor_tissues_str = ", ".join(donor_tissues)
    print(f'donor_tissues <- c({donor_tissues_str})')
    print(f'lab_x <- "Tissue"; lab_y <- "Allele length"; lab_main <- "{donor}"; counter <- 1;')
    print(f'png("{donor}_extracted_locusID_{lid}.png", height=3000, width={len(donor_tissues)}*500, res=300)')
    print(f'plot(NA, ylim=c(0, max(unlist(alens))), xlim=c(0,length(donor_tissues)), xaxt ="n", xlab=lab_x, ylab=lab_y, main=lab_main, cex.lab=1.5, cex.main=2); \n'
          f'for (tis in donor_tissues) {{ \n'
          f'    y <- alens[[tis]];  \n'
          f'    x <- runif(length(y), counter-0.25, counter+0.25); \n'
          f'    r <- runif(1, 0, 255); g <- runif(1, 0, 255); b <- runif(1, 0, 255); \n'
          f'    points(x, y, pch=19, cex=2, col=rgb(r, g, b, 128, maxColorValue=255)); \n'
          f'    counter = counter+1  }}; \n'
          f'axis(side = 1, at = 1:length(donor_tissues), labels = donor_tissues)\n'
          f'dev.off()\n'
        )

    if need_close:
        reader.close()


def name2tissue(sample_name: str):
    # expected format: DONOR-TISSUE-CORE-SEX/AGE-TECH-GCC
    # example SMHT001-3A-XX-M42-B001-broad
    return tissue_name_converter[sample_name.split("-")[1]]


def name2tech(sample_name: str):
    # expected format: DONOR-TISSUE-CORE-SEX/AGE-TECH-GCC
    # example SMHT001-3A-XX-M42-B001-broad
    return tissue_tech_converter[sample_name.split("-")[4]]


def process_locus(lid:str, locus: pl.DataFrame, vaf_state: str, debug: bool = False):
    """
    lid: LocusID
    locus: polars DataFrame that has both ALT and REF-like counts. We count as ALT the 
           allele of interest and REF the remaning reads. If more than one allele exisits
           REF is the sum of the read counts of the remaining alleles.
    vaf_state: Filter applied to the locus, either GERMLINE or PASS, the second being 
               like mosaic, unstable or error
    """
    # INFO tags to be added to the VCF entry
    infos = {}
    # For the output update
    allele_count = 0
    # Overall coverage => ALT + REF per tissue/core/sample
    coverage = locus.select(pl.all().sum()).to_numpy()
    # Number of tissue/core/samples
    n_samples = locus.width
    # VAF-like value computes as ALT/(REF+ALT), currently 20% delimits the 
    # difference between GERMLINE and MOSAIC. This value is arbitrary
    vaf_limit_germline_mosaic = 0.20
    locus_infos = []
    for row in locus.iter_rows():
        # counts and vaf-like 
        np_row = np.array(row)
        row_p = np.array([ r/c if c > 0 else 0 for r,c in zip(row, coverage.tolist()[0]) ])
        # support vector: counted as present if sample has +2 reads or, overall reads in locus is +3
        # will be added as INFO filed 'SUPP_VEC'
        supp_vec = [1 if (x > 1 or (x == 1 and sum(row) >= 2)) else 0 for x in row]
        # sample names according to the support vector
        in_samples = [smpl for smpl, supp in zip(locus.columns, supp_vec) if 1 == supp]
        # number of NON-ZERO entries
        n_suppv_tissue = np.nansum([1 if x > 0 else 0 for x in row])
        # samples with mosaic-like counts, cannot be zero and should be at least two reads support
        # will be added as INFO filed 'SUPP_MOSAIC'
        supp_mosaic = [1 if x else 0 for x in ((row_p < vaf_limit_germline_mosaic) & (np_row > 1)).tolist()]
        # samples with germline-like counts, will be added as INFO filed 'SUPP_GERMLINE'
        supp_germ = [1 if x else 0 for x in (row_p >= vaf_limit_germline_mosaic).tolist()]
        # will be added as INFO filed 'TR_ANNOT'
        tr_annot = ""
        tr_annot_add = ""
        if np.nansum(supp_vec) == 0:
            # none
            tr_annot = "LOW_QC"
        elif n_suppv_tissue == 1:
            # core
            if sum(row) == 1:
                tr_annot = "LOW_QC"
            else:
                tiss = [name2tissue(t) for t in in_samples].pop()
                tr_annot = f'CORE'
                tr_annot_add = f'{tiss}'
            log.dev(f'{allele_count}|core', debug)
        elif n_suppv_tissue == n_samples:
            # in all
            and_germ = ""
            if sum(supp_germ) == n_samples:
                and_germ = "GERMLINE,"
            tr_annot = f'{and_germ}ALL'
            log.dev(f'{allele_count}|in_all', debug)
        elif n_suppv_tissue >= n_samples * 0.80:
            # almost all
            and_germ = ""
            if sum(supp_germ) == n_suppv_tissue:
                and_germ = "GERMLINE,"
            tr_annot = f'{and_germ}UBIQ'
            log.dev(f'{allele_count}|the almost: {n_suppv_tissue} {n_samples}', debug)
        elif 1 < n_suppv_tissue <= 3:
            # maybe tissue
            tiss = set([name2tissue(t) for t in in_samples])
            from_tech = [name2tech(t) for t in in_samples]        
            tr_annot = f'TISSUE'
            if len(tiss) == 1:
                # tissue
                tr_annot_add = f'{tiss.pop()}'
            else:
                # fibroblast++
                tiss_str = ",".join([name2tissue(t) for t in in_samples])
                if "Fibroblast,Skin,Calf" == tiss_str:  
                    tr_annot_add = f'{tiss_str}'
                elif "Blood" in tiss_str and len(tiss) == 2:
                    tr_annot_add = f'{tiss_str}'
                elif len(set(from_tech)) == 1:
                    from_tech_str = set(from_tech).pop()
                    tr_annot = f'TECH_BIAS'
                    tr_annot_add = f'{from_tech_str}'
                else:
                    # maybe derms
                    derm = [tissue_derms[name2tissue(t)] for t in in_samples]
                    derm_str = ",".join(derm)
                    uderm = set(derm)
                    tr_annot = f'DERM'
                    if len(uderm) == 1 or (len(uderm) == 2 and "Fibroblast" in tiss) or (len(uderm) == 2 and "Blood,Whole" in tiss):
                        if len(uderm) == 2 and "Fibroblast" in tiss:
                            tr_annot_add = f'{[x for x in uderm if x != "Other"].pop()}'
                        elif len(uderm) == 2 and "Blood,Whole" in tiss:
                            tr_annot_add = f'{[x for x in uderm if x != "Mesoderm"].pop()}'
                        else:
                            tr_annot_add = f'{uderm.pop()}'
                    else:
                        this_vaf = "MOSAIC,"
                        if sum(supp_germ) > sum(supp_mosaic):
                            this_vaf = "GERMLINE,"
                        tr_annot = f'{this_vaf}UBIQ'
            log.dev(f'{allele_count}|maybe tissue: {n_suppv_tissue} {n_samples} {tiss}', debug)
        elif n_suppv_tissue > 3:
            # maybe derm
            tiss = [name2tissue(t) for t in in_samples]
            derm = [tissue_derms[name2tissue(t)] for t in in_samples]
            from_tech = [name2tech(t) for t in in_samples]
            uderm = set(derm)
            derm_str = ",".join(derm)
            tiss_str = ",".join(tiss)
            tr_annot = f'DERM'
            if len(set(from_tech)) == 1:
                from_tech_str = set(from_tech).pop()
                tr_annot = f'TECH_BIAS'
                tr_annot_add = f'{from_tech_str}'
            elif len(uderm) == 1 or (len(uderm) == 2 and "Fibroblast" in tiss) or (len(uderm) == 2 and "Blood,Whole" in tiss):
                if len(uderm) == 2 and "Fibroblast" in tiss:
                    tr_annot_add = f'{[x for x in uderm if x != "Other"].pop()}'
                elif len(uderm) == 2 and "Blood,Whole" in tiss:
                    tr_annot_add = f'{[x for x in uderm if x != "Mesoderm"].pop()}'
                else:
                    tr_annot_add = f'{uderm.pop()}'
            else:
                # TODO :AQUI
                if "GERMLINE" == vaf_state:
                    tr_annot = "GERMLINE,UBIQ"
                else:
                    this_vaf = "MOSAIC,"
                    if sum(supp_germ) > sum(supp_mosaic):
                        this_vaf = "GERMLINE,"
                    tr_annot = f'{this_vaf}UBIQ'
            log.dev(f'{allele_count}|maybe derm: {n_suppv_tissue} {n_samples} {derm_str}|{tiss_str}', debug)
        else:
            # ??
            tr_annot = "UNKNOWN"
            log.warn(f'{allele_count}|?')
        supp_vec_str = "".join([str(s) for s in supp_vec])
        supp_mosaic_str = "".join([str(s) for s in supp_mosaic])
        supp_germ_str = "".join([str(s) for s in supp_germ])
        if sum(supp_vec) > 0 and sum(supp_mosaic) == 0 and "GERMLINE" not in tr_annot:
            tr_annot = f'LOW_SUPP,{tr_annot}'
        tr_annot = tr_annot if "" != tr_annot else None
        tr_annot_add = tr_annot_add if "" != tr_annot_add else None
        infos[f'{lid}.{allele_count}'] = [
            ("SUPP_VEC", supp_vec_str),
            ("SUPP_GERM", supp_germ_str),
            ("SUPP_MOSAIC", supp_mosaic_str),
            ("TR_ANNOT", tr_annot),
            ("TR_ANNOT_INFO", tr_annot_add)
        ]
        if tr_annot is not None and "TECH_BIAS" in tr_annot:
            if "" != tr_annot_add:
                locus_infos.append((f'{lid}.{allele_count}', f'({allele_count},{tr_annot},{tr_annot_add})'))
            else:
                locus_infos.append((f'{lid}.{allele_count}', f'({allele_count},{tr_annot})'))
        allele_count += 1

    for ac in range(0, allele_count, 1):
        add_info = ",".join([txt for aid, txt in locus_infos if aid != f'{lid}.{ac}'])
        add_info = add_info if "" != add_info else None
        infos[f'{lid}.{ac}'].append(("TR_ANNOT_LOCUS", add_info))
    return infos


def make_df(vcf_tr: VCFSample):
    df_schema = {tissue.name: pl.UInt16 for tissue in vcf_tr.values()}
    return pl.DataFrame({tissue.name: tissue.get("AD")[0] for tissue in vcf_tr.values()})


def make_support(vcf_tr: VCFSample):
    df_schema = {tissue.name: pl.UInt8 for tissue in vcf_tr.values()}
    return pl.DataFrame({tissue.name: 1 if tissue.get("AD")[0] > 0 else 0 for tissue in vcf_tr.values()})


def load_catalog(bed_file: str):
    bed = BED(bed_file)
    catalog_dict = {}
    for entry in bed.fetch():
        contig, pos, end, info = entry.split("\t")
        lid = info.split(";").pop(0).split("=").pop()
        catalog_dict[lid] = (contig, pos, end)
    bed.close()
    return catalog_dict


def get_gene(tr_id, genes, catalog):
    lid = tr_id.split(".").pop(0)
    chro, pos, end = catalog[lid]
    overlap_genes = []
    for entry in genes.fetch(region=f'{chro}:{pos}-{end}'):
        _, _, _, gene_info = entry.split("\t")
        g = [item.split('=').pop() for item in gene_info.split(';') if '=' in item and ("gene_name" in item or "gene_type" in item) and "pseudogene" not in gene_info]
        if len(g) > 0:
            gg = ",".join(g)
            overlap_genes.append(gg)
    return ",".join(overlap_genes)


def parse_vcf(config: TRConfig):
    assert os.path.isfile(config.tr_bed), f'File {config.tr_bed} does not exixist'
    assert os.path.isfile(config.genes), f'File {config.genes} does not exixist'
    vcf_file = config.infile
    vcf = VCF(vcf_file)
    vcf_info = VCFInfo()
    vcf.header.add_meta(key=vcf_info.info, items=vcf_info.supp_vec)
    vcf.header.add_meta(key=vcf_info.info, items=vcf_info.supp_germ)
    vcf.header.add_meta(key=vcf_info.info, items=vcf_info.supp_mosc)
    vcf.header.add_meta(key=vcf_info.info, items=vcf_info.tr_annot)
    vcf.header.add_meta(key=vcf_info.info, items=vcf_info.tr_annot_add)
    vcf.header.add_meta(key=vcf_info.info, items=vcf_info.tr_annot_lid)
    vcf.header.add_meta(key=vcf_info.info, items=vcf_info.tr_annot_gene)
    vcf_split_name = vcf_file.split(".")
    use_index_name = len(vcf_split_name) -1 if "gz" not in vcf_file else len(vcf_split_name) -2
    vcf_base_name = ".".join(vcf_split_name[:use_index_name])
    genes = BED(config.genes)
    catalog = load_catalog(config.tr_bed)
    out = VCF(f'{vcf_base_name}.annot.vcf', "w", header=vcf.header)
    tr: Var
    locus = pl.DataFrame()
    locus_id = ""
    current_lid = ""
    nloci = 0
    counter = 0
    tr_state = ""
    class _InfoValue(Any):
        pass
    results: dict[str, list[tuple[str, _InfoValue]]] = {}
    for tr in vcf.fetch():
        tr_filter = tr.filter.keys().pop(0)
        if tr_filter != "PASS":
            continue
        counter += 1
        assert tr.id is not None
        tr_state = str(tr.filter.keys().pop())
        # if tr.filter.get("PASS") is None:
        #     continue
        # if tr.info.get("CAND_CLASS") not in ["smear"]:
        #     continue
        [locus_id, n_allele] = tr.id.split(".") 
        region = ({tr.contig}, {tr.pos})
        if "" == current_lid:
            # 1st
            current_lid = locus_id
            locus = make_df(tr.samples)
        elif current_lid == locus_id:
            # same
            locus = locus.vstack(pl.DataFrame({tissue.name: tissue.get("AD")[1] for tissue in tr.samples.values()}))
            pass
        else:
            nloci += 1
            # process, dump ...
            results.update(process_locus(current_lid, locus, tr_state))
            # and re-start
            current_lid = locus_id
            locus = make_df(tr.samples)
    results.update(process_locus(current_lid, locus, tr_state))

    for tr in vcf.fetch():
        if tr.id in results:
            for key, val in results[tr.id]:
                if val is not None:
                    tr.info[key] = val
            add_gene = get_gene(tr.id, genes, catalog)
            if "" !=  add_gene:
                tr.info["GENE"] = get_gene(tr.id, genes, catalog)
            out.write(tr)

    out.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmd", required=True, type=str)     # command to use
    parser.add_argument("--infile", required=True, type=str)  # VCF input
    # Extra
    parser.add_argument("--donor", required=False, type=str, default="")   # donor name/id   <= needed in parse_table
    parser.add_argument("--tr-bed", required=False, type=str, default="")  # TR catalog      <= needed in parse_vcf
    parser.add_argument("--genes", required=False, type=str, default="")   # genes for annot <= needed in parse_vcf

    return TRConfig(parser.parse_args())


def main(config: TRConfig):
    log.info(f'Using command "{config.cmd}", with args: {config.infile} {config.donor}')
    if "vcf" == config.cmd:
        parse_vcf(config)
    elif "tab" == config.cmd:
        parse_table(config)


if "__main__" == __name__:
    config = parse_args()
    assert config.cmd in {'vcf', 'tab'}, f'Unknown command {config.cmd}, available: "vcf", "tab" '
    assert os.path.isfile(config.infile), f'File {config.infile} does not exixist'
    main(config)
