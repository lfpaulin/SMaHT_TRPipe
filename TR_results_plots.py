#!/usr/bin/env python
import os
import sys
import json
import argparse
import numpy as np
import polars as pl
from pathlib import Path
from pysam import VariantFile as VCF
from pysam import VariantRecord as Var


aqui = "TODO"


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
    "Colon_Asc": "Endoderm",
    "Colon_Desc": "Endoderm",
    "Liver_Sample": "Endoderm",
    "Adrena_GL": "Ectoderm",
    "Adrena_GR": "Ectoderm",
    "Aorta_Abd": "Mesoderm",
    "Lung": "Endoderm",
    "Heart_LV": "Mesoderm",
    "Testis_L": "Other",
    "Testis_R": "Other",
    "Ovary_L": "Other",
    "Ovary_R": "Other",
    "Fibroblast": "Other",
    "Skin_Calf": "Ectoderm",
    "Skin_Abdomen": "Ectoderm",
    "Muscle": "Mesoderm",
    "Brain": "Ectoderm",
    "Front_Lobe_BrLH": "Ectoderm",
    "Temp_Lobe_BrLH": "Ectoderm",
    "Cereb_BrLH": "Ectoderm",
    "Hippoc_BrLH": "Ectoderm",
    "Hippoc_BrRH": "Ectoderm"    
}


tissue_codes = {
    "3A": "BLOO",
    "3B": "BUCC",
    "3C": "ESOP",
    "3E": "COAS",
    "3G": "COAD",
    "3I": "LIVR",
    "3K": "ADGL",
    "3M": "ADGR",
    "3O": "AORT",
    "3Q": "LUNG",
    "3S": "HART",
    "3U": "TESL",
    "3W": "TESR",
    "3Y": "OVAL",
    "3AA": "OVAR",
    "3AC": "FBRO",
    "3AD": "SKSE",
    "3AF": "SKNE",
    "3AH": "MUSC",
    "3AK": "BRFL",
    "3AL": "BRTL",
    "3AM": "BRCE",
    "3AN": "BRHL",
    "3AO": "BRHR",
    "BLOO": "3A",
    "BUCC": "3B",
    "ESOP": "3C",
    "COAS": "3E",
    "COAD": "3G",
    "LIVR": "3I",
    "ADGL": "3K",
    "ADGR": "3M",
    "AORT": "3O",
    "LUNG": "3Q",
    "HART": "3S",
    "TESL": "3U",
    "TESR": "3W",
    "OVAL": "3Y",
    "OVAR": "3AA",
    "FBRO": "3AC",
    "SKSE": "3AD",
    "SKNE": "3AF",
    "MUSC": "3AH",
    "BRFL": "3AK",
    "BRTL": "3AL",
    "BRCE": "3AM",
    "BRHL": "3AN",
    "BRHR": "3AO"
}


PROTOCOLS = {
    "3A": {"tissue": "Blood, Whole", "preservation": "Snap Frozen", "layer": "Clinically_accessible",
           "tissue_abv": "BLOO", "DedupKey": "D01", "color": "#d22a2b"},
    "3B": {"tissue": "Buccal Swab", "preservation": "Fresh", "layer": "Clinically_accessible",
           "tissue_abv": "BUCC", "DedupKey": "D02", "color": "#e97451"},
    "3C": {"tissue": "Esophagus", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "ESOP", "DedupKey": "D03", "color": "#d9e58d"},
    "3D": {"tissue": "Esophagus", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "ESOP", "DedupKey": "D03", "color": "#d9e58d"},
    "3E": {"tissue": "Colon, Ascending", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "COAS", "DedupKey": "D04", "color": "#a7d16b"},
    "3F": {"tissue": "Colon, Ascending", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "COAS", "DedupKey": "D04", "color": "#a7d16b"},
    "3G": {"tissue": "Colon, Descending", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "COAD", "DedupKey": "D04", "color": "#2ab463"},
    "3H": {"tissue": "Colon, Descending", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "COAD", "DedupKey": "D04", "color": "#2ab463"},
    "3I": {"tissue": "Liver", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "LIVR", "DedupKey": "D05", "color": "#0a984f"},
    "3J": {"tissue": "Liver", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "LIVR", "DedupKey": "D05", "color": "#0a984f"},
    "3K": {"tissue": "Adrenal Gland, Left", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "ADGL", "DedupKey": "D06", "color": "#c2a5cf"},
    "3L": {"tissue": "Adrenal Gland, Left", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "ADGL", "DedupKey": "D06", "color": "#c2a5cf"},
    "3M": {"tissue": "Adrenal Gland, Right", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "ADGR", "DedupKey": "D06", "color": "#9970ab"},
    "3N": {"tissue": "Adrenal Gland, Right", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "ADGR", "DedupKey": "D06", "color": "#9970ab"},
    "3O": {"tissue": "Aorta, Abdominal", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "AORT", "DedupKey": "D07", "color": "#762a83"},
    "3P": {"tissue": "Aorta, Abdominal", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "AORT", "DedupKey": "D07", "color": "#762a83"},
    "3Q": {"tissue": "Lung", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "LUNG", "DedupKey": "D08", "color": "#006937"},
    "3R": {"tissue": "Lung", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "LUNG", "DedupKey": "D08", "color": "#006937"},
    "3S": {"tissue": "Heart, LV", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "HART", "DedupKey": "D09", "color": "#3f124b"},
    "3T": {"tissue": "Heart, LV", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "HART", "DedupKey": "D09", "color": "#3f124b"},
    "3U": {"tissue": "Testis, Left", "preservation": "Snap Frozen", "layer": "Germinal",
           "tissue_abv": "TESL", "DedupKey": "D10", "color": "#b2d8d8"},
    "3V": {"tissue": "Testis, Left", "preservation": "Fixed", "layer": "Germinal",
           "tissue_abv": "TESL", "DedupKey": "D10", "color": "#b2d8d8"},
    "3W": {"tissue": "Testis, Right", "preservation": "Snap Frozen", "layer": "Germinal",
           "tissue_abv": "TESR", "DedupKey": "D10", "color": "#66b2b2"},
    "3X": {"tissue": "Testis, Right", "preservation": "Fixed", "layer": "Germinal",
           "tissue_abv": "TESR", "DedupKey": "D10", "color": "#66b2b2"},
    "3Y": {"tissue": "Ovary, Left", "preservation": "Snap Frozen", "layer": "Germinal",
           "tissue_abv": "OVAL", "DedupKey": "D11", "color": "#ffc2cd"},
    "3Z": {"tissue": "Ovary, Left", "preservation": "Fixed", "layer": "Germinal",
           "tissue_abv": "OVAL", "DedupKey": "D11", "color": "#ffc2cd"},
    "3AA": {"tissue": "Ovary, Right", "preservation": "Snap Frozen", "layer": "Germinal",
            "tissue_abv": "OVAR", "DedupKey": "D11", "color": "#ff93ac"},
    "3AB": {"tissue": "Ovary, Right", "preservation": "Fixed", "layer": "Germinal",
            "tissue_abv": "OVAR", "DedupKey": "D11", "color": "#ff93ac"},
    "3AC": {"tissue": "Dermal Fibroblast", "preservation": "Cultured Cells",  "layer": "Cultured_Cells",
            "tissue_abv": "FBRO", "DedupKey": "D12", "Notes": "Isolated from fresh calf skin", "color": "#ffeead"},
    "3AD": {"tissue": "Skin, Calf", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "SKSE", "DedupKey": "D13", "color": "#bf812b"},
    "3AE": {"tissue": "Skin, Calf", "preservation": "Fixed", "layer": "Ectoderm",
            "tissue_abv": "SKSE", "DedupKey": "D13", "color": "#bf812b"},
    "3AF": {"tissue": "Skin, Abdomen", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "SKNE", "DedupKey": "D13", "color": "#dfc17c"},
    "3AG": {"tissue": "Skin, Abdomen", "preservation": "Fixed", "layer": "Ectoderm",
            "tissue_abv": "SKNE", "DedupKey": "D13", "color": "#dfc17c"},
    "3AH": {"tissue": "Muscle", "preservation": "Snap Frozen", "layer": "Mesoderm",
            "tissue_abv": "MUSC", "DedupKey": "D14", "color": "#ffcc5c"},
    "3AI": {"tissue": "Muscle", "preservation": "Fixed", "layer": "Mesoderm",
            "tissue_abv": "MUSC", "DedupKey": "D14", "color": "#ffcc5c"},
    "3AJ": {"tissue": "Brain", "preservation": "Fresh", "layer": "Ectoderm",
            "tissue_abv": "BRAI", "DedupKey": "D15", "color": "#000000"},
    "3AK": {"tissue": "Frontal Lobe, Brain, Left hemisphere", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRFL", "DedupKey": "D15", "color": "#e9f2f4"},
    "3AL": {"tissue": "Temporal Lobe, Brain, Left hemisphere", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRTL", "DedupKey": "D15", "color": "#addae9"},
    "3AM": {"tissue": "Cerebellum, Brain, Left hemisphere", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRCE", "DedupKey": "D15", "color": "#75add1"},
    "3AN": {"tissue": "Hippocampus, Brain, Left hemisphere", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRHL", "DedupKey": "D15", "color": "#4676b5"},
    "3AO": {"tissue": "Hippocampus, Brain, Right hemisphere", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRHR", "DedupKey": "D15", "color": "#323a94"},
}


DEDUP_LOOKUP = {
    'D01': {'tissue': 'Blood', 'tissue_abv': 'BLO'},
    'D02': {'tissue': 'Buccal Swab', 'tissue_abv': 'BUC'},
    'D03': {'tissue': 'Esophagus', 'tissue_abv': 'ESO'},
    'D04': {'tissue': 'Colon', 'tissue_abv': 'COA'},
    'D05': {'tissue': 'Liver', 'tissue_abv': 'LIV'},
    'D06': {'tissue': 'Adrenal Gland', 'tissue_abv': 'ADG'},
    'D07': {'tissue': 'Aorta', 'tissue_abv': 'AOR'},
    'D08': {'tissue': 'Lung', 'tissue_abv': 'LUN'},
    'D09': {'tissue': 'Heart', 'tissue_abv': 'HAR'},
    'D10': {'tissue': 'Testis', 'tissue_abv': 'TES'},
    'D11': {'tissue': 'Ovary', 'tissue_abv': 'OVA'},
    'D12': {'tissue': 'Fibroblast', 'tissue_abv': 'FBR'},
    'D13': {'tissue': 'Skin', 'tissue_abv': 'SKI'},
    'D14': {'tissue': 'Muscle', 'tissue_abv': 'MUS'},
    'D15': {'tissue': 'Brain', 'tissue_abv': 'BRA'},
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
    "B001": "PB",
    "D001": "ONT",
    "D003": "ONT",
    "B004": "PB"
}


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

log = MyLog()


class TRPlotConfig:
    def __init__(self, args: argparse.Namespace ):
        self.donor = args.donor
        self.vcf = args.donor_vcf
        self.pq = args.candidate_pq
        self.use_loci = args.selected_loci
        self.plot_only = args.plot_only
        self.merge_tissue = args.merge_tissue
        self.trplot_r = "TRPipe/TR_plotting.R"


from datetime import datetime
class DevTiming(object):
    def __init__(self, timing_process: str = ""):
        self.start_ts: float
        self.timing_process = timing_process if timing_process != "" else "process"

    def start(self):
        # Output
        log.info(f'Starting {self.timing_process} ...')
        # time start for process
        start_dt = datetime.now()
        self.start_ts = float(datetime.timestamp(start_dt))

    def end(self, add_sep=True, return_time=False):
        end_dt = datetime.now()
        end_ts = float(datetime.timestamp(end_dt))
        separator = "# ########################## #"
        time_spent = np.round(end_ts - self.start_ts, 3)
        log.info(f'Timing {time_spent}\n{separator}') if add_sep \
            else log.info(f'Timing {time_spent}')
        if return_time:
            return time_spent


class PlotOut:
    def __init__(self, donor) -> None:
        pass
        self.donor = donor
        self.rmeta = f'{donor}_trplot_meta.json'
        self.plot = ""
        self.tab = ""
        self.meta = ""

    def update(self, lid):
        self.plot = f'{self.donor}_trplot_{lid:08d}.png'
        self.tab = f'{self.donor}_trplot_{lid:08d}.tsv'
        self.meta = f'{self.donor}_trplot_{lid:08d}.txt'


def name_translation(sample_name: str, a2b: str) -> str:
    """
    sample_name: str
    expected format: DONOR-TISSUE-CORE-SEX/AGE-TECH-GCC
    example SMHT001-3A-XX-M42-B001-broad
    """
    if "tissue" == a2b:
        return tissue_name_converter[sample_name.split("-")[1]]
    elif "tech" == a2b:
        return tissue_tech_converter[sample_name.split("-")[4]]
    elif "abv" == a2b:
        return tissue_codes[sample_name.split("-")[1]]
    else:
        return ""


def from_vcf_2_plot_lids(trplot_config: TRPlotConfig) -> tuple[set[int], dict[int,list[tuple[str,str]]]]:
    use_dev = True
    selec_tr_annot = {"CORE", "TISSUE", "DERM"}
    vcf = VCF(config.vcf)
    var: Var
    select_tr_byID = set()
    select_tr_byID_info = dict()
    try:
        with open(config.use_loci, "r") as file:
            use_tr_loci = [int(line.rstrip('\n')) for line in file]
    except FileNotFoundError:
        log.info(f'Using all loci, no "selected loci" file provided {config.use_loci}')
        use_tr_loci = []
    for var in vcf.fetch():
        tr_annot = var.info.get("TR_ANNOT")
        if tr_annot in selec_tr_annot or use_dev:
            lid = var.info.get("LOCUSID")
            if len(use_tr_loci) > 0:
                if lid not in use_tr_loci:
                    continue
            select_tr_byID.add(lid)
            tr_len = var.info.get("SVLEN")
            tr_annot_info = var.info.get("TR_ANNOT_INFO")
            tr_annot_locus = var.info.get("TR_ANNOT_LOCUS")
            gene = var.info.get("GENE")
            if lid not in select_tr_byID_info:
                select_tr_byID_info[lid] = [(tr_annot, tr_annot_info, tr_annot_locus, gene, var.id, tr_len)]
            else:
                select_tr_byID_info[lid].append((tr_annot, tr_annot_info, tr_annot_locus, gene, var.id, tr_len))
    return select_tr_byID, select_tr_byID_info, 


def derm_order(tissues: list[str]) -> list[str]:
    use_order = [
        "BLOO", "AORT", "HART", "MUSC", "ESOP", "LUNG", 
        "LIVR", "COAS", "COAD", "ADGL", "ADGR", "SKSE", 
        "SKNE", "BRFL", "BRTL", "BRCE", "BRHL", "BRHR", 
        "FBRO", "BUCC", "TESL", "TESR", "OVAL", "OVAR"
    ]
    torder_idx = set()
    torder_tiss = {}
    for tiss in tissues:
        tabv_oder_index = use_order.index(name_translation(tiss, "abv"))
        torder_idx.add(tabv_oder_index)
        if tabv_oder_index not in torder_tiss:
            torder_tiss[tabv_oder_index] = [tiss]
        else:
            torder_tiss[tabv_oder_index].append(tiss)

    return [ tiss for ti in torder_idx for tiss in torder_tiss[ti] ]


def from_dict_to_table(lid: int, tissue_counts: dict, largets_counts: int, tissue_list: list, plot: PlotOut, merge_tissue: bool = False):
    # needs to fill empty with NAs
    # use largets_counts is the height of the table
    df_counts = {}
    df_haps = {}
    sname_list = []
    largets_counts_merge = 0
    for k in tissue_list:
        sname, counts, haps = tissue_counts[k]
        if merge_tissue:
            if sname not in df_counts:
                df_counts[sname] = counts
                sname_list.append(sname)
                df_haps[sname] = haps
                largets_counts_merge = len(counts) if largets_counts_merge < len(counts) else largets_counts_merge
            else:
                df_counts[sname] += counts
                df_haps[sname] += haps
                largets_counts_merge = len(df_counts[sname]) if largets_counts_merge < len(df_counts[sname]) else largets_counts_merge
        else:
            use_sname = f'{sname}_{name_translation(k, "tech")}'
            if len(counts) < largets_counts:
                counts += [-1] * (largets_counts - len(counts))
                haps += [-1] * (largets_counts - len(haps))
            if use_sname not in df_counts:
                df_counts[use_sname] = counts
                sname_list.append(use_sname)
                df_haps[use_sname] = haps
            else:
                new_sname = f'{use_sname}2'
                df_counts[new_sname] = counts
                sname_list.append(new_sname)
                df_haps[new_sname] = haps
    if merge_tissue:
        for k in tissue_list:
            sname, counts, haps = tissue_counts[k]
            if len(df_counts[sname]) < largets_counts_merge:
                df_counts[sname] += [-1] * (largets_counts_merge - len(df_counts[sname]))
                df_haps[sname] += [-1] * (largets_counts_merge - len(df_haps[sname]))

    df = pl.DataFrame({k: pl.Series(df_counts[k]) for k in sname_list})
    df.write_csv(plot.tab, separator="\t")
    df = pl.DataFrame({k: pl.Series(df_haps[k]) for k in sname_list})
    df.write_csv(f'{plot.tab}2', separator="\t")


def make_meta(lid, metadata, plot: PlotOut):
    gidx = 3
    gene = metadata[0][gidx] if metadata[0][gidx] is not None else ""
    trlidx = 5
    trlen = ", ".join([ str(m[trlidx]) if m[trlidx] is not None else "0" for m in metadata ])
    out = open(plot.meta, "w")
    out.write(f'lid <- "{lid}"\n')
    out.write(f'gene <- "{gene}"\n')
    out.write(f'trlen <- c({trlen})\n')
    out.close()


def srun(cmd):
    import subprocess
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        log.err(result.stderr)
    # log.info(result.stdout)


def main(config: TRPlotConfig):
    plot = PlotOut(config.donor)
    final_lids, final_lids_info = from_vcf_2_plot_lids(config)
    with open(plot.rmeta, 'w') as f:
        json.dump(final_lids_info, f)
    tab_canidate = pl.read_parquet(config.pq)
    nfinal_lid = len(final_lids)
    lid_count = 1
    for lid in sorted(list(final_lids)):
        larger_count = 0
        tissue_counts = {}
        tissue_list = []
        plot.update(lid)
        log.info(f'Starting analysis for LocusID {lid} ({lid_count}/{nfinal_lid})')
        lid_count += 1
        if config.plot_only and os.path.exists(plot.tab):
            pass
        else:
            for sname in tab_canidate["SampleName"].unique().to_list():
                prettyName = name_translation(sname, "abv")
                counts = tab_canidate.filter((pl.col("LocusID") == lid) & (pl.col("SampleName") == sname))["ReadLen"].to_list()
                haps = tab_canidate.filter((pl.col("LocusID") == lid) & (pl.col("SampleName") == sname))["HP"].to_list()
                larger_count = len(counts) if len(counts) > larger_count else larger_count
                tissue_counts[sname] = (prettyName, counts, haps)
                tissue_list.append(sname)
            tissue_list = derm_order(tissue_list)
            from_dict_to_table(lid, tissue_counts, larger_count, tissue_list, plot, config.merge_tissue)
            make_meta(lid, final_lids_info[lid], plot)
        flag_merged = "--merget" if config.merge_tissue else ""
        srun(f'{config.trplot_r}  --tab {plot.tab}  --out {plot.plot} --meta {plot.meta} --donor {config.donor} {flag_merged}')
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--donor", required=True, type=str)            # donor name/id
    parser.add_argument("--donor-vcf", required=True, type=str)        # VCF output from pipeline
    parser.add_argument("--candidate-pq", required=True, type=str)     # pq file with candidates (step2)
    parser.add_argument("--selected-loci", required=False, type=str, default="")  # use only these loci
    parser.add_argument("--merge-tissue", required=False, default=False, action="store_true")  # merge samples from same tissue
    parser.add_argument("--plot-only", required=False, default=False, action="store_true")
    return TRPlotConfig(parser.parse_args())


if "__main__" == __name__:
    config = parse_args()
    assert os.path.isfile(config.pq), f'File {config.pq} does not exixist'
    assert os.path.isfile(config.vcf), f'File {config.vcf} does not exixist'
    main(config)
