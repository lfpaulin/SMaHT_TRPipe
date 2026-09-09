#!/usr/bin/env python
import os
import sys
import json
import gzip
import pickle
import numpy as np
import argparse
from sklearn.cluster import MeanShift

os.environ["POLARS_MAX_THREADS"] = "8"
import polars as pl


tissue_name_converter = {
    "3A": "Blood,Whole",
    "3B": "Buccal_Swab",
    "3C": "Esophagus",
    "3D": "Esophagus",
    "3E": "Colon,Ascending",
    "3F": "Colon,Ascending",
    "3G": "Colon,Descending",
    "3H": "Colon,Descending",
    "3I": "Liver_Sample",
    "3J": "Liver_Sample",
    "3K": "Adrenal_Gland,Left",
    "3L": "Adrenal_Gland,Left",
    "3M": "Adrenal_Gland,Right",
    "3N": "Adrenal_Gland,Right",
    "3O": "Aorta,Abdominal",
    "3P": "Aorta,Abdominal",
    "3Q": "Lung",
    "3R": "Lung",
    "3S": "Heart,LV",
    "3T": "Heart,LV",
    "3U": "Testis,Left",
    "3V": "Testis,Left",
    "3W": "Testis,Right",
    "3X": "Testis,Right",
    "3Y": "Ovary,Left",
    "3Z": "Ovary,Left",
    "3AA": "Ovary,Right",
    "3AB": "Ovary,Right",
    "3AC": "Dermal_Fibroblast",
    "3AD": "Skin,Calf",
    "3AE": "Skin,Calf",
    "3AF": "Skin,Abdomen",
    "3AG": "Skin,Abdomen",
    "3AH": "Muscle",
    "3AI": "Muscle",
    "3AJ": "Brain",
    "3AK": "Frontal_Lobe,Brain,Left_hemisphere",
    "3AL": "Temporal_Lobe,Brain,Left_hemisphere",
    "3AM": "Cerebellum,Brain,Left_hemisphere",
    "3AN": "Hippocampus,Brain,Left_hemisphere",
    "3AO": "Hippocampus,Brain,Right_hemisphere",
    "3AP": "Frontal_Lobe,Brain,Left_hemisphere",
    "3AQ": "Temporal_Lobe,Brain,Left_hemisphere",
    "3AR": "Cerebellum,Brain,Left_hemisphere",
    "3AS": "Hippocampus,Brain,Left_hemisphere",
    "3AT": "Hippocampus,Brain,Right_hemisphere"
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
        sys.exit(1)

log = MyLog()


class TRConfig:
    def __init__(self):
        self.donors: str = ""
        self.donors_meta: str = ""
        self.data: str = ""
        self.res: str = ""
        self.date: str = ""
        self.loci: str = ""
        self.tr_bed: str = ""
        self.noqc: bool = False
        self.min_reads_sample = 5

    def set(self, args: argparse.Namespace):
        self.donors = args.donors_list
        self.donors_meta = args.donors_meta
        self.data = args.data_dir
        self.res = args.res_dir
        self.date = args.date
        self.loci = args.selected_loci
        self.tr_bed = args.tr_bed
        self.noqc = args.noqc

    def setu(self, arg_tupple: tuple[str,str,str,str,str,str,str]):
        #donors_list, donors_meta, data_dir, res_dir, date, selected_loci, tr_bed
        self.donors, self.donors_meta, self.data, self.res, self.date, self.loci, self.tr_bed = arg_tupple


class TRPaths:
    def __init__(self, config: TRConfig):
        self.tr_read_lens = config.data
        # these are TSV tables
        self.read_merge_path: list = []
        self.read_trgt_path: list = []
        self.read_medaka_path: list = []
        # schema for the TSV loading in polars, for meme efficiency
        self.complete_schema = {
            "Donor": pl.String,
            "Tissue": pl.String,
            "Core": pl.String,
            "Contig": pl.String,
            "Pos": pl.UInt32,
            "End": pl.UInt32,
            "LocusID": pl.UInt32,
            "HP": pl.UInt8,
            "ReadLen": pl.Int32
        }
        self.adotto_v2_flank_filter = 602361
        self.adotto_v2_ff_tr_ids = 1795951
        self.adotto_v2_conitgs = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"]
        self.adotto_v2_conitgs_wlen = [("chr1", 248956422), ("chr2", 242193529), ("chr3", 198295559), ("chr4", 190214555), ("chr5", 181538259), ("chr6", 170805979), ("chr7", 159345973), ("chr8", 145138636), ("chr9", 138394717), ("chr10", 133797422), ("chr11", 135086622), ("chr12", 133275309), ("chr13", 114364328), ("chr14", 107043718), ("chr15", 101991189), ("chr16", 90338345), ("chr17", 83257441), ("chr18", 80373285), ("chr19", 58617616), ("chr20", 64444167), ("chr21", 46709983), ("chr22", 50818468), ("chrX", 156040895), ("chrY", 57227415)]
        self.adotto_v2_path = config.tr_bed

    def load_read_files(self):
        self.read_merge_path = [f for f in os.listdir(self.tr_read_lens) if os.path.isfile(os.path.join(self.tr_read_lens, f))]
        dont_exist = [f for f in os.listdir(self.tr_read_lens) if not os.path.isfile(os.path.join(self.tr_read_lens, f))]
        if len(dont_exist) > 0:
            log.info(f'These files do not exist or these are not files: {dont_exist}')
        assert len(self.read_merge_path) - len(dont_exist) > 0, f'No files found in {self.tr_read_lens}: {self.read_merge_path}'


padding = {
    "trgt": 50,
    "medaka": 250
}


def load_adotto_bed(bed_path: str):
    tr_catalog = {}
    bed_file = gzip.open(bed_path, "rt") if "gz" in bed_path else open(bed_path, "r")
    for line in bed_file:
        contig, pos, end, info = line.rstrip("\n").split("\t")
        trid, motif, *_ = info.split(";")
        _, trid = trid.split("=")
        _, motif = motif.split("=")
        tr_catalog[trid] = (contig, int(pos), int(end), motif)
    bed_file.close()
    return tr_catalog


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


def clusteirng(allele_lens: np.ndarray, bandwidth: int = 20, return_uval_raw: bool = False):
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(allele_lens)
    # Extract cluster labels and center coordinates
    uval_, ucount_ = np.unique(ms.labels_, return_counts=True)
    uvalv_ = [ms.cluster_centers_[v][0] for v in uval_]
    # ulabel_schema = {"Mean": pl.Float64, "Counts": pl.UInt16}
    # ulabel = pl.DataFrame({"Mean": uval, "Counts": ucount}, schema=ulabel_schema)
    # filters
    has_negative_len_ = sum([v < 0 for v in uvalv_]) > 0
    with_single_read_ = sum([c <= 2 for c in ucount_])
    if not return_uval_raw:
        uval_ = uvalv_

    return ms.cluster_centers_, ms.labels_, uval_, ucount_, has_negative_len_, with_single_read_


def find_candidates(tr_reads: pl.DataFrame, bandwidth = 20, check_many_clusters=7):
    candidate = []
    candidate_stats = {}
    loci_list = tr_reads["LocusID"].unique().to_list()
    for locus in loci_list:
        lid_alleles = tr_reads.filter(pl.col("LocusID") == locus)["ReadLen"].to_numpy()
        lid_alleles = lid_alleles.reshape(-1, 1)
        cluster_centers, labels, uval, ucount, has_negative_len, with_single_read = clusteirng(lid_alleles, bandwidth)
        num_clusters = len(np.unique(labels))
        # strict germline
        if num_clusters > 1 and not has_negative_len:
            # filter only two alleles
            if num_clusters == 2:
                [x, y] = ucount
                if x/(x+y) - 0.5 < 0.05:
                    continue
            # ################################################
            # re-cluster those with high cluster number
            if num_clusters >= check_many_clusters:
                [x, y] = [uval[i] for i in range(2)]
                bw = max(int(abs(x - y)/10), bandwidth)
                log.info(f"Re-clustering candidates")
                log.info(f'PRE: {locus} {bandwidth} v {bw}, {uval}, {ucount}, {num_clusters}')
                if bw > bandwidth:
                    bandwidth = bw
                    cluster_centers, lables, uval, ucount, has_negative_len, with_single_read = clusteirng(lid_alleles, bandwidth)
                    num_clusters = len(np.unique(labels))
                    log.info(f'AFT: {locus} {bandwidth}, {uval}, {ucount}, {num_clusters}')
                    # strict 2 clusters again
                    if num_clusters > 1 and not has_negative_len:
                        # filter only two alleles
                        if num_clusters == 2:
                            [x, y] = ucount
                            if x/(x+y) - 0.5 < 0.05:
                                continue
            # ################################################

            candidate.append(locus)
            candidate_stats[locus] = (bandwidth, cluster_centers, uval, ucount, with_single_read)
            if (num_clusters - with_single_read) >= 5:
                log.info(f'Candidate: {locus} with {num_clusters}.')
    return candidate, candidate_stats


class TRFilters:
    def __init__(self):
        # single len per hap, both
        self.double_single_allele = []
        # single len in one hap, 2 in the other
        self.single_allele = []
        # single len in one allele, +3 in the other
        self.single_allele_smear = []
        # 2 alleles 85%
        self.germline_allele = []
        # 2 alleles 85%, +1 lower VAF
        self.mosaic_allele = []
        # 3 alleles 85%
        self.third_allele = []
        self.third_allele_mosaic = []
        # +4 alleles 85%
        self.smear = []
        # none of the above
        self.unknown = []
        # read threshold < 10 reads for all mosaic/smear alleles
        self.low_support_global = []

    def __repr__(self):
        # single len per hap, both
        double_single_allele = len(self.double_single_allele)
        # single len in one hap, 2 in the other
        single_allele = len(self.single_allele)
        # single len in one allele, +3 in the other
        single_allele_smear = len(self.single_allele_smear)
        # 2 alleles 85%
        germline_allele = len(self.germline_allele)
        # 2 alleles 85%, +1 lower VAF
        mosaic_allele = len(self.mosaic_allele)
        # 3 alleles 85%
        third_allele = len(self.third_allele)
        third_allele_mosaic = len(self.third_allele_mosaic)
        # +4 alleles 85%
        smear = len(self.smear)
        # none of the above
        unknown = len(self.unknown)
        # read threshold < 10 reads for all mosaic/smear alleles
        low_support_global = len(self.low_support_global)
        return (
            f"Double-Single length-allele per haplotype: {double_single_allele}\n"
            f"Single length-allele in one haplotype, 2 in the other: {single_allele}\n"
            f"Single length-allele in one haplotype, +3 in the other: {single_allele_smear}\n"
            f"Two alleles represent 85% VAF: {germline_allele}\n"
            f"Two alleles represent 85% VAF, +1 w/lower VAF: {mosaic_allele}\n"
            f"Three alleles represent 85%: {third_allele}\n"
            f"Three alleles represent 85% + mosaic: {third_allele_mosaic}\n"
            f"+4 alleles represent 85%: {smear}\n"
            f"None of the above: {unknown}\n"
            f"\n"
            f"Low support global: {low_support_global}"
        )

    def update(self, tr_candidate_list):
        self.double_single_allele.extend(tr_candidate_list.double_single_allele)
        self.single_allele.extend(tr_candidate_list.single_allele)
        self.single_allele_smear.extend(tr_candidate_list.single_allele_smear)
        self.germline_allele.extend(tr_candidate_list.germline_allele)
        self.mosaic_allele.extend(tr_candidate_list.mosaic_allele)
        self.third_allele.extend(tr_candidate_list.third_allele)
        self.third_allele_mosaic.extend(tr_candidate_list.third_allele_mosaic)
        self.smear.extend(tr_candidate_list.smear)
        self.unknown.extend(tr_candidate_list.unknown)
        self.low_support_global.extend(tr_candidate_list.low_support_global)

    def filters(self, tr_class):
        filter_pass = {"single_allele_smear", "mosaic_allele", "third_allele_mosaic", "smear"}
        if hasattr(self, tr_class):
            if tr_class in filter_pass:
                return "PASS"
            elif tr_class == "unknown":
                return "UNKNOWN"
            else:
                return "GERMLINE"
        else:
            log.err(f'"{tr_class}" not found in TRFilters')
            return "UNKNOWN"

    def get(self, tr_class):
        if hasattr(self, tr_class):
            return getattr(self, tr_class, [])
        else:
            log.err(f'"{tr_class}" not found in TRFilters')
            return []

    def keys(self):
        return vars(self).keys()

    def jsonfy(self, outputf: str):
        with open(outputf, "w") as fout:
            json.dump(vars(self), fout, indent=4)
        fout.close()

    def from_json(self, inputf: str):
        with open(inputf, 'r') as fin:
            # Parse the JSON file directly
            data = json.load(fin)
        fin.close()
        for k in data.keys():
            setattr(self, k, data[k])


def filter_candidates(tr_candidate_reads: pl.DataFrame, min_reads_overall: int, save_germline: bool = False, cand_stats: dict = {}, min_reads_global_somatic = 10, min_reads_global_mosaic = 3):
    applied_filter = TRFilters()
    use_raw_uval = True
    loci_list = tr_candidate_reads["LocusID"].unique().to_list()
    def_bandwidth = 20
    for locus in loci_list:
        bandwidth = def_bandwidth if len(cand_stats) == 0 else cand_stats[locus][0]
        tr_locus = tr_candidate_reads.filter(pl.col("LocusID") == locus)
        lid_alleles = tr_locus["ReadLen"].to_numpy()
        if len(lid_alleles.tolist()) < min_reads_overall:
            applied_filter.unknown.append(locus)
            continue
        lid_alleles = lid_alleles.reshape(-1, 1)
        cluster_centers, labels, uval, ucount, has_negative_len, with_single_read = clusteirng(lid_alleles, bandwidth, use_raw_uval)
        num_clusters = len(np.unique(labels))
        if max(uval) >= len(ucount):
            uval, ucount = fill_spaces(uval, ucount)
        uval = np.array([cluster_centers[v][0] for v in uval])

        tr_locus = tr_locus.with_columns(pl.Series("Label", labels), pl.Series("Length-ish", ucount[labels]))

        hap_multi_allele_count = 0
        haps = []
        hp_single = []
        hapm = []
        hp_multi = []
        alleles_by_hap = []
        seen_haps = tr_locus["HP"].unique().to_list()
        for hap in seen_haps:
            alleles_by_hap = tr_locus.filter(pl.col("HP") == hap)["Label"].unique().to_list()
            if len(alleles_by_hap) == 1:
                if hap == 0:
                    hap_multi_allele_count += 1
                else:
                    hp_single.append(uval[alleles_by_hap])
                    haps.append(hap)
            else:
                hp_multi.append(uval[alleles_by_hap])
                hapm.append(hap)
                hap_multi_allele_count += 1
        n_seen_alleles = len(uval)
        add_meta = ""
        trclasstr = ""
        allele_index = 0
        if hap_multi_allele_count == 0:
            # 1 allele size per hap, and no unphased
            if save_germline:
                applied_filter.double_single_allele.append(locus)
                log.info(f"LocusID = {locus}: Double Single-length haplotype: {haps},{hp_single} | {n_seen_alleles}")
            else:
                applied_filter.double_single_allele.append(1)
            continue
        elif hap_multi_allele_count == 1 and n_seen_alleles == 1 and 0 in seen_haps:
            # 1 allele size per hap, and some unphased
            if save_germline:
                applied_filter.double_single_allele.append(locus)
                log.info(f"LocusID = {locus}: Double Single-length haplotype: {haps},{hp_single} | {n_seen_alleles}")
            else:
                applied_filter.double_single_allele.append(1)
            continue
        elif hap_multi_allele_count == 1:
            if 1 in hapm or 2 in hapm:
                # 1 allele size in hapA, 2 alleles in hapB
                applied_filter.single_allele.append(locus)
                trclasstr = f'Single-length haplotype'
                add_meta = ", Mosaic" if n_seen_alleles <= 3 else ", Smear"
                allele_index = 2  # 0 and 1 SHOULD always the germline alleles
                applied_filter.mosaic_allele.append(locus)
                if sum(ucount[allele_index:]) < min_reads_global_somatic or sum(ucount[allele_index:] < min_reads_global_mosaic) == len(ucount[allele_index:]):
                    applied_filter.low_support_global.append(locus)
                    add_meta += ", LOW SUPP GLOBAL"
            elif 0 in hapm:
                alleles_as_sinlge_hap = len(set([a[0] for a in hp_single]))
                if alleles_as_sinlge_hap == 0:
                    trclasstr = f'No phase info'
                    add_meta = ", in hap0, nope"
                    applied_filter.unknown.append(locus)
                elif alleles_as_sinlge_hap == 1:
                    # 1 allele size in hapA, 2 alleles in hapB
                    trclasstr = f'Double Single-length haplotype'
                    add_meta = ", in hap0, nope"
                    if save_germline:
                        applied_filter.double_single_allele.append(locus)
                    else:
                        applied_filter.double_single_allele.append(1)
                elif alleles_as_sinlge_hap == 2:
                    # 1 allele size in hapA, 1 alleles in hapB
                    trclasstr = f'Single-length haplotype'
                    if save_germline:
                        applied_filter.single_allele.append(locus)
                    else:
                        applied_filter.single_allele.append(1)
                else:
                    assert alleles_as_sinlge_hap > 2, f"LocusID = {locus}: Should not happen ==>> {haps},{hp_single} | {hapm},{hp_multi}{add_meta} | {n_seen_alleles}"
            else:
                if len(hapm) == 0 and [1, 2] == haps:
                    if save_germline:
                        applied_filter.single_allele.append(locus)
                    else:
                        applied_filter.single_allele.append(1)
                else:
                    log.warn(f'"else.slh", {locus} {n_seen_alleles}, {alleles_by_hap}, {seen_haps} | {haps},{hp_single} | {hapm},{hp_multi}{add_meta}')
            hp_multi = [float(v) for v in hp_multi[0]] if len(hp_multi) > 0 else hp_multi
            log.info(f"LocusID = {locus}: {trclasstr}: {haps},{hp_single} | {hapm},{hp_multi}{add_meta} | {n_seen_alleles}")
            continue
        else:
            pass
        germline_allele = 0.0
        germline_proportion = 0.85
        n_alleles= 0
        allele_counts = sum(ucount)
        while germline_allele < germline_proportion:
            germline_allele += ucount[allele_index] / allele_counts
            allele_index += 1
            n_alleles += 1
        if n_alleles <= 2:
            if allele_index == 1:
                log.info(f"LocusID = {locus}: Germline allele(s): {uval[:allele_index]}. Mosaic allele(s): {[float(v) for v in uval[allele_index:]]} VAF = {1-germline_allele} | {n_seen_alleles}")
                applied_filter.mosaic_allele.append(locus)
                if sum(ucount[allele_index:]) < min_reads_global_somatic or sum(ucount[allele_index:] < min_reads_global_mosaic) == len(ucount[allele_index:]):
                    applied_filter.low_support_global.append(locus)
                    add_meta += ", LOW SUPP GLOBAL"
            elif allele_index == 2 and germline_allele > 0.99:
                if save_germline:
                    applied_filter.germline_allele.append(locus)
                    log.info(f"LocusID = {locus}: Germline allele(s): {uval[:allele_index]}.")
                else:
                    applied_filter.germline_allele.append(1)
            elif allele_index == 2 < len(ucount):
                log.info(f"LocusID = {locus}: Germline allele(s): {uval[:allele_index]}. Mosaic allele(s): {[float(v) for v in uval[allele_index:]]} VAF = {1-germline_allele} | {n_seen_alleles}")
                applied_filter.mosaic_allele.append(locus)
                if sum(ucount[allele_index:]) < min_reads_global_somatic or sum(ucount[allele_index:] < min_reads_global_mosaic) == len(ucount[allele_index:]):
                    applied_filter.low_support_global.append(locus)
                    add_meta += ", LOW SUPP GLOBAL"
            else:
                pass
            continue
        else:
            if n_alleles == 3:
                if n_seen_alleles == 3:
                    log.info(f"LocusID = {locus}: Third allele: {uval[:allele_index]} | {n_seen_alleles}")
                    applied_filter.third_allele.append(locus)
                    if sum(ucount[allele_index:]) < min_reads_global_somatic or sum(ucount[allele_index:] < min_reads_global_mosaic) == len(ucount[allele_index:]):
                        applied_filter.low_support_global.append(locus)
                        add_meta += ", LOW SUPP GLOBAL"
                else:
                    log.info(f"LocusID = {locus}: Third allele: {uval[:allele_index]}, 4+ VAF = {1-germline_allele} {[float(v) for v in uval[allele_index:]]}| {n_seen_alleles}")
                    applied_filter.third_allele_mosaic.append(locus)
                    if sum(ucount[allele_index:]) < min_reads_global_somatic or sum(ucount[allele_index:] < min_reads_global_mosaic) == len(ucount[allele_index:]):
                        applied_filter.low_support_global.append(locus)
                        add_meta += ", LOW SUPP GLOBAL"
                continue
            else:
                log.info(f"LocusID = {locus}: Smear ?: {[float(v) for v in uval]} | {n_seen_alleles}")
                applied_filter.smear.append(locus)
                if sum(ucount[allele_index:]) < min_reads_global_somatic or sum(ucount[allele_index:] < min_reads_global_mosaic) == len(ucount[allele_index:]):
                    applied_filter.low_support_global.append(locus)
                    add_meta += ", LOW SUPP GLOBAL"
                continue
        # It stays here in case
        log.warn(f"LocusID = {locus}: ? ? {uval} | {n_seen_alleles}")
        applied_filter.unknown.append(locus)

    return applied_filter


def fill_spaces(uval_, ucount_):
    uval_ = list(uval_)
    new_uval, new_ucount = [], []
    for c in range(0, max(uval_)+1):
        if c in uval_:
            new_uval.append(np.int64(c))
            new_ucount.append(ucount_[uval_.index(c)])
        else:
            new_uval.append(np.int64(c))
            new_ucount.append(np.int64(0))
    return np.array(new_uval), np.array(new_ucount)


def load_data(config: TRConfig):
    tr_data = TRPaths(config)
    tr_data.load_read_files()
    samples_by_donor = {}
    tr_data.read_merge_path.sort()

    for f in tr_data.read_merge_path:
        donor_name = f.split("_")[0]
        if donor_name not in samples_by_donor:
            if f.endswith("tsv.gz"):
                samples_by_donor[donor_name] = os.path.join(tr_data.tr_read_lens, f)
        else:
            if f.endswith("tsv.gz"):
                log.warn(f'Multiple files for donor {donor_name}: {f}')

    log.info(str(samples_by_donor))
    assert len(samples_by_donor) > 0, f'Nothing found'
    return tr_data, samples_by_donor


def get_donors(donors_file_list):
    return open(donors_file_list).read().splitlines()


def summary_candidates(filtered_candidates: TRFilters, candidate_pq: str, include_germline: bool = False, bandwidth: int = 20):
    do_once = True
    uid = []
    labels = []
    loci = []
    asizes = []
    sample_names = []
    filters_list = []
    read_tissue_sum = {}
    cand_list = []
    results: pl.DataFrame = None
    for cand_type in filtered_candidates.keys():
        log.info(f'Candidate type: {cand_type}')
        cand_filter = filtered_candidates.filters(cand_type)
        if "GERMLINE" == cand_filter and not include_germline:
            continue
        candidate_loci = filtered_candidates.get(cand_type)
        if  0 == len(candidate_loci):
            continue
        tr_reads_lazy = pl.scan_parquet(candidate_pq).filter(pl.col("LocusID").is_in(candidate_loci))
        tr_candidates = tr_reads_lazy.collect(engine="streaming")
        if 0 == len(sample_names):
            sample_names = sorted(tr_candidates["SampleName"].unique().to_list())
        for lid in candidate_loci:
            if lid in filtered_candidates.low_support_global:
                continue
            tr_locus = tr_candidates.filter(pl.col("LocusID") == lid)
            l = tr_locus["ReadLen"].to_numpy()
            if 0 == len(l):
                continue
            ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
            ms.fit(l.reshape(-1, 1))
            cluster_centers = ms.cluster_centers_
            tr_locus = tr_locus.with_columns(pl.Series("Label", ms.labels_))
            l_groups = tr_locus["Label"].unique().to_list()
            ncol = len(sample_names)
            if max(l_groups) >= len(l_groups):
                l_groups, _ = fill_spaces(l_groups, [0]*len(l_groups))
            nrow = len(l_groups)
            uid.extend([f'{lid}_{l}' for l in l_groups])
            labels.extend(l_groups)
            asizes.extend([np.round(cluster_centers.tolist()[l][0]) for l in l_groups])
            loci.extend([lid] * nrow)
            filters_list.extend([cand_filter] * nrow)
            cand_list.extend([cand_type] * nrow)
            data = {f"col_{i}": [0] * nrow for i in range(ncol)}
            tmp_res = pl.DataFrame(data)
            tmp_res.columns = sample_names
            read_tissue = []
            for sname in sample_names:
                labs, lcount = np.unique(tr_locus.filter(pl.col("SampleName") == sname)["Label"].to_list(), return_counts=True)
                read_tissue.append(sum(lcount))
                for l, c in zip(labs, lcount):
                    tmp_res[int(l), sname] = c
            if do_once:
                results = tmp_res
                do_once = False
            else:
                results = pl.concat([results, tmp_res])
            read_tissue_sum[lid] = ",".join(str(x) for x in read_tissue)

    if results is not None:
        results = results.select(
            pl.sum_horizontal(pl.all()).alias("total"),
            pl.all()
        )
        results.insert_column(0, pl.Series("Filter", list(filters_list)))
        results.insert_column(0, pl.Series("CandidateClass", list(cand_list)))
        results.insert_column(0, pl.Series("AlleleSize", list(asizes)))
        results.insert_column(0, pl.Series("Labels", list(labels)))
        results.insert_column(0, pl.Series("LocusID", list(loci)))
        results.insert_column(0, pl.Series("UID", list(uid)))

    return results, read_tissue_sum


def vcf_header(tnames: list, contig_info: list):
    tnames_str = "\t".join(tnames)
    the_header = ["##fileformat=VCFv4.2", "##source=TR_pipeline v0.1"]
    for contig, contig_len in contig_info:
        the_header.append(f'##contig=<ID={contig},length={contig_len}>')
    the_header.append('##INFO=<ID=LOCUSID,Number=1,Type=Integer,Description="TR ID from catalog">')
    the_header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="TR allele length">')
    the_header.append('##INFO=<ID=MOTIF,Number=1,Type=String,Description="TR Motif">')
    the_header.append('##INFO=<ID=MOTIF_SIZE,Number=1,Type=Integer,Description="Size of the TR Motif">')
    the_header.append('##INFO=<ID=CAND_CLASS,Number=1,Type=String,Description="TR candidate classification">')
    the_header.append('##FILTER=<ID=PASS,Description="Filters passed">')
    the_header.append('##FILTER=<ID=GERMLINE,Description="It is considered a germline variant">')
    the_header.append('##FILTER=<ID=UNKNOWN,Description="Failed filter assignment">')
    the_header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    the_header.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth: ALT, REF">')
    the_header.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{tnames_str}")
    return "\n".join(the_header)


def vcf_maker(results_pl: pl.DataFrame, read_tissue_sum: dict, locus_catalog: dict, tr_data: TRPaths, each_donor: str, result_path: str = "."):
    tissue_names = results_pl.columns[7:]
    qual = "."
    # vcf_lines = []
    vcf_out = open(f'{result_path}/{each_donor}_results_tr.vcf', 'w')
    vcf_out.write(vcf_header(tissue_names, tr_data.adotto_v2_conitgs_wlen))
    vcf_out.write("\n")
    for row in results_pl.sort(by="LocusID").iter_rows(named=True):
        lid = row["LocusID"]
        tissue_sums = [int(x) for x in read_tissue_sum[lid].split(",")]
        vcf_entry = []
        chrom, pos, _, motif = locus_catalog[f'{lid}']
        vcf_entry.append(chrom)  # CHROM
        vcf_entry.append(str(pos))    # POS
        motif_size = len(motif)
        lab = row["Labels"]
        trid = f'{lid}.{lab}'
        tr_filter = row["Filter"]
        vcf_entry.append(trid)       # ID
        vcf_entry.append("N")        # REF
        vcf_entry.append("N")        # ALT
        vcf_entry.append(qual)        # QUAL
        vcf_entry.append(tr_filter)  # FILTER
        svlen = row["AlleleSize"]
        cand_class = row["CandidateClass"]
        vcf_entry.append(f'LOCUSID={lid};MOTIF={motif};MOTIF_SIZE={motif_size};SVLEN={svlen};CAND_CLASS={cand_class}')    # INFO
        vcf_entry.append(f'GT:AD')    # FORMAT(GT:AD)
        tc = 0
        for tissue in tissue_names:
            alt = row[tissue]
            tot = tissue_sums[tc]
            tc += 1
            ref = tot - alt
            gt = "./."
            if tot > 0:
                altp = alt/float(tot)
                gt = "0/0" if 0.15 > altp else "1/1" if altp > 0.85 else "0/1"
            vcf_entry.append(f'{gt}:{ref},{alt}')
        vcf_line = "\t".join(vcf_entry)
        vcf_out.write(f'{vcf_line}\n')

    log.info("Count finish")
    vcf_out.close()
    return True


def keys_to_int(d):
    return {int(k): v for k, v in d.items()}


def main(tr_config: TRConfig):
    tr_data, samples_by_donor = load_data(tr_config)
    # pre-concat
    results_dir = f'{config.res}/{config.date}' if "" != config.date else config.res

    donors_list = get_donors(tr_config.donors)
    with open(tr_config.donors_meta, 'r', encoding='utf-8') as file:
        sample_IDs = json.load(file)

    # By LocusID
    log.info(str(tr_data.adotto_v2_ff_tr_ids))
    # Smaller chunk size makes analysis faster
    divide_by = 100
    chunk_size = int(tr_data.adotto_v2_ff_tr_ids / divide_by)
    # results
    # donor_results = {d: "" for d in donors_list}

    for each_donor in donors_list:
        ##################################################
        # STEP 1
        # Candidate selection
        t = DevTiming(f'{each_donor}, step1: candidate search')
        t.start()
        output1 = f'{results_dir}/{each_donor}_step1-candidates_locusID.txt'
        cand_o1 = f'{results_dir}/{each_donor}_step1-candidates_stats.pkl'
        log.info(f'TR analysis of donor: {each_donor}')
        candidates = []
        candidates_stats = {}
        chunk_counter = 0
        if os.path.isfile(output1):
            log.info(f'"Output file {output1} exist, skipping..."')
            candidates = [int(lid) for lid in open(output1).read().splitlines()]
            with open(cand_o1, 'rb') as fin:
                # Parse the Pickle file directly
                candidates_stats = pickle.load(fin)
        elif os.path.isfile(tr_config.loci):
            candidates = [int(lid) for lid in open(tr_config.loci).read().splitlines()]
        else:
            for chunk in range(1, tr_data.adotto_v2_ff_tr_ids, chunk_size):
                chunk_counter += 1
                log.info(f"Loading data in chunks IDs: {chunk} - {chunk + chunk_size} ({chunk_counter}/{divide_by})")
                tr_reads_chunks = (
                     pl.scan_csv(samples_by_donor[each_donor], separator='\t', has_header=False, schema=tr_data.complete_schema)
                    .filter((pl.col("LocusID") >= chunk) & (pl.col("LocusID") < chunk + chunk_size))
                )
                tr_reads = tr_reads_chunks.collect(engine="streaming")
                log.info(f"Candidate selection start...")
                loci_list = tr_reads["LocusID"].unique().to_list()
                n_loci_select = len(loci_list)
                if 0 == n_loci_select:
                    continue
                loci_chunk_size = max(1, int(n_loci_select /10))
                if n_loci_select > 0 and loci_chunk_size > 0:
                    for chunkl in range(0, n_loci_select, loci_chunk_size):
                        these_loci = loci_list[chunkl:chunkl+loci_chunk_size]
                        cand, cand_stats = find_candidates(tr_reads.filter(pl.col("LocusID").is_in(these_loci)))
                        candidates.extend(cand)
                        candidates_stats.update(cand_stats)
                log.info(f"Candidates added to list, total: {len(candidates)}")

            # Save candidate LocusIDs as text
            log.info(f'Saving candidates results for sample {each_donor}')
            with open(f"{output1}", "w") as fout:
                for cand in candidates:
                    fout.write(f'{cand}\n')
            fout.close()
            with open(cand_o1, "wb") as fout:
                pickle.dump(candidates_stats, fout)
            fout.close()

        t.end(True)
        ##################################################

        ##################################################
        # STEP 2
        # Candidate cleaning
        t = DevTiming(f'{each_donor}, step2: update table')
        t.start()

        if len(candidates) == 0:
            lid_candidate_list = [int(lid) for lid in open(output1).read().splitlines()]
        else:
            lid_candidate_list = candidates

        tr_step2_tsv = f'{results_dir}/{each_donor}_step2-candidates.tsv'
        output2 = f'{results_dir}/{each_donor}_step2-candidates.pq'
        if os.path.isfile(output2):
            log.info(f'Output file "{output2}" exist, skipping...')
            tr_reads = pl.read_parquet(output2)
        else:
            # By LocusID load all
            tr_reads_chunks = (
                pl.scan_csv(samples_by_donor[each_donor], separator='\t', has_header=False, schema=tr_data.complete_schema)
                .filter((pl.col("LocusID").is_in(lid_candidate_list)))
            )
            tr_reads = tr_reads_chunks.collect(engine="streaming")

            # re-save?
            sample_names = []
            for k in tr_reads.with_columns(
                SampleID = pl.concat_str(["Donor", "Tissue", "Core"], separator="-")
            )["SampleID"].to_list():
                sample_names.append(sample_IDs[k])

            # Save candidate table
            tr_reads = tr_reads.with_columns(pl.Series("SampleName", sample_names))
            tr_reads.write_parquet(output2)
            tr_reads.write_csv(tr_step2_tsv, separator="\t", include_header=True)
        t.end(True)
        ##################################################

        ##################################################
        # STEP 3
        # Candidate filtering
        t = DevTiming(f'{each_donor}, step3: candidate filter')
        t.start()
        output3 = f'{results_dir}/{each_donor}_step3-candidates_class.json'
        filtered_candidates = TRFilters()
        if os.path.isfile(output3):
            log.info(f'Output file "{output3}" exist, skipping...')
            tr_reads = pl.read_parquet(output2)
        else:
            min_reads = tr_reads.width*config.min_reads_sample
            # Candidate filtering also by chunks
            loci_list = tr_reads["LocusID"].unique().to_list()
            n_loci_cand = len(loci_list)
            cand_chunk_size = max(int(n_loci_cand /100), 1)
            for chunk in range(0, n_loci_cand, cand_chunk_size):
                this_loci = loci_list[chunk:chunk+cand_chunk_size]
                filtered_candidates.update(filter_candidates(tr_reads.filter(pl.col("LocusID").is_in(this_loci)), min_reads, config.noqc, candidates_stats))
            log.info(f'Donor {each_donor} candidate results:\n{filtered_candidates}')
            filtered_candidates.jsonfy(output3)
        t.end(True)
        ##################################################

        ##################################################
        # STEP 4
        # Make results table
        t = DevTiming(f'{each_donor}, step4: candidate summary')
        t.start()
        output4 = f'{results_dir}/{each_donor}_step4-results_table_tissue.pq'
        output42 = f'{results_dir}/{each_donor}_step4-reads_per_tissue_sums.json'
        if os.path.isfile(output4) and os.path.isfile(output42):
            log.info(f'Output file "{output4}" exist, skipping...')
            results = pl.read_parquet(output4)
            with open(output42, 'r') as fin:
                # Parse the JSON file directly
                read_tissue_sum = json.load(fin, object_hook=keys_to_int)
            fin.close()
        else:
            filtered_candidates.from_json(output3)
            results, read_tissue_sum = summary_candidates(filtered_candidates, output2, True)
            # Save output as pq
            if results is not None:
                results.sort(by="LocusID").write_parquet(output4)
            else:
                log.warn("No results writen")
            if len(read_tissue_sum) > 0:
                with open(output42, "w") as fout:
                    json.dump(read_tissue_sum, fout, indent=4)
                fout.close()
            else:
                log.warn("No results writen2")
        t.end()
        ##################################################

        ##################################################
        # STEP 5
        # catalog locus data
        t = DevTiming(f'{each_donor}, step5: VCF output')
        t.start()
        locus_catalog = load_adotto_bed(tr_data.adotto_v2_path)
        # Save output as VCF
        if vcf_maker(results, read_tissue_sum, locus_catalog, tr_data, each_donor, f'{results_dir}'):
            log.info("VCF writen")
        else:
            log.err("No results writen")
        t.end()
        output_done = f'{results_dir}/{each_donor}_steps.DONE'
        with open(output_done, "w") as _: pass



def parse_args() -> TRConfig:
    parser = argparse.ArgumentParser()
    parser.add_argument("--donors-list", required=True, type=str)          # donor name(s)/id(s)
    parser.add_argument("--donors-meta", required=True, type=str)          # donor metadata (donors.json)
    parser.add_argument("--data-dir", required=True, type=str)             # base directory
    parser.add_argument("--res-dir", required=True, type=str)              # results directory
    parser.add_argument("--tr-bed", required=True, type=str)               # bed file withthe TR catalog
    # extra
    parser.add_argument("--date", required=False, type=str, default="")                  # date for results directory
    parser.add_argument("--selected-loci", required=False, type=str, default="")  # use only these loci
    parser.add_argument("--noqc", required=False, action='store_true', default=False)  # do not remove germline TRs
    trconf = TRConfig()
    trconf.set(parser.parse_args())
    return trconf


if "__main__" == __name__:
    config = parse_args()
    assert os.path.isfile(config.donors), f'File {config.donors} does not exixist'
    assert os.path.isfile(config.tr_bed), f'File {config.tr_bed} does not exixist'
    assert os.path.isfile(config.donors_meta), f'File {config.donors_meta} does not exixist'
    assert os.path.isdir(config.data), f'Directory {config.data} does not exixist'
    if not os.path.isdir(config.res):
        log.info(f'Directory {config.res} does not exixist')
        try:
            os.makedirs(config.res, exist_ok=True)
        except:
            log.err("??")
    if "" != config.date and not os.path.isdir(f'{config.res}/{config.date}'):
        log.info(f'Directory {config.res}/{config.date} does not exixist')
        try:
            os.makedirs(f'{config.res}/{config.date}', exist_ok=True)
        except:
            log.err(f"cannot make dirs: {config.res}/{config.date}")
    main(config)
