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


class MyLog:
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
        self.data: str = ""
        self.tr_bed: str = ""
        self.min_reads_sample = 5

    def set(self, args: argparse.Namespace):
        self.donors = args.donors_list
        self.data = args.data_dir
        self.tr_bed = args.tr_bed


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

    def load_adotto_bed(self):
        tr_catalog = {}
        bed_file = gzip.open(self.adotto_v2_path, "rt") if "gz" in self.adotto_v2_path \
            else open(self.adotto_v2_path, "r")
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


def clustering(tr_reads_df: pl.DataFrame, bandwidth: int = 20):
    allele_lens = tr_reads_df["ReadLen"].to_numpy().reshape(-1, 1)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(allele_lens)
    # Extract cluster labels and center coordinates
    lab_index, lab_count = np.unique(ms.labels_, return_counts=True)
    return len(lab_index) <= 2
    lab_values = np.array([ms.cluster_centers_[v][0] for v in lab_index])
    with_single_read_ = sum([c <= 2 for c in lab_count])

    bias_status = ""
    if len(lab_index) <= 3:
        bias_status = "TBD"
    else:
        rlen_pb = tr_reads_df.filter(pl.col("Core").str.contains("B00"))["ReadLen"].to_numpy().reshape(-1, 1)
        rlen_ont = tr_reads_df.filter(pl.col("Core").str.contains("D00"))["ReadLen"].to_numpy().reshape(-1, 1)
        if len(rlen_pb) == 0 or len(rlen_ont) == 0:
            bias_status = "NA_PB" if 0 == len(rlen_pb) else "NA_ONT"
        elif len(rlen_pb) < 50 and len(rlen_ont) < 50:
            bias_status = "LOW"
        elif 0 < len(rlen_pb) < 50:
            bias_status = "LOW_PB"
        elif 0 < len(rlen_ont) < 50:
            bias_status = "LOW_ONT"
        else:
            # Sometimes Mean-shift get a single value instead of array, although not used provided,
            # We need to catch that error
            bias_statusl = []
            try:
                mspb = MeanShift(bandwidth=bandwidth, seeds=lab_values.reshape(-1, 1))
                mspb.fit(rlen_pb)
                n_lab_pb = len(np.unique(mspb.labels_))
            except ValueError as epb:
                bias_statusl.append("MS_ERR_PB")
            try:
                msont = MeanShift(bandwidth=bandwidth, seeds=lab_values.reshape(-1, 1))
                msont.fit(rlen_ont)
                n_lab_ont = len(np.unique(msont.labels_))
            except ValueError as eont:
                bias_statusl.append("MS_ERR_ONT")

            if len(bias_statusl) == 0:
                n_lab_donor = len(lab_index)
                if n_lab_donor == n_lab_pb == n_lab_ont:
                    bias_status = "NONE0"  # same size => same values (?)
                elif n_lab_donor > n_lab_pb == n_lab_ont:
                    bias_status = "NONE1"  # same size tech => NOT same values i.e  A,B,C,D     A,B,C     A,B,D
                elif n_lab_donor >= n_lab_ont > n_lab_pb:
                    bias_status = "ONT"  # bias ONT  A,B,C,D,E     A,C,D,E     A,B
                elif n_lab_donor >= n_lab_pb > n_lab_ont:
                    bias_status = "PB"  # bias PB
                else:
                    pass  # wtf
            else:
                bias_status = ",".join(bias_statusl)  # error in mean-shift that we sort of cannot control

    return ms.cluster_centers_, ms.labels_, lab_index, lab_values, lab_count, has_negative_len_, with_single_read_, bias_status


def remove_germline(tr_reads: pl.DataFrame, bandwidth = 20):
    germline = []
    non_germline = []
    loci_list = tr_reads["LocusID"].unique().to_list()
    for locus in loci_list:
        allele_lens = tr_reads.filter(pl.col("LocusID") == locus)["ReadLen"].to_numpy().reshape(-1, 1)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(allele_lens)
        # Extract cluster labels and count them
        lab_index, lab_count = np.unique(ms.labels_, return_counts=True)
        if len(lab_index) <= 2:
            germline.append(locus)
        else:
            non_germline.append(locus)
    return germline, non_germline


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
        # seq-tech bias
        self.tech_bias = []

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
        # tech bias
        tech_bias = len(self.tech_bias)
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
            f"Low support global: {low_support_global}\n"
            f"Likely tech bias: {tech_bias}"
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
        self.tech_bias.extend(tr_candidate_list.tech_bias)

    def filters(self, tr_class):
        filter_pass = {"single_allele_smear", "mosaic_allele", "third_allele_mosaic", "smear"}
        if hasattr(self, tr_class):
            if tr_class in filter_pass:
                return "PASS"
            elif tr_class == "unknown":
                return "UNKNOWN"
            elif tr_class == "tech_bias":
                return "TECH_BIAS"
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


def main(tr_config: TRConfig):
    tr_data, samples_by_donor = load_data(tr_config)
    # pre-concat
    donors_list = get_donors(tr_config.donors)

    # By LocusID
    log.info(str(tr_data.adotto_v2_ff_tr_ids))
    # Smaller chunk size makes analysis faster
    divide_by = 40
    chunk_size = int(tr_data.adotto_v2_ff_tr_ids / divide_by)
    # results
    # donor_results = {d: "" for d in donors_list}

    for each_donor in donors_list:
        ##################################################
        # STEP 1
        # Candidate selection
        t = DevTiming(f'{each_donor}, step1: candidate search')
        t.start()
        chunk_counter = 0
        germline = []
        not_germline = []
        for chunk in range(1, tr_data.adotto_v2_ff_tr_ids, chunk_size):
            chunk_counter += 1
            log.info(f"Loading data in chunks IDs: {chunk} - {chunk + chunk_size} ({chunk_counter}/{divide_by})")
            t2 = DevTiming(f'{each_donor}, step1.1: data loader')
            t2.start()
            tr_reads_chunks = (
                 pl.scan_csv(samples_by_donor[each_donor], separator='\t', has_header=False, schema=tr_data.complete_schema)
                .filter((pl.col("LocusID") >= chunk) & (pl.col("LocusID") < chunk + chunk_size))
            )
            tr_reads = tr_reads_chunks.collect(engine="streaming")
            t2.end()
            loci_list = tr_reads["LocusID"].unique().to_list()
            n_loci_select = len(loci_list)
            if 0 == n_loci_select:
                continue
            loci_chunk_size = max(1, int(n_loci_select /10))
            if n_loci_select > 0 and loci_chunk_size > 0:
                for chunkl in range(0, n_loci_select, loci_chunk_size):
                    t3 = DevTiming(f'{each_donor}, step1.2: chunk loader {chunkl}')
                    t3.start()
                    these_loci = loci_list[chunkl:chunkl+loci_chunk_size]
                    germ, not_germ = remove_germline(tr_reads.filter(pl.col("LocusID").is_in(these_loci)))
                    germline.extend(germ)
                    not_germline.extend(not_germ)
                    t3.end()
            log.info(f"Germline added to list, total: {len(germline)}")
            print(f'DATA: germline {germline}')
            log.info(f"NOT germline added to list, total: {len(not_germline)}")
            print(f'DATA: non-germ: {not_germline}')

        t.end()


def parse_args() -> TRConfig:
    parser = argparse.ArgumentParser()
    parser.add_argument("--donors-list", required=True, type=str)          # donor name(s)/id(s)
    parser.add_argument("--data-dir", required=True, type=str)             # base directory
    parser.add_argument("--tr-bed", required=True, type=str)               # bed file with the TR catalog
    # extra
    trconf = TRConfig()
    trconf.set(parser.parse_args())
    return trconf


if "__main__" == __name__:
    config = parse_args()
    assert os.path.isfile(config.donors), f'File {config.donors} does not exist'
    assert os.path.isfile(config.tr_bed), f'File {config.tr_bed} does not exist'
    assert os.path.isdir(config.data), f'Directory {config.data} does not exist'
    main(config)
