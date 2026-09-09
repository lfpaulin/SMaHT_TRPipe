#!/usr/bin/env python
import os
import sys
import json
import gzip
import pickle
import polars as pl
import numpy as np
import argparse
from sklearn.cluster import MeanShift


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
    "B004": "PacBio Revio HiFi Fiber-seq",
    "D001": "ONT PromethION 24 WG",
    "D003": "ONT PromethION 24 Ultra-Long WGS"
}


class MyLog():
    def __init__(self) -> None:
        pass

    @staticmethod
    def out(text: str):
        sys.stdout.write(f'{text}\n')

    @staticmethod
    def info(text: str):
        sys.stdout.write(f'#INFO: {text}\n')

    @staticmethod
    def warn(text: str):
        sys.stdout.write(f'#WARNING: {text}\n')

    @staticmethod
    def err(text: str):
        sys.stderr.write(f'#ERROR: {text}\n')
        sys.exit(1)

    @staticmethod
    def dev(text: str):
        sys.stderr.write(f'#DEBUG: {text}\n')

log = MyLog()


class TRConfig:
    def __init__(self):
        self.candidate_loci: str = ""
        self.donors_meta: str = ""
        self.data: str = ""
        self.res: str = ""
        self.min_reads_sample = 2
        self.tr_pq_schema = {
            "Donor": pl.String,
            "Tissue": pl.String,
            "Core": pl.String,
            "Contig": pl.String,
            "Pos": pl.UInt32,
            "End": pl.UInt32,
            "LocusID": pl.UInt32,
            "HP": pl.UInt8,
            "ReadLen": pl.Int32,
            "SampleName": pl.String
        }

    def set(self, args: argparse.Namespace):
        self.candidate_loci = args.candidate_loci
        self.donors_meta = args.donors_meta
        self.data = args.data_dir
        self.res = args.res_dir

    def setu(self, arg_tuple: tuple[str,str,str,str]):
        # donors_list, donors_meta, data_dir, res_dir
        self.candidate_loci, self.donors_meta, self.data, self.res = arg_tuple


class TRPaths:
    def __init__(self, config: TRConfig):
        self.tr_parsed = config.data
        # these are TSV tables
        self.donors: list = []
        self.step1_stats: list = []
        self.step2_candidates: list = []
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
            "ReadLen": pl.Int32,
            "SampleName": pl.String
        }
        self.tr_loci_adotto = 1800070

    def load_files(self):
        files_in_path = [
            f for f in os.listdir(self.tr_parsed) if os.path.isfile(os.path.join(self.tr_parsed, f))
        ]
        assert len(files_in_path) > 0, f'No files found in {self.tr_parsed}: {files_in_path}'
        self.step1_stats = sorted([ f for f in files_in_path if "pkl" in f ])
        self.step2_candidates = sorted([ f for f in files_in_path if "pq" in f ])
        assert len(self.step1_stats) > 0, f'No files from "step1" found in {files_in_path}: {self.step1_stats}'
        assert len(self.step2_candidates) > 0, f'No files from "step2" found in {files_in_path}: {self.step2_candidates}'
        assert len(self.step1_stats) == len(self.step2_candidates), f'Same name of files needed in both "step1" and "step2": {len(self.step1_stats)} =/= {len(self.step2_candidates)}'
        self.donors = [d.split("_").pop(0) for d in self.step1_stats]


padding = {
    "trgt": 50,
    "medaka": 250
}


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


def load_data(config: TRConfig):
    tr_data = TRPaths(config)
    tr_data.load_files()
    return tr_data


def keys_to_int(d):
    return {int(k): v for k, v in d.items()}


def load_candidates(candidate_loci_file: str) -> dict[int, list]:
    if os.path.isfile(f'{candidate_loci_file}.json'):
        with open(f'{candidate_loci_file}.json', 'r', encoding='utf-8') as file:
            cand_loci = json.load(file, object_hook=keys_to_int)
    else:
        with open(candidate_loci_file, 'r', encoding='utf-8') as file:
            if "json" in candidate_loci_file:
                cand_loci = json.load(file, object_hook=keys_to_int)
            else:
                cand_loci = {}
                for line in file:
                    lid, count = line.rstrip("\n").split("\t")
                    if count not in cand_loci:
                        cand_loci[count] = [lid]
                    else:
                        cand_loci[count].append(lid)
                with open(f'{candidate_loci_file}.json', "w") as fout:
                    json.dump(cand_loci, fout, indent=4)
        file.close()
    return cand_loci


def clusteirng(allele_lens: np.ndarray, bandwidth: int, use_seeds):
    if len(allele_lens) == 0:
        return np.array([]), np.array([])
    ms = MeanShift(bandwidth=bandwidth, seeds=use_seeds)
    ms.fit(allele_lens)
    # Extract cluster labels and center coordinates
    uval_, ucount_ = np.unique(ms.labels_, return_counts=True)
    uvalv_ = [ms.cluster_centers_[v][0] for v in uval_]
    # ulabel_schema = {"Mean": pl.Float64, "Counts": pl.UInt16}
    # ulabel = pl.DataFrame({"Mean": uval, "Counts": ucount}, schema=ulabel_schema)
    # filters
    has_negative_len_ = sum([v < 0 for v in uvalv_]) > 0
    with_single_read_ = sum([c <= 2 for c in ucount_])

    # return ms.cluster_centers_, ms.labels_, uval_, ucount_, has_negative_len_, with_single_read_
    return ucount_, ms.cluster_centers_.reshape(1,-1)


def center_norm(cl_center:np.ndarray, mod:int = 13):
    return cl_center - (cl_center % mod)


def detect_smears(donors_list: list, candidate_loci: dict, donor_stats: dict) -> dict[int, tuple[int, int, int, int]]:
    t = DevTiming(f'Detect Smears: process/test')
    t.start()
    """
    The stats include the counts per cluster, we use them to select for 
    loci that have four or more alleles with min 2 reads, expexting 1-2 
    to be germline and the rest non-germline. 
    """
    min_clusters_2reads = 4
    min_reads_cluster = 2
    donor_range_low, donor_range_upp = (3, 10)
    smeary_loci = {}

    for ndonor in candidate_loci:
        if ndonor >= donor_range_low:
            for lid in candidate_loci[ndonor]:
                lid = int(lid)
                donor_count = []
                donor_count_strict = []
                donor_count_single = []
                clust_tech_bias = []
                for each_donor in donors_list:
                    try:
                        bw, cl_center, _, ucount, _ = donor_stats[each_donor][lid]
                        clst_over_min = sum(1 for c in ucount if c >= min_reads_cluster)
                        dcount = 1 if clst_over_min > min_clusters_2reads else 0
                        donor_count.append(dcount)
                        dcountstx = 1 if sum(1 for c in ucount if c >= min_reads_cluster*4) > min_clusters_2reads+2 else 0
                        donor_count_strict.append(dcountstx)
                        sgl = sum(1 for c in ucount if c == 1)
                        ncl = len(ucount)
                        dcountsgl = 1 if (clst_over_min > 2) or (clst_over_min >= 1 and sgl not in {1,2,3}) else 0  #BUG
                        donor_count_single.append(dcountsgl)
                        # if True: #sum(donor_count_strict) > 0:
                        #     uc = ", ".join([str(c) for c in ucount])
                        #     log.dev(f'{lid}|{each_donor}|{ndonor}|{dcount}|{dcountstx}|{dcountsgl}: {bw} | {uc}')
                    except KeyError:
                        # not all donors have the candidate, by design
                        pass
                if donor_range_low <= sum(donor_count_single) <= donor_range_upp or sum(donor_count_strict) > 0:
                    # log.out(lid, ndonor, sum(donor_count), sum(donor_count_single), sum(donor_count_strict))
                    smeary_loci[lid] = (ndonor, sum(donor_count), sum(donor_count_single), sum(donor_count_strict))
        else:
            continue
    t.end()
    return smeary_loci


def assess_smear_biased(reads: pl.DataFrame, each_donor: str, donor_stats: dict) -> dict[int, tuple[str, str]]:
    """
    No bias: value ==>> 0
        same name of cluseters
    Some bias: value ==>> 1
        dual unique, each tech shows one or more unique 
        cluster sizes, same number of cluseters from each 
        tech but less than the overall
    Clear bias ONT: value ==>> 2
        more clusters seen on ONT compared to PB
    Clear bias PB: value ==>> 3
        more clusters seen on PB compared to ONT
    Missed pattern: value ==>> -1

    We assess if the number of discrepancies/biases is in a single tissue, 
    or there is an abundandance across the tissues
    """
    t = DevTiming(f'Assess Smear Biased: process/test')
    t.start()
    tech_bias = {}
    for lid in reads["LocusID"].unique():
        bw, cl_center, _, ucount, _ = donor_stats[each_donor][lid]
        if len(ucount) <= 3:
            clust_tech_bias = 91
            print(lid, clust_tech_bias, ucount)
            tech_bias[lid] = (each_donor, clust_tech_bias)
            continue
        
        rlen_pb = reads.filter(  (pl.col("Core").str.contains("B00")) & (pl.col("LocusID") == lid) )["ReadLen"].to_numpy()
        rlen_ont = reads.filter( (pl.col("Core").str.contains("D00")) & (pl.col("LocusID") == lid) )["ReadLen"].to_numpy()

        if len(rlen_pb) < 50 or len(rlen_ont) < 50:
            clust_tech_bias = 90
            print(lid, clust_tech_bias, ucount, "|", rlen_pb.tolist(), rlen_ont.tolist())
            continue

        rlen_pb = rlen_pb.reshape(-1, 1)
        rlen_ont = rlen_ont.reshape(-1, 1)

        clust_tech_bias = 0
        cluster_pb_counts, cluster_pb_centers = clusteirng(rlen_pb, bw, cl_center)
        cluster_ont_counts, cluster_ont_centers = clusteirng(rlen_ont, bw, cl_center)


        if len(cluster_pb_centers) == 0 or len(cluster_ont_centers) == 0:
            clust_tech_bias = 92
            print(lid, clust_tech_bias, ucount)
            tech_bias[lid] = (each_donor, clust_tech_bias)
            continue

        cldon = [int(cl) for cl in center_norm(cl_center.reshape(1,-1), bw)[0].tolist()  ]
        clpcb = [int(cl) for cl in center_norm(cluster_pb_centers, bw)[0].tolist() ]
        clont = [int(cl) for cl in center_norm(cluster_ont_centers, bw)[0].tolist()]

        # log.dev(f'{each_donor}\n    D: {cldon}\n    P: {clpcb}\n    O: {clont}')
        if len(cldon) == len(clpcb) == len(clont):
            clust_tech_bias = 0  # same size => same values (?)
        elif len(cldon) > len(clpcb) == len(clont):
            clust_tech_bias = 1  # same size tech => NOT same values
        elif len(cldon) >= len(clont) > len(clpcb):
            clust_tech_bias = 2  # bias ONT
        elif len(cldon) >= len(clpcb) > len(clont):
            clust_tech_bias = 3  # bias PB
        else:
            clust_tech_bias = -1  # wtf
            log.warn(f'{lid} {each_donor} -1: {cldon} {clpcb} {clont}')

        # bias_id, counts_id = np.unique(clust_tech_bias, return_counts=True)
        print(lid, clust_tech_bias, cldon, clpcb, clont)
        tech_bias[lid] = (each_donor, clust_tech_bias)

    t.end()
    return tech_bias


def main(tr_config: TRConfig):
    t = DevTiming(f'Main: load data')
    t.start()
    tr_data = load_data(tr_config)
    t.end()

    # candidate loci
    t = DevTiming(f'Main: load candidates')
    t.start()
    candidate_loci = load_candidates(tr_config.candidate_loci)
    t.end()


    t = DevTiming(f'Main: read pickle')
    t.start()
    donor_stats = {}
    for donor, pkl in zip(tr_data.donors, tr_data.step1_stats):
        with open(pkl, 'rb') as fin:
            # Parse the Pickle file directly
            donor_stats[donor] = pickle.load(fin)
    t.end()

    smeary_loci = detect_smears(tr_data.donors, candidate_loci, donor_stats)
    
    t = DevTiming(f'Main: read pqs')
    t.start()
    # donor_qps = { donor: pq for donor, pq in zip(tr_data.donors, tr_data.step2_candidates) }
    donor_qps = { donor: pl.scan_parquet(pq, schema=tr_data.complete_schema) for donor, pq in zip(tr_data.donors, tr_data.step2_candidates) }
    t.end()


    reads: pl.DataFrame
    # divide_by = 100
    # chunk_counter = 0
    # chunk_size = int(tr_data.tr_loci_adotto/divide_by)
    tech_bias_lids = {lid: [] for lid in smeary_loci.keys()}
    keep_cols = ["Donor", "Core", "LocusID", "HP", "ReadLen"]
    for each_donor in tr_data.donors:
        log.info(f'{each_donor}')
        reads = donor_qps[each_donor].collect(engine="streaming").filter(pl.col("LocusID").is_in(smeary_loci.keys()))[keep_cols]
        tmp = assess_smear_biased(reads, each_donor, donor_stats)
        try:
            for lid in smeary_loci:
                tech_bias_lids[lid].append(tmp[lid])
        except KeyError:
            # not all donors have the candidate, by design
            pass

    with open("test-bias_tech.pkl", "wb") as fout:
        pickle.dump(tech_bias_lids, fout)
    fout.close()
    
    # for lid in smeary_loci:
    #     tech_bias = tech_bias_lids[lid]
    #     ndonor, donor_count, donor_count_single, donor_count_strict = smeary_loci
    #     print(lid, donor_count, ndonor, donor_count_strict)





def parse_args() -> TRConfig:
    parser = argparse.ArgumentParser()
    # required
    parser.add_argument('-c', "--candidate-loci", required=True, type=str)  # use only these loci
    parser.add_argument('-d', "--data-dir", required=True, type=str)        # base directory
    parser.add_argument('-r', "--res-dir", required=True, type=str)               # results directory
    # extra
    parser.add_argument('-m', "--donors-meta", required=False, type=str)     # donor metadata (donors.json)
    # parse
    trconf = TRConfig()
    trconf.set(parser.parse_args())
    return trconf


if "__main__" == __name__:
    config = parse_args()
    assert os.path.isfile(config.candidate_loci), f'File {config.candidate_loci} does not exixist'
    assert os.path.isfile(config.donors_meta), f'File {config.donors_meta} does not exixist'
    assert os.path.isdir(config.data), f'Directory {config.data} does not exixist'
    if not os.path.isdir(config.res):
        log.info(f'Directory {config.res} does not exixist')
        try:
            os.makedirs(config.res, exist_ok=True)
        except:
            log.err("mkdir failed")
    main(config)
