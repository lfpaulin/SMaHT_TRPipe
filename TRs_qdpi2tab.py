#!/usr/bin/env python
import os
import sys
import json
import gzip
import pysam
import polars as pl
import numpy as np
import argparse
from TRs_polars_allele_sizes_from_reads import main as TRAlleleSize
from TRs_polars_allele_sizes_from_reads import TRConfig

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


class TRQDPIConfig:
    def __init__(self, args: argparse.Namespace ):
        self.donor = args.donor
        self.qdpi_list = args.qdpi_list
        self.tr_bed = args.tr_bed
        self.output = args.output


def load_adotto_bed(bed_path: str, use_loci: list[int] = []):
    tr_catalog = {}
    bed_file = gzip.open(bed_path, "rt") if "gz" in bed_path else open(bed_path, "r")
    for line in bed_file:
        contig, pos, end, info = line.rstrip("\n").split("\t")
        trid, motif, _ = info.split(";")
        _, trid = trid.split("=")
        _, motif = motif.split("=")
        tr_catalog[trid] = (contig, int(pos), int(end), motif)
    bed_file.close()
    return tr_catalog


def get_files(qdpi_file_list):
    qdpi_list = []
    for qdpi_f in open(qdpi_file_list).read().splitlines():
        if qdpi_f.endswith("gz") and os.path.isfile(qdpi_f):
            qdpi_list.append(qdpi_f)
        else:
            log.warn(f'File {qdpi_f} does not exisit, skippig')
    assert len(qdpi_list) > 0, f'No files found'
    return qdpi_list


def get_sample_names(file_list):
    # expected SMHT042-3AD-002A2-M79-D001-bcm-SMAFILFO2FXV-(...).01_GRCh38.aligned.sorted.dacHP.bed.gz
    fields = 6
    return ["-".join(os.path.basename(f).split("-")[:fields]) for f in file_list]


def load_data(config: TRQDPIConfig):
    files_list = get_files(config.qdpi_list)
    sample_list = get_sample_names(files_list)
    
    if len(sample_list) != len(set(sample_list)):
        log.warn(f'Multiple files the a sample')

    return sample_list, files_list


def parse_rlen(hp0, hp1, hp2, tr_len_base):
    rln0, rln1, rln2 = [], [], []
    if hp0 != ".":
        rln0 = [int(r)+tr_len_base for r in hp0.split(",")]
    if hp1 != ".":
        rln1 = [int(r)+tr_len_base for r in hp1.split(",")]
    if hp2 != ".":
        rln2 = [int(r)+tr_len_base for r in hp2.split(",")]
    return rln0, rln1, rln2


def main(config: TRQDPIConfig):
    # get files
    names_by_tissue, samples_by_tissue = load_data(config)
    # get handlers
    handlers_by_tissue = [pysam.TabixFile(f) for f in samples_by_tissue]

    tr_bed = gzip.open(config.tr_bed, "rt")
    out = gzip.open(config.output, "wt")
    for bed_entry in tr_bed:
        chrom, pos, end, info = bed_entry.split("\t")
        region = f'{chrom}:{pos}-{end}'
        lid = info.split(";")[0].split("=")[1]
        tr_len_base = int(end) - int(pos)
        
        for sname, qdpi_handl in zip(names_by_tissue, handlers_by_tissue):
            donor, tissue, core ,_, tech, _ = sname.split("-")
            for qdpi_entry in qdpi_handl.fetch(region=region):
                chrom, pos, end, cov, hp0, hp1, hp2 = qdpi_entry.split("\t")
                hp0_rl, hp1_rl, hp2_rl = parse_rlen(hp0, hp1, hp2, tr_len_base)
                hap = 0
                for rl in hp0_rl:
                    out.write(f'{donor}\t{tissue}\t{core}|{tech}\t{chrom}\t{pos}\t{end}\t{lid}\t{hap}\t{rl}\n')
                hap = 1
                for rl in hp1_rl:
                    out.write(f'{donor}\t{tissue}\t{core}|{tech}\t{chrom}\t{pos}\t{end}\t{lid}\t{hap}\t{rl}\n')
                hap = 2
                for rl in hp2_rl:
                    out.write(f'{donor}\t{tissue}\t{core}|{tech}\t{chrom}\t{pos}\t{end}\t{lid}\t{hap}\t{rl}\n')

    # closing
    _ = [fh.close() for fh in handlers_by_tissue]
    out.close()
    # done
    output_done = f'{config.output}.DONE'
    with open(output_done, "w") as _: pass



    # some_args = ("",)
    # tr_config = TRConfig()
    # tr_config.setu(some_args)
    # TRAlleleSize(tr_config)

    # donors_list = get_donors(tr_config.donors)
    # with open(tr_config.donors_meta, 'r', encoding='utf-8') as file:
    #     sample_IDs = json.load(file)

    # # By LocusID
    # log.info(str(tr_data.adotto_v2_ff_tr_ids))
    # # Smaller chunk size makes analysis faster
    # chunk_size = int(tr_data.adotto_v2_ff_tr_ids / 100)
    # # results
    # # donor_results = {d: "" for d in donors_list}

    # for each_donor in donors_list:
    #     ##################################################
    #     # STEP 1
    #     # Candidate selection
    #     t = DevTiming(f'{each_donor}, step1: candidate search')
    #     t.start()
    #     output1 = f'{results_dir}/{each_donor}_step1-candidates_locusID.txt'
    #     log.info(f'TR analysis of donor: {each_donor}')
    #     candidates = []
    #     if os.path.isfile(output1):
    #         log.info(f'"Output file {output1} exist, skipping..."')
    #         candidates = [int(lid) for lid in open(output1).read().splitlines()]
    #     elif os.path.isfile(tr_config.loci):
    #         candidates = [int(lid) for lid in open(tr_config.loci).read().splitlines()]
    #     else:
    #         for chunk in range(1, tr_data.adotto_v2_ff_tr_ids, chunk_size):
    #             log.info(f"Loading data in chunks IDs: {chunk} - {chunk + chunk_size}")
    #             tr_reads_chunks = (
    #                  pl.scan_csv(samples_by_donor[each_donor], separator='\t', has_header=False, schema=tr_data.complete_schema)
    #                 .filter((pl.col("LocusID") >= chunk) & (pl.col("LocusID") < chunk + chunk_size))
    #             )
    #             tr_reads = tr_reads_chunks.collect(engine="streaming")
    #             log.info(f"Candidate selection start...")
    #             loci_list = tr_reads["LocusID"].unique().to_list()
    #             n_loci_select = len(loci_list)
    #             if 0 == n_loci_select:
    #                 continue
    #             loci_chunk_size = max(1, int(n_loci_select /10))
    #             if n_loci_select > 0 and loci_chunk_size > 0:
    #                 for chunkl in range(0, n_loci_select, loci_chunk_size):
    #                     these_loci = loci_list[chunkl:chunkl+loci_chunk_size]
    #                     cand = find_candidates(tr_reads.filter(pl.col("LocusID").is_in(these_loci)))
    #                     candidates.extend(cand)
    #             log.info(f"Candidates added to list, total: {len(candidates)}")

    #         # Save candidate LocusIDs as text
    #         log.info(f'Saving candidates results for sample {each_donor}')
    #         with open(f"{output1}", "w") as fout:
    #             for cand in candidates:
    #                 fout.write(f'{cand}\n')
    #         fout.close()
    #     t.end(True)
    #     ##################################################

    #     ##################################################
    #     # STEP 2
    #     # Candidate cleaning
    #     t = DevTiming(f'{each_donor}, step2: update table')
    #     t.start()

    #     if len(candidates) == 0:
    #         lid_candidate_list = [int(lid) for lid in open(output1).read().splitlines()]
    #     else:
    #         lid_candidate_list = candidates

    #     output2 = f'{results_dir}/{each_donor}_step2-candidates.tsv'
    #     tr_results_pq = f'{results_dir}/{each_donor}_step2-candidates.pq'
    #     if os.path.isfile(output2):
    #         log.info(f'Output file "{output2}" exist, skipping...')
    #         tr_reads = pl.read_parquet(tr_results_pq)
    #     else:
    #         # By LocusID load all 
    #         tr_reads_chunks = (
    #             pl.scan_csv(samples_by_donor[each_donor], separator='\t', has_header=False, schema=tr_data.complete_schema)
    #             .filter((pl.col("LocusID").is_in(lid_candidate_list)))
    #         )
    #         tr_reads = tr_reads_chunks.collect(engine="streaming")

    #         # re-save?
    #         sample_names = []
    #         for k in tr_reads.with_columns(
    #             SampleID = pl.concat_str(["Donor", "Tissue", "Core"], separator="-")
    #         )["SampleID"].to_list():
    #             sample_names.append(sample_IDs[k])

    #         # Save candidate table
    #         tr_reads = tr_reads.with_columns(pl.Series("SampleName", sample_names))
    #         tr_reads.write_parquet(tr_results_pq)
    #         tr_reads.write_csv(output2, separator="\t", include_header=True)
    #     t.end(True)
    #     ##################################################

    #     ##################################################
    #     # STEP 3
    #     # Candidate filtering
    #     t = DevTiming(f'{each_donor}, step3: candidate filter')
    #     t.start()
    #     output3 = f'{results_dir}/{each_donor}_step3-candidates_class.json'
    #     filtered_candidates = TRFilters()
    #     if os.path.isfile(output3):
    #         log.info(f'Output file "{output3}" exist, skipping...')
    #         tr_reads = pl.read_parquet(tr_results_pq)
    #     else:
    #         # Candidate filtering also by chunks
    #         loci_list = tr_reads["LocusID"].unique().to_list()
    #         n_loci_cand = len(loci_list)
    #         cand_chunk_size = max(int(n_loci_cand /100), 1)
    #         for chunk in range(0, n_loci_cand, cand_chunk_size):
    #             this_loci = loci_list[chunk:chunk+cand_chunk_size]
    #             filtered_candidates.update(filter_candidates(tr_reads.filter(pl.col("LocusID").is_in(this_loci)), config.noqc))
    #         log.info(f'Donor {each_donor} candidate results:\n{filtered_candidates}')
    #         filtered_candidates.jsonfy(output3)
    #     t.end(True)
    #     ##################################################

    #     ##################################################
    #     # STEP 4
    #     # Make results table
    #     t = DevTiming(f'{each_donor}, step4: candidate summary')
    #     t.start()
    #     output4 = f'{results_dir}/{each_donor}_step4-results_table_tissue.pq'
    #     output42 = f'{results_dir}/{each_donor}_step4-reads_per_tissue_sums.json'
    #     if os.path.isfile(output4) and os.path.isfile(output42):
    #         log.info(f'Output file "{output4}" exist, skipping...')
    #         results = pl.read_parquet(output4)
    #         with open(output42, 'r') as fin:
    #             # Parse the JSON file directly
    #             read_tissue_sum = json.load(fin, object_hook=keys_to_int)
    #         fin.close()
    #     else:
    #         filtered_candidates.from_json(output3)
    #         results, read_tissue_sum = summary_candidates(filtered_candidates, tr_results_pq, True)
    #         # Save output as pq
    #         if results is not None:
    #             results.sort(by="LocusID").write_parquet(output4)
    #         else:
    #             log.warn("No results writen")
    #         if len(read_tissue_sum) > 0:
    #             with open(output42, "w") as fout:
    #                 json.dump(read_tissue_sum, fout, indent=4)
    #             fout.close()
    #         else:
    #             log.warn("No results writen2")
    #     t.end()
    #     ##################################################

    #     ##################################################
    #     # STEP 5
    #     # catalog locus data
    #     t = DevTiming(f'{each_donor}, step5: VCF output')
    #     t.start()
    #     locus_catalog = load_adotto_bed(tr_data.adotto_v2_path)
    #     # Save output as VCF
    #     if vcf_maker(results, read_tissue_sum, locus_catalog, tr_data, each_donor, f'{results_dir}'):
    #         log.info("VCF writen")
    #     else:
    #         log.err("No results writen")
    #     t.end()
    #     output_done = f'{results_dir}/{each_donor}_steps.DONE'
    #     with open(output_done, "w") as _: pass



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--donor", required=True, type=str)        # donor name(s)/id(s)
    parser.add_argument("--qdpi-list", required=True, type=str)    # donor metadata (donors.json)
    parser.add_argument("--tr-bed", required=True, type=str)       # bed file withthe TR catalog
    parser.add_argument("--output", required=True, type=str)       # output file table
    # extra
    parser.add_argument("--selected-loci", required=False, type=str, default="")  # use only these loci
    return TRQDPIConfig(parser.parse_args())


if "__main__" == __name__:
    config = parse_args()
    assert os.path.isfile(config.qdpi_list), f'File {config.qdpi_list} does not exixist'
    assert os.path.isfile(config.tr_bed), f'File {config.tr_bed} does not exixist'
    main(config)
