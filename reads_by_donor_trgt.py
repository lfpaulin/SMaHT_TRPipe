#!/usr/bin/env python
import sys
import os
import gzip
import argparse
from pysam import AlignmentFile as BAM
from pysam import TabixFile as BED
from pysam import AlignedSegment as Read


class TRGTReadsConfig:
    def __init__(self, args: argparse.Namespace ):
        self.bams = args.bam_list
        self.donor = args.donor
        self.bed = args.bed_use
        self.trgt_flank = args.trgt_flank


def read_name_info(read_name: str|None):
    # m84 <- Revio
    rname, contig, pos, end, strand, hap, ps = "", "", "", "", "", "", ""
    if read_name is not None:
        read_name_list = read_name.split("_")
        if len(read_name_list) == 11:
            rname, contig, pos, end, _, _, _, strand, hap, ps, _ = read_name.split("_")
        elif len(read_name_list) == 12:
            rname, _, contig, pos, end, _, _, _, strand, hap, ps, _ = read_name.split("_")
        else:
            assert len(read_name_list) == 11 and len(read_name_list) == 12, read_name
    return rname, contig, pos, end, strand, hap, ps


def make_output(read_in: Read, donor_source: str, padding: int, tr_info: tuple):
    trid, _, pos, end,  = tr_info
    rname = read_in.query_name
    contig = read_in.reference_name
    hap = read_in.get_tag("HP") if read_in.has_tag("HP") else 0
    # EXPECTING 6 fileds, only using 6
    donor, tissue, core, _, _, gcc = donor_source.split("-")[:6]
    core = core if "XX" != core else f'{core}p'
    # return f'{donor}\t{tissue}\t{tissue}_{core}\t{donor_source}:{rname}\t{contig}\t{pos}\t{trid}\t{contig}:{pos}-{end}\thap{hap}\t{read_in.query_length-padding}\n'
    return f'{donor}\t{tissue}\t{core}\t{contig}\t{pos}\t{end}\t{trid}\t{hap}\t{read_in.query_length-padding}\n'


def fetch_paths(trgt_bam_list: str):
    trgt_bam_path_dict = {}
    trgt_bam_list_handler = open(trgt_bam_list)
    for line in trgt_bam_list_handler:
        # SMHT001-3AC-XX-M42-B001-uwsc.trgt.spanning.sorted.bam
        donor_path = line.rstrip("\n")
        donor = os.path.basename(donor_path.split("-")[0])
        if donor not in trgt_bam_path_dict:
            trgt_bam_path_dict[donor] = [donor_path]
        else:
            trgt_bam_path_dict[donor].append(donor_path)
    return trgt_bam_path_dict


def get_source(bam_path: str):
    # SMHT001-3AC-XX-M42-B001-uwsc.trgt.spanning.sorted.bam
    return os.path.basename(bam_path.split(".")[0])


def load_regions(use_regions: str):
    tr_regions = BED(use_regions)
    tr_regions_dict = {}
    for reg in tr_regions.fetch():
        c,p,e,i = reg.split("\t")
        trid = i.split(";")[0].split("=")[1]
        tr_regions_dict[f'{c}:{p}-{e}'] = (trid, c, p, e)
    return tr_regions_dict


def main(config: TRGTReadsConfig):
    donor = config.donor
    tr_regions = load_regions(config.bed)
    used_padding = config.trgt_flank*2
    out_prefix = "trgt"
    read: Read
    trgt_bam_donor_paths = fetch_paths(config.bams)
    assert donor in trgt_bam_donor_paths.keys(), f'Provided donor "{donor}" not in the TR file list {config.bams}\n'
    sys.stderr.write(f'{donor}: {len(trgt_bam_donor_paths[donor])}\n')
    donor_out = gzip.open(f'{out_prefix}_{donor}_readlens.tsv.gz', "wt")
    for bam_path in trgt_bam_donor_paths[donor]:
        donor_source = get_source(bam_path)
        bam_trim = BAM(bam_path)
        sys.stderr.write(f'{donor_source} {bam_trim.mapped}\n')
        counter = 0
        for reg in tr_regions.keys():
            for read in bam_trim.fetch(region=reg):
                donor_out.write(make_output(read, donor_source, used_padding, tr_regions[reg]))
                counter += 1
                if counter%1000 == 0:
                    sys.stderr.write(f'{counter}\r')
        bam_trim.close() 
        sys.stderr.write(f'\n')
    donor_out.close()


def parse_args():
    adottov2_lite_flanks = "adotto.v2.trgt.lite.flank_filt.bed.gz"
    bed_default = adottov2_lite_flanks
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam-list", required=True, type=str)     # path to list of all basm from donor
    parser.add_argument("--donor", required=True, type=str)        # donor name/ID
    parser.add_argument("--bed-use", required=False, type=str, default=bed_default)     # bed file to use
    parser.add_argument("--trgt-flank", required=False, type=int, default=50)     # size of the anchoring regions in TRGT from parameter '--flank-len 200'

    return TRGTReadsConfig(parser.parse_args())


if "__main__" == __name__:
    config = parse_args()
    assert os.path.isfile(config.bams), f'File {config.bams} does not exixist'
    assert "" != config.donor, f'Donort name not provided'
    main(config)

