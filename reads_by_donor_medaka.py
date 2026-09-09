#!/usr/bin/env python
import sys
import os
import gzip
from pysam import AlignmentFile as BAM
from pysam import TabixFile as BED
from pysam import AlignedSegment as Read


def read_name_info(read_name: str|None):
    # f48d5317-2ff4-4e8f-9cc1-50f186ee7f4b_
    # chr11_
    # 55347516_
    # 55347639_
    # pad_
    # 55347266_
    # 55347889_
    # fwd_
    # hap1_
    # phased-set54525563_
    # ploidy2
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


def make_output_old(read_in: Read, donor_source: str, padding: int, tr_info: tuple):
    trid, _, pos, end,  = tr_info
    rname, contig, pos, end, strand, hap, ps = read_name_info(read_in.query_name)
    donor, tissue, core, _, _, gcc = donor_source.split("-")
    core = core if "XX" != core else f'{core}o'
    #return f'{donor}\t{tissue}\t{tissue}_{core}\t{donor_source}:{rname}\t{contig}\t{pos}\t{end}\t{contig}:{pos}-{end}\t{hap}\t{read_in.query_length-padding}\n'
    return f'{donor}\t{tissue}\t{core}\t{contig}\t{pos}\t{end}\t{trid}\t{hap}\t{read_in.query_length-padding}\n'


def make_output(read_in: Read, donor_source: str, padding: int, tr_info: tuple):
    trid, _, pos, end,  = tr_info
    rname = read_in.query_name
    contig = read_in.reference_name
    hap = read_in.get_tag("HP") if read_in.has_tag("HP") else 0
    donor, tissue, core, _, _, gcc = donor_source.split("-")
    core = core if "XX" != core else f'{core}o'
    # return f'{donor}\t{tissue}\t{tissue}_{core}\t{donor_source}:{rname}\t{contig}\t{pos}\t{trid}\t{contig}:{pos}-{end}\thap{hap}\t{read_in.query_length-padding}\n'
    return f'{donor}\t{tissue}\t{core}\t{contig}\t{pos}\t{end}\t{trid}\t{hap}\t{read_in.query_length-padding}\n'


def fetch_paths(medaka_bam_list: str):
    medaka_bam_path_dict = {}
    medaka_bam_list_handler = open(medaka_bam_list)
    for line in medaka_bam_list_handler:
        # SMHT001-3AC-XX-M42-D003-uwsc.medaka/trimmed_reads_to_poa.bam
        line = line.rstrip("\n")
        line_list = line.split("/")
        line_len = len(line_list)
        donor_path, trim_bam = line_list[(line_len-2):]
        donor = os.path.basename(donor_path.split("-")[0])
        if donor not in medaka_bam_path_dict:
            medaka_bam_path_dict[donor] = [line]
        else:
            medaka_bam_path_dict[donor].append(line)
    return medaka_bam_path_dict


def get_source(bam_path: str):
    # SMHT001-3AC-XX-M42-D003-uwsc.medaka
    path_list = bam_path.split("/")
    path_len = len(path_list)
    donor_path, trim_bam = path_list[(path_len-2):]
    return donor_path.split(".")[0]


def load_regions(tr_bed: str):
    tr_regions = BED(tr_bed)
    tr_regions_dict = {}
    for reg in tr_regions.fetch():
        c,p,e,i = reg.split("\t")
        trid = i.split(";")[0].split("=")[1]
        tr_regions_dict[f'{c}:{p}-{e}'] = (trid, c, p, e)
    return tr_regions_dict


def main(tr_bam_list: str, donor: str, tr_bed: str):
    tr_regions = load_regions(tr_bed)
    used_padding = 250*2
    out_prefix = "medaka"
    read: Read
    medaka_bam_donor_paths = fetch_paths(tr_bam_list)
    if donor not in medaka_bam_donor_paths.keys():
        sys.stderr.write(f'Provided donor "{donor}" not in the TR file list {tr_bam_list}\n')
        sys.exit(1)
    sys.stderr.write(f'{donor}: {len(medaka_bam_donor_paths[donor])}\n')
    donor_out = gzip.open(f'{out_prefix}_{donor}_readlens.tsv.gz', "wt")
    for bam_path in medaka_bam_donor_paths[donor]:
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


if "__main__" == __name__:
    assert len(sys.argv) == 4, f'Three paramters needed in the following order: <BAM_list> <TR_bed> < donorID>'
    _, tr_bam_list, tr_bed, donor = sys.argv
    assert os.path.isfile(tr_bed), f'File {tr_bed} does not exisit'
    main(tr_bam_list, tr_bed, donor)


