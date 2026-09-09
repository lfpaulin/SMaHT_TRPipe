#!/usr/bin/env python
import os
import sys
import gzip
import json


bias_decoder2 ={
    91: "No bias assessed, only three clusters",
    "91": "No bias assessed, only three clusters",
    92: "No bias assessed, no replicates available in one tech",
    "92": "No bias assessed, no replicates available in one tech",
    0: "No bias detected, both techs have the same number of clusters",
    "0": "No bias detected, both techs have the same number of clusters",
    1: "Some detected, both techs have 1+ additional clusters",
    "1": "Some detected, both techs have 1+ additional clusters",
    2: "Bias detected in ONT, it has more clusters",
    "2": "Bias detected in ONT, it has more clusters",
    3: "Bias detected in PB, it has more clusters",
    "3": "Bias detected in PB, it has more clusters"
}

bias_decoder ={
    91: "NoClueB",
    "91": "NoClueB",
    92: "NoClueA",
    "92": "NoClueA",
    0: "OK",
    "0": "OK",
    1: "NoClueC",
    "1": "NoClueC",
    2: "BAD:ONT",
    "2": "BAD:ONT",
    3: "BAD:PB",
    "3": "BAD:PB"
}


def decode_bias(biases: dict) -> str:
    bias_str = ",".join([f'{bias_decoder[k]}:{v}' for k,v in biases.items()])
    return bias_str


def load_bed(bed_plus: list):
    bed, out_bed = {}, False
    if len(bed_plus) > 0:
        bed_file = bed_plus.pop(0)
        bed_handler = open(bed_file, "r") if "gz" not in bed_file else gzip.open(bed_file, "rt")
        out_bed = True
        for line in bed_handler:
            chro, pos, end, tr_info = line.rstrip('\n').split('\t')
            lid_bed = [i.split("=").pop() for i in tr_info.split(";") if "ID=" in i].pop()
            bed[lid_bed] = (chro, pos, end, tr_info)

    return bed, out_bed


def parse_smears(outfile: str):
    if os.path.isfile(outfile):
        with open(outfile, 'r') as fin:
            # Parse the JSON file directly
            infos = json.load(fin)    
    else:
        donor: str = ""
        infos = {}
        for line in sys.stdin:
            if line.startswith("#"):
                donor = line.rstrip("\n").split(" ")[1]
            else:
                assert "" != donor, f'ERROR: Donor cannot be empty {donor} (?)'
                lid, bias, *_ = line.rstrip("\n").split(" ")
                assert lid != "", f'ERROR: entry has locusID: {line}'
                if lid not in infos:
                    infos[lid] = [(donor, bias)]
                else:
                    infos[lid].append((donor, bias))

        with open(outfile, "w") as fout:
            json.dump(infos, fout, indent=4)
        fout.close()

    return infos


def main(outfile: str, bed_plus: list):
    infos = parse_smears(outfile)
    bed, out_bed = load_bed(bed_plus)  #{}, False

    lid_done = set()
    chro, pos, end, tr_info = "", "", "", ""
    for lid in infos.keys():
        bias_data = {}
        bias_data_ext = {}
        in_n_donor = 0
        nbias = 0
        for donor, bias in infos[lid]:
            in_n_donor += 1
            if bias not in bias_data:
                nbias += 1
                bias_data[bias] = 1
                bias_data_ext[bias] = [donor]
            else:
                bias_data[bias] += 1
                bias_data_ext[bias].append(donor)

        if out_bed:
            chro, pos, end, tr_info = bed[lid]

        decoded_bias = decode_bias(bias_data)

        if nbias == 1:
            bias = list(bias_data.keys()).pop()
            bias_out = bias_decoder[bias] if "OK" not in decoded_bias else "OK"
            if out_bed:
                print(f'{chro}\t{pos}\t{end}\t{tr_info}\t{bias_out}\t{decoded_bias}')
            else:
                print(f'Locus {lid} results: {bias_out}|{bias_decoder2[bias]} | {decoded_bias}')
        else:
            bias, bias_count = sorted([(k,v) for k,v in bias_data.items()], key=lambda bias: bias[1], reverse=True).pop(0)
            if "OK" in decoded_bias:
                ok_donor = ",".join(bias_data_ext["0"])
                bias_out = "OK"
                decoded_bias = f'{decoded_bias}|{ok_donor}'
            else:
                bias_out = bias_decoder[bias]
            if bias_count/float(in_n_donor) >= 0.5:
                if out_bed:
                    print(f'{chro}\t{pos}\t{end}\t{tr_info}\t{bias_out}\t{decoded_bias}')
                else:
                    print(f'Locus {lid} results: {bias_out}|{bias_decoder2[bias]} | {decoded_bias}')
            else:
                if out_bed:
                    if "OK" in decoded_bias:
                        ok_donor = ",".join(bias_data_ext["0"])
                        bias_out = "OK"
                        decoded_bias = f'{decoded_bias}|{ok_donor}'
                    else:
                        bias_out = "NoClue0"
                    print(f'{chro}\t{pos}\t{end}\t{tr_info}\t{bias_out}\t{decoded_bias}')

        # for bias in bias_data:
        #     bias_prop = bias_data[bias]/float(in_n_donor)
        #     if bias_prop > 0.5 and nbias > 1:
        #         lid_done.add(lid)
        #         if out_bed:
        #             chro, pos, end, tr_info = bed[lid]
        #             print(f'{chro}\t{pos}\t{end}\t{tr_info}\t{bias_decoder[bias]}\t{decode_bias(bias_data)}')
        #         else:
        #             print(f'Locus {lid} results: {bias_decoder2[bias]} | {decode_bias(bias_data)}')
        #     elif bias_prop > 0.5 and nbias == 1:
        #         lid_done.add(lid)
        #         if out_bed:
        #             chro, pos, end, tr_info = bed[lid]
        #             print(f'{chro}\t{pos}\t{end}\t{tr_info}\t{bias_decoder[bias]}\t{decode_bias(bias_data)}')
        #         else:
        #             print(f'Locus {lid} results: {bias_decoder2[bias]} | {decode_bias(bias_data)}')
        #     else:
        #         pass
        # if lid not in lid_done:
        #     lid_done.add(lid)
        #     if out_bed:
        #         chro, pos, end, tr_info = bed[lid]
        #         print(f'{chro}\t{pos}\t{end}\t{tr_info}\tNoClue0\t{decode_bias(bias_data)}')


if "__main__" == __name__:
    assert len(sys.argv) in {2, 3}, f'ERROR: output file is needed, input is from STDIN, additional to bed file. Given: {sys.argv}'
    _, out, *bed = sys.argv
    main(out, bed)
