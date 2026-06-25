#!/usr/bin/env python3
import json
from collections import defaultdict


def extract_parent_segments(parents, inf_hapi_data):
    the_parent_segments = defaultdict(list)
    prev_chrom = None
    most_recent_p_site = None
    p_sites = defaultdict(list)

    for marker_idx, (inf_codes, cur_chrom) in enumerate(zip(inf_hapi_data[parents]["codes"], inf_hapi_data["chr"])):
        if cur_chrom != prev_chrom:
            if prev_chrom is not None:
                assert most_recent_p_site is not None
                first_marker = most_recent_p_site - inf_hapi_data["chrstr"][prev_chrom]
                last_marker = inf_hapi_data["chrstr"][cur_chrom] - 1 - inf_hapi_data["chrstr"][prev_chrom]
                if last_marker - first_marker >= 100:
                    the_parent_segments[prev_chrom].append([first_marker, last_marker])
            most_recent_p_site = None
            prev_chrom = cur_chrom

        if inf_codes is not None and (inf_codes[0] == 'P' or inf_codes[0] == 'PC' or inf_codes[0] == 'PA'):
            p_sites[str(cur_chrom)].append(marker_idx)
            if most_recent_p_site is not None:
                first_marker = most_recent_p_site - inf_hapi_data["chrstr"][cur_chrom]
                last_marker = marker_idx - 1 - inf_hapi_data["chrstr"][cur_chrom]
                if last_marker - first_marker >= 100:
                    the_parent_segments[prev_chrom].append([first_marker, last_marker])
            most_recent_p_site = marker_idx

    assert most_recent_p_site is not None
    first_marker = most_recent_p_site - inf_hapi_data["chrstr"][prev_chrom]
    last_marker = marker_idx - 1 - inf_hapi_data["chrstr"][prev_chrom]
    if last_marker - first_marker >= 100:
        the_parent_segments[prev_chrom].append([first_marker, last_marker])

    return the_parent_segments, p_sites


path = r'C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/hapi2_trio/hapi2_out/all.json'
with open(path, 'r', encoding='utf-8') as handle:
    data = json.load(handle)

parents = 'FAM1:100001-FAM1:100002'
segments, p_sites = extract_parent_segments(parents, data)
chrom = '22'
chr_start = data['chrstr'][chrom]
chr_end = None
chrom_names = list(data['chrstr'].keys())
chrom_order = sorted(chrom_names, key=lambda x: data['chrstr'][x])
idx = chrom_order.index(chrom)
if idx + 1 < len(chrom_order):
    chr_end = data['chrstr'][chrom_order[idx + 1]] - 1
else:
    chr_end = len(data['chr']) - 1

print(f'Chr{chrom} global start index: {chr_start}')
print(f'Chr{chrom} global end index:   {chr_end}')
print(f'Chr{chrom} marker count:       {chr_end - chr_start + 1}')
print()
print('P/PC/PA sites on chr22:')
for global_idx in p_sites[chrom]:
    rel_idx = global_idx - chr_start
    bp = data['physpos'][global_idx]
    code = data[parents]['codes'][global_idx]
    print(f'  global={global_idx} rel={rel_idx} bp={bp} code={code}')
print()
print('Extracted HAPI parent segments on chr22 (>=100 markers only):')
for seg_idx, (start_rel, end_rel) in enumerate(segments[chrom], start=1):
    global_start = chr_start + start_rel
    global_end = chr_start + end_rel
    start_bp = data['physpos'][global_start]
    end_bp = data['physpos'][global_end]
    print(
        f'  seg{seg_idx}: rel={start_rel}-{end_rel} '
        f'global={global_start}-{global_end} bp={start_bp}-{end_bp} '
        f'markers={end_rel - start_rel + 1}'
    )

print()
print('Gaps between extracted segments:')
prev_end = None
for seg_idx, (start_rel, end_rel) in enumerate(segments[chrom], start=1):
    if prev_end is not None:
        gap_start = prev_end + 1
        gap_end = start_rel - 1
        if gap_end >= gap_start:
            gs = chr_start + gap_start
            ge = chr_start + gap_end
            print(
                f'  gap before seg{seg_idx}: rel={gap_start}-{gap_end} '
                f'bp={data["physpos"][gs]}-{data["physpos"][ge]} markers={gap_end-gap_start+1}'
            )
    prev_end = end_rel

if segments[chrom]:
    first_start = segments[chrom][0][0]
    if first_start > 0:
        gs = chr_start
        ge = chr_start + first_start - 1
        print(
            f'  gap before seg1: rel=0-{first_start-1} '
            f'bp={data["physpos"][gs]}-{data["physpos"][ge]} markers={first_start}'
        )
    last_end = segments[chrom][-1][1]
    chr_last_rel = chr_end - chr_start
    if last_end < chr_last_rel:
        gs = chr_start + last_end + 1
        ge = chr_end
        print(
            f'  tail after last seg: rel={last_end+1}-{chr_last_rel} '
            f'bp={data["physpos"][gs]}-{data["physpos"][ge]} markers={chr_last_rel-last_end}'
        )
