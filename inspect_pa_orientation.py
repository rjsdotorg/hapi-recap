#!/usr/bin/env python3
import json

p = r'C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/hapi2_trio/hapi2_out/all.json'
d = json.load(open(p))
key = 'FAM1:100001-FAM1:100002'
chr_start = d['chrstr']['22']
h0 = d[key]['parhaps'][0]
h1 = d[key]['parhaps'][1]

print("Codes and parhap orientation at PA boundaries (rel 6120-6240)")
print(f"{'rel':>6} {'bp':>12} {'code':<14} {'h0[0]':>5} {'h0[1]':>5} {'h1[0]':>5} {'h1[1]':>5}")
for rel in range(6120, 6240):
    gidx = chr_start + rel
    code = d[key]['codes'][gidx]
    if code is None:
        continue
    print(f"{rel:>6} {d['physpos'][gidx]:>12} {str(code):<14} "
          f"{str(h0[0][gidx]):>5} {str(h0[1][gidx]):>5} "
          f"{str(h1[0][gidx]):>5} {str(h1[1][gidx]):>5}")

print()
print("Checking whether haplotype assignment actually SWITCHES at each PA boundary:")
prev_h0_0 = None
for rel in range(6120, 6240):
    gidx = chr_start + rel
    code = d[key]['codes'][gidx]
    if code is None or code[0] not in ('P', 'PC', 'PA'):
        continue
    cur = (h0[0][gidx], h0[1][gidx], h1[0][gidx], h1[1][gidx])
    switch = "SWITCH" if (prev_h0_0 is not None and cur[0] != prev_h0_0) else "same"
    prev_h0_0 = cur[0]
    print(f"rel={rel} bp={d['physpos'][gidx]} code={code[0]} h0[0]={cur[0]} {switch}")
