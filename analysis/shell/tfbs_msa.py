#!/usr/bin/env python3
"""Build full letter-level alignments for SuppFig_E_MSA.

For element E (peak_108181027) and up to three other peaks that are only
partially conserved outside Bdor/Bcor (shortest first), extract the chained
homologue of peak +/-50 bp from each species and align all six with
MAFFT L-INS-i. Uses the chained chain endpoints recorded in
fulllength_identity.tsv.

Writes into <out>/msa/: <peak>.fa, <peak>.mafft.fa, meta.tsv
Usage: tfbs_msa.py <out-dir>
"""
import sys, os, csv, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

OUT = sys.argv[1]
MDIR = f'{OUT}/msa'
os.makedirs(MDIR, exist_ok=True)
HIC = '/path/to/data/Hic'
MAFFT = '/path/to/conda/envs/genera/bin/mafft'
NP = f'{HIC}/Bdor/results_hicexplorer/Bd_MO.pval0.01.300K.narrowPeak'
FLANK = 50

# ---------- inputs ----------
sp_locus = {}
for l in open(f'{OUT}/species_locus.tsv'):
    v = l.rstrip('\n').split('\t')
    if not v or not v[0]:
        continue
    sp_locus[v[0]] = dict(fa=v[1], chr=v[2], lo=int(v[3]), hi=int(v[4]),
                          strand=v[5], n1=int(v[6]))
OTH = [s for s in sp_locus if s != 'Bdor']

queries = {}          # name -> dict with core coords
for l in open(f'{OUT}/queries.bed'):
    v = l.rstrip('\n').split('\t')
    queries[v[3]] = dict(s=int(v[1]), e=int(v[2]))

# MACS2 summit (best -log10q row per interval, same pick as queries.bed)
summit = {}
for l in open(NP):
    v = l.split('\t')
    if v[0] != 'NC_064303.1':
        continue
    key = (int(v[1]), int(v[2]))
    if key not in summit or float(v[8]) > summit[key][0]:
        summit[key] = (float(v[8]), int(v[9]))

fl = {}
for r in csv.DictReader(open(f'{OUT}/fulllength_identity.tsv'), delimiter='\t'):
    fl[(r['query'], r['species'])] = r
hom = {}
for r in csv.DictReader(open(f'{OUT}/homolog_table.tsv'), delimiter='\t'):
    hom[(r['query'], r['species'])] = r['homolog'] == 'Y'

peaks = [q for q in queries if q.startswith('peak_')]
# pick: E plus the partially-conserved peaks (homologs outside Bcor but not all
# five species >=80% full-length), shortest first, at most 3
def pick():
    partial = []
    for q in peaks:
        if q == 'peak_108181027':
            continue
        hs = [s for s in OTH if hom.get((q, s))]
        if len(hs) < 2:                     # nothing beyond Bcor
            continue
        vals = [float(fl[(q, s)]['fl_ident'] or 0) for s in hs]
        if all(v >= 80 for v in vals):
            continue                        # fully conserved -> not "partial"
        partial.append(q)
    partial.sort(key=lambda q: queries[q]['e'] - queries[q]['s'])
    return ['peak_108181027'] + partial[:3]

chosen = pick()
print('msa peaks:', chosen)

# ---------- load sequences ----------
bdor = next(SeqIO.parse(f'{OUT}/locus_Bdor.fa', 'fasta'))
BLO = int(bdor.id.split(':')[-1].split('-')[0])
BSEQ = bdor.seq.upper()
locus_seq = {sp: next(SeqIO.parse(f'{OUT}/locus_{sp}.fa', 'fasta')) for sp in OTH}

meta = open(f'{MDIR}/meta.tsv', 'w')
meta.write('query\tcore_s\tcore_e\tsummit\tflank_s\tflank_e\tspecies\n')
for q in chosen:
    cs, ce = queries[q]['s'], queries[q]['e']
    fs, fe = cs - FLANK, ce + FLANK
    sm = summit.get((cs, ce), (0, -1))[1]
    smpos = cs + sm if sm >= 0 else ''
    seqs = {'Bdor': str(BSEQ[fs - BLO:fe - BLO])}
    qlen = (ce + 100) - (cs - 100)          # the +-100 query length
    tlo, thi = (fs - (cs - 100)) + 1, fe - (cs - 100)
    for sp in OTH:
        r = fl.get((q, sp))
        if not r or not r['fl_ident'] or not hom.get((q, sp)):
            continue
        slo = int(r['s_start']); shi = int(r['s_end'])
        qlo, qhi = int(r['q_lo']), int(r['q_hi'])
        if r['strand'] == '+':
            lo = slo - (qlo - tlo); hi = shi + (thi - qhi)
        else:
            lo = slo - (thi - qhi); hi = shi + (qlo - tlo)
        sl = sp_locus[sp]
        lo = max(lo, sl['lo']); hi = min(hi, sl['hi'])
        i0, i1 = lo - sl['lo'], hi - sl['lo']
        s = locus_seq[sp].seq[i0:i1]
        if r['strand'] == '-':
            s = s.reverse_complement()
        seqs[sp] = str(s).upper()
    with open(f'{MDIR}/{q}.fa', 'w') as f:
        for sp in ['Bdor'] + OTH:
            if sp in seqs:
                f.write(f'>{sp}\n{seqs[sp]}\n')
    subprocess.run([MAFFT, '--localpair', '--maxiterate', '1000',
                    f'{MDIR}/{q}.fa'],
                   stdout=open(f'{MDIR}/{q}.mafft.fa', 'w'),
                   stderr=open(f'{MDIR}/{q}.mafft.log', 'w'), check=True)
    meta.write(f'{q}\t{cs}\t{ce}\t{smpos}\t{fs}\t{fe}\t'
               f'{",".join(s for s in ["Bdor"] + OTH if s in seqs)}\n')
meta.close()
print('msa written to', MDIR)
