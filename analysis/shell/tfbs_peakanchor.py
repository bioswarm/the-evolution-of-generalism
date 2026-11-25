#!/usr/bin/env python3
"""Peak-anchored homolog + TFBS analysis for the Gr28b locus (round 2).
Reads:  <out>/species_locus.tsv, <out>/queries_flank100.bed, <out>/queries_flank100.fa,
        <out>/blast/locus_<sp>.tsv (blastn outfmt 6)
Writes into <out>: homolog_table.tsv, all_hsps.tsv, E_homologs.fa/.mafft.fa,
        control_homologs.fa/.mafft.fa, identities_peakanchored.tsv,
        Gr28b_E_TFBS_peakanchored.tsv, Gr28b_E_TFBS_counts_peakanchored.tsv,
        sites_on_aln_peakanchored.tsv, bmin_atac_homologs.tsv, genome_hits.tsv
"""
import sys, os, re, csv, subprocess, glob

OUT = sys.argv[1]
TDIR = os.path.dirname(OUT)                      # results/locus/tfbs
GS = '/path/to/data/genesets'
HIC = '/path/to/data/Hic'
BLAST = '/opt/ncbi-blast-2.10.0+/bin'
MAFFT = '/path/to/conda/envs/genera/bin/mafft'
MEME = '/path/to/conda/envs/anvio-9/bin'
CHIP = '/path/to/conda/envs/chipseq/bin/python'
JASPAR = f'{TDIR}/jaspar_insects.meme'

# ---------- inputs ----------
sp_locus = {}
for l in open(f'{OUT}/species_locus.tsv'):
    v = l.rstrip('\n').split('\t')
    if not v or not v[0]: continue
    sp_locus[v[0]] = dict(fa=v[1], chr=v[2], lo=int(v[3]), hi=int(v[4]),
                          strand=v[5], n1=int(v[6]))
queries = []
for l in open(f'{OUT}/queries_flank100.bed'):
    v = l.rstrip('\n').split('\t')
    queries.append(dict(chrom=v[0], s=int(v[1]), e=int(v[2]), name=v[3],
                        peak_id=v[4], mlog10q=v[5], core_s=int(v[6]), core_e=int(v[7])))
QLEN = {q['name']: q['e'] - q['s'] for q in queries}
OTH = [s for s in sp_locus if s != 'Bdor']

# chained-HSP query intervals from fulllength_identity.tsv (produced by
# tfbs_fulllength.py, which runs first): used to score 'no_homology' for sites
# and homologous-core identities -- a Bdor position outside a species' chain
# has no detectable similarity.
qints = {}
if os.path.exists(f'{OUT}/fulllength_identity.tsv'):
    with open(f'{OUT}/fulllength_identity.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            ints = []
            if r.get('chain_qints'):
                for p in r['chain_qints'].split(';'):
                    a_, b_ = p.split('-')
                    ints.append((int(a_), int(b_)))
            qints[(r['query'], r['species'])] = ints
TAG2Q = {'E': 'peak_108181027', 'control': 'control_108180599'}

# ---------- parse locus HSPs ----------
def read_hsps(f):
    out = []
    if not os.path.exists(f): return out
    for l in open(f):
        v = l.rstrip('\n').split('\t')
        if len(v) < 12: continue
        out.append(dict(q=v[0].split('::')[0], s=v[1], pident=float(v[2]), alen=int(v[3]),
                        qlen=int(v[4]), slen=int(v[5]), qs=int(v[6]), qe=int(v[7]),
                        ss=int(v[8]), se=int(v[9]), ev=float(v[10]), bs=float(v[11])))
    return out

hsps = {sp: read_hsps(f'{OUT}/blast/locus_{sp}.tsv') for sp in OTH}
with open(f'{OUT}/all_hsps.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['species', 'query', 'sseqid', 'pident', 'alen', 'qcov', 'qs', 'qe',
                'ss', 'se', 'evalue', 'bitscore'])
    for sp in OTH:
        for h in hsps[sp]:
            w.writerow([sp, h['q'], h['s'], h['pident'], h['alen'],
                        round(h['alen'] / h['qlen'], 3), h['qs'], h['qe'],
                        h['ss'], h['se'], h['ev'], h['bs']])

# ---------- synteny cross-check ----------
syn_beds = {'Bcor': 'Bdor_vs_Bcor', 'Bmin': 'Bdor_vs_Bmin', 'Zcuc': 'Bdor_vs_Bcuc',
            'Ztau': 'Bdor_vs_Btau', 'Btsu': None}
SYN = {}
for sp, b in syn_beds.items():
    SYN[sp] = []
    if b:
        for i, l in enumerate(open(f'{HIC}/Bdor/results_hicexplorer/{b}.gui-gapped1000.genomic.bed')):
            if i == 0 or not l.strip(): continue
            v = l.split('\t')
            if v[0] == 'NC_064303.1':
                SYN[sp].append((int(v[1]), int(v[2])))
def in_syn(sp, s, e):
    return 'Y' if any(a < e and b > s for a, b in SYN[sp]) else ('NA' if syn_beds.get(sp) is None else 'N')

# ---------- best HSP + collinearity ----------
# order rank: queries sorted by Bdor midpoint (E = element peak)
q_by_mid = sorted(queries, key=lambda q: (q['core_s'] + q['core_e']) / 2)
q_rank = {q['name']: i + 1 for i, q in enumerate(q_by_mid)}

homolog = {}     # (query, sp) -> dict
for sp in OTH:
    sl = sp_locus[sp]
    per_q = {}
    for h in hsps[sp]:
        if h['ev'] >= 1e-3: continue
        per_q.setdefault(h['q'], []).append(h)
    best = {q: max(hs, key=lambda h: (h['bs'], -h['ev'])) for q, hs in per_q.items()}
    # order: sort the best hits by gene-relative midpoint (transcription orientation)
    rel = {}
    for q, h in best.items():
        mid = (h['ss'] + h['se']) / 2 + sl['lo'] - 1
        rel[q] = (mid - sl['n1']) if sl['strand'] == '+' else (sl['n1'] - mid)
    hit_rank = {q: i + 1 for i, q in enumerate(sorted(rel, key=rel.get))}
    for q in [q['name'] for q in queries]:
        h = best.get(q)
        ent = dict(query=q, species=sp, homolog='N', sseqid='', s_start='', s_end='',
                   strand='', evalue='', pident='', qcov='', rel_mid_kb='',
                   q_rank=q_rank[q], hit_rank='', strand_ok='', order_ok='',
                   synteny=in_syn(sp, *[[x['core_s'], x['core_e']] for x in queries if x['name'] == q][0]))
        if h:
            strand_ok = (h['ss'] < h['se']) == (sl['strand'] == '+')
            ent.update(sseqid=h['s'], evalue=f'{h["ev"]:.2g}', pident=f'{h["pident"]:.1f}',
                       qcov=f'{h["alen"]/h["qlen"]:.2f}',
                       rel_mid_kb=f'{rel[q]/1000:.3f}', hit_rank=hit_rank[q],
                       strand_ok='Y' if strand_ok else 'N')
            lo_c = min(h['ss'], h['se']) + sl['lo'] - 1
            hi_c = max(h['ss'], h['se']) + sl['lo'] - 1
            ent.update(s_start=lo_c, s_end=hi_c,
                       strand='+' if h['ss'] < h['se'] else '-')
            # order consistent if rank among this species' hits == rank among queries (within +/-1 of tied spread)
            ent['order_ok'] = 'Y' if hit_rank[q] == q_rank[q] else 'N'
            if strand_ok:
                ent['homolog'] = 'Y'
                ent['hsp'] = h
        homolog[(q, sp)] = ent

# ---------- genome-wide fallback for queries without a locus homolog ----------
from Bio import SeqIO
from Bio.Seq import Seq
qseq = {}
for r in SeqIO.parse(f'{OUT}/queries_flank100.fa', 'fasta'):
    qseq[r.id.split('::')[0].split()[0]] = str(r.seq)
locus_seq = {sp: next(SeqIO.parse(f'{OUT}/locus_{sp}.fa', 'fasta')).seq for sp in OTH}
gw = []
for q in queries:
    for sp in OTH:
        if homolog[(q['name'], sp)]['homolog'] == 'Y': continue
        h = homolog[(q['name'], sp)]
        if not h['evalue']:          # no HSP at all in the locus
            qfa = f'{OUT}/_q_{q["name"]}.fa'
            with open(qfa, 'w') as f:
                f.write('>' + q['name'] + '\n' + qseq[q['name']] + '\n')
            r = subprocess.run([f'{BLAST}/blastn', '-query', qfa,
                                '-subject', sp_locus[sp]['fa'],
                                '-task', 'blastn', '-word_size', '7', '-evalue', '1e-2',
                                '-dust', 'no', '-soft_masking', 'false',
                                '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore'],
                               capture_output=True, text=True)
            hits = [l_.split('\t') for l_ in r.stdout.splitlines() if l_.strip()]
            hits.sort(key=lambda v: float(v[10]))
            if hits:
                v = hits[0]
                gw.append([q['name'], sp, v[1], v[2], v[3], v[8], v[9], v[10]])
            os.remove(qfa)
with open(f'{OUT}/genome_hits.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['query', 'species', 'sseqid', 'pident', 'alen', 'sstart', 'send', 'evalue'])
    w.writerows(gw)

# ---------- extract E/control homologs (HSP extended to full query span) ----------


def extract(sp, h):
    """subject subsequence matching the full query span; '+'-oriented to Bdor."""
    sl = sp_locus[sp]
    seq = locus_seq[sp]
    if h['ss'] < h['se']:
        lo = h['ss'] - (h['qs'] - 1); hi = h['se'] + (h['qlen'] - h['qe'])
        lo = max(lo, 1); hi = min(hi, len(seq))
        return str(seq[lo - 1:hi]), lo + sl['lo'] - 1, hi + sl['lo'] - 1, '+'
    else:
        lo = h['se'] - (h['qlen'] - h['qe']); hi = h['ss'] + (h['qs'] - 1)
        lo = max(lo, 1); hi = min(hi, len(seq))
        return str(Seq(str(seq[lo - 1:hi])).reverse_complement()), lo + sl['lo'] - 1, hi + sl['lo'] - 1, '-'

def ext_span(sp, h):
    sl = sp_locus[sp]; L = len(locus_seq[sp])
    if h['ss'] < h['se']:
        lo = max(h['ss'] - (h['qs'] - 1), 1); hi = min(h['se'] + (h['qlen'] - h['qe']), L)
    else:
        lo = max(h['se'] - (h['qlen'] - h['qe']), 1); hi = min(h['ss'] + (h['qs'] - 1), L)
    return lo + sl['lo'] - 1, hi + sl['lo'] - 1
for q in queries:
    for sp in OTH:
        e = homolog[(q['name'], sp)]
        if e['homolog'] == 'Y':
            e['ext_start'], e['ext_end'] = ext_span(sp, e['hsp'])

for tag in ('peak_108181027', 'control_108180599'):
    fa_out = open(f'{OUT}/{"E" if tag.startswith("peak") else "control"}_homologs.fa', 'w')
    q = [x for x in queries if x['name'] == tag][0]
    fa_out.write(f'>Bdor\n{qseq[tag]}\n')
    for sp in OTH:
        ent = homolog[(tag, sp)]
        if ent['homolog'] == 'Y':
            s, lo, hi, st = extract(sp, ent['hsp'])
            ent['ext_start'], ent['ext_end'] = lo, hi
            fa_out.write(f'>{sp}\n{s}\n')
    fa_out.close()
    subprocess.run([MAFFT, '--localpair', '--maxiterate', '1000', '--thread', '4',
                    f'{OUT}/{"E" if tag.startswith("peak") else "control"}_homologs.fa'],
                   stdout=open(f'{OUT}/{"E" if tag.startswith("peak") else "control"}_homologs.mafft.fa', 'w'),
                   stderr=open(f'{OUT}/mafft_{tag}.log', 'w'), check=True)

# ---------- identities over the E alignment ----------
def load_aln(f):
    a, n = {}, None
    for l in open(f):
        if l.startswith('>'): n = l[1:].split()[0]; a[n] = ''
        else: a[n] += l.strip()
    return a
for tag, aln_f, core in (('E', 'E_homologs.mafft.fa', (101, 528)),
                         ('control', 'control_homologs.mafft.fa', (101, 528))):
    if not os.path.exists(f'{OUT}/{aln_f}'): continue
    a = load_aln(f'{OUT}/{aln_f}')
    bd = a['Bdor']; p2c, k = {}, 0
    for c, ch in enumerate(bd):
        if ch != '-': k += 1; p2c[k] = c
    cols = [p2c[p] for p in range(core[0], core[1] + 1)]
    qn = TAG2Q[tag]
    with open(f'{OUT}/identities_{tag}.tsv', 'w') as f:
        f.write('species\tid_core\tncols_both\tid_homologous\tncols_homologous\tstatus_in_aln\n')
        for sp in a:
            s = a[sp]
            both = [c for c in cols if s[c] != '-' and bd[c] != '-']
            m = sum(1 for c in both if s[c] == bd[c])
            hcols = both
            if sp != 'Bdor':
                ints = qints.get((qn, sp), [])
                # restrict to alignment columns whose Bdor query position is
                # inside the species' chained HSPs (no force-aligned positions)
                p_of = {c: pp for pp, c in p2c.items()}
                hcols = [c for c in both
                         if any(x <= p_of[c] <= y for x, y in ints)]
            m2 = sum(1 for c in hcols if s[c] == bd[c])
            f.write(f'{sp}\t{m/len(cols):.4f}\t{len(both)}\t'
                    f'{m2/len(hcols):.4f}\t{len(hcols)}\taligned\n')
    with open(f'{OUT}/identities_{tag}.tsv', 'a') as f:
        for sp in OTH:
            if sp not in a:
                f.write(f'{sp}\tNA\t0\tNA\t0\tno_detectable_homolog\n')

# ---------- Bmin ATAC at each homolog ----------
BW = glob.glob('/public_data/ATAC-seq/deprecate/Bmin/run_Bmin/atac/*/call-macs2_signal_track/shard-0/execution/Bm_MO.*fc.signal.bigwig')[0]
bmin_peaks = [(int(v[1]), int(v[2])) for l in open(f'{HIC}/Bmin/results_hicexplorer/Bm_MO.pval0.01.300K.narrowPeak')
              for v in [l.split('\t')] if v[0] == 'BM1']
sl = sp_locus['Bmin']
# 428-bp non-overlapping tiles over the Bmin locus +/-50 kb (chrom BmChr1)
regions, meta = [], []
gs, ge = sl['lo'], sl['hi']
for s in range(gs, ge - 428, 428):
    regions.append(f'BmChr1:{s}-{s + 428}'); meta.append(('tile', s))
for q in queries:
    e = homolog[(q['name'], 'Bmin')]
    if e['homolog'] == 'Y':
        regions.append(f'BmChr1:{e["ext_start"]}-{e["ext_end"]}'); meta.append(('q', q['name']))
r = subprocess.run([CHIP, 'shell/tfbs_bwmean.py', BW] + regions, capture_output=True, text=True)
bwv = {l.split('\t')[0]: float(l.split('\t')[1]) for l in r.stdout.splitlines() if '\t' in l}
tile_vals = sorted(bwv[regions[i]] for i, m in enumerate(meta) if m[0] == 'tile')
def bwmean_reg(lo, hi):
    return bwv.get(f'BmChr1:{lo}-{hi}', float('nan'))

with open(f'{OUT}/bmin_atac_homologs.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['query', 'bmin_homolog', 'ext_start', 'ext_end', 'mean_fc',
                'peak_overlap', 'fc_percentile_vs_428bp_tiles'])
    for q in queries:
        ent = homolog[(q['name'], 'Bmin')]
        if ent['homolog'] != 'Y':
            w.writerow([q['name'], 'N', '', '', '', '', ''])
            continue
        lo, hi = ent['ext_start'], ent['ext_end']
        fc = bwmean_reg(lo, hi)
        pk = 'Y' if any(a < hi and b > lo for a, b in bmin_peaks) else 'N'
        pct = 100 * sum(1 for t in tile_vals if t <= fc) / len(tile_vals)
        w.writerow([q['name'], 'Y', lo, hi, f'{fc:.4f}', pk, f'{pct:.1f}'])

# ---------- homolog table ----------
with open(f'{OUT}/homolog_table.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['query', 'bd_start', 'bd_end', 'species', 'homolog', 'sseqid',
                's_start', 's_end', 'strand', 'evalue', 'pident', 'qcov',
                'rel_mid_kb_to_N1exon', 'query_rank', 'hit_rank', 'strand_ok',
                'order_ok', 'bd_in_synteny_seg'])
    for q in q_by_mid:
        for sp in OTH:
            e = homolog[(q['name'], sp)]
            w.writerow([q['name'], q['core_s'], q['core_e'], sp, e['homolog'],
                        e['sseqid'], e['s_start'], e['s_end'], e['strand'],
                        e['evalue'], e['pident'], e['qcov'], e['rel_mid_kb'],
                        e['q_rank'], e['hit_rank'], e['strand_ok'], e['order_ok'],
                        e['synteny']])

# ---------- FIMO on E/control homolog sequences ----------
for tag in ('E', 'control'):
    fa = f'{OUT}/{tag}_homologs.fa'
    if not os.path.exists(fa): continue
    for r in SeqIO.parse(fa, 'fasta'):
        sp = r.id
        with open(f'{OUT}/_tmp_{sp}.fa', 'w') as f:
            f.write(f'>{sp}\n{r.seq}\n')
        for T in ('1e-4', '1e-3'):
            d = f'{OUT}/fimo/{tag}_{sp}_{T}'
            os.makedirs(d, exist_ok=True)
            subprocess.run([f'{MEME}/fimo', '--oc', d, '--thresh', T,
                            '--bgfile', f'{TDIR}/bg_{sp}.txt', JASPAR,
                            f'{OUT}/_tmp_{sp}.fa'],
                           capture_output=True)
        os.remove(f'{OUT}/_tmp_{sp}.fa')

# ---------- classification on the new E alignment ----------
motif_name = {}
for m in re.finditer(r'MOTIF (\S+)(?: +(\S+))?', open(JASPAR).read()):
    motif_name[m.group(1)] = m.group(2) or m.group(1)

def fimo_hits(tag, sp, thr):
    f = f'{OUT}/fimo/{tag}_{sp}_{thr}/fimo.tsv'
    if not os.path.exists(f): return []
    out = []
    for l in open(f):
        if l.startswith('#') or l.startswith('motif_id'): continue
        v = l.rstrip('\n').split('\t')
        if len(v) < 9: continue
        out.append(dict(motif=v[0], tf=motif_name.get(v[0], v[1]), a=int(v[3]),
                        b=int(v[4]), strand=v[5], p=float(v[7]),
                        q=float(v[8]) if v[8] not in ('', '-') else 1.0))
    return out

def classify_region(tag):
    aln = load_aln(f'{OUT}/{tag}_homologs.mafft.fa')
    bd = aln['Bdor']; p2c, k = {}, 0
    for c, ch in enumerate(bd):
        if ch != '-': k += 1; p2c[k] = c
    pos2col = {'Bdor': p2c}
    for sp, s in aln.items():
        if sp == 'Bdor': continue
        p, kk = {}, 0
        for c, ch in enumerate(s):
            if ch != '-': kk += 1; p[kk] = c
        pos2col[sp] = p
    core_cols = {p2c[p] for p in range(101, 529)}          # 428-bp core of the 628-bp query
    H1 = {sp: fimo_hits(tag, sp, '1e-4') for sp in aln}
    H3 = {sp: fimo_hits(tag, sp, '1e-3') for sp in aln}
    qn = TAG2Q[tag]
    def chain_cov(sp, a_, b_):
        """fraction of Bdor segment positions a_..b_ inside sp's chained HSPs."""
        ints = qints.get((qn, sp), [])
        if not ints:
            return 0.0
        c = sum(min(b, b_) - max(a, a_) + 1 for a, b in ints
                if a <= b_ and b >= a_)
        return c / (b_ - a_ + 1)
    def cols_of(sp, a_, b_):
        return {pos2col[sp][p] for p in range(a_, b_ + 1) if p in pos2col[sp]}
    def status(sp, motif, S):
        best_p, best_ov = None, 0.0
        for h in H3[sp]:
            if h['motif'] != motif: continue
            cs = cols_of(sp, h['a'], h['b'])
            ov = len(cs & S) / len(S)
            if ov > 0 and (best_p is None or h['p'] < best_p):
                best_p, best_ov = h['p'], ov
        if best_p is None: return 'lost', None, 0.0
        if best_p < 1e-4 and best_ov >= 0.5: return 'present', best_p, best_ov
        return 'weak', best_p, best_ov
    rows = []
    for h in sorted([x for x in H1['Bdor'] if len(cols_of('Bdor', x['a'], x['b']) & core_cols)], key=lambda x: x['p']):
        S = cols_of('Bdor', h['a'], h['b'])
        st = {}
        for sp in OTH:
            if sp not in aln:
                st[sp] = ('no_homolog', None, 0.0)
            elif chain_cov(sp, h['a'], h['b']) < 0.5:
                st[sp] = ('no_homology', None, chain_cov(sp, h['a'], h['b']))
            else:
                st[sp] = status(sp, h['motif'], S)
        sts = (st['Bmin'][0], st['Btsu'][0])
        if 'present' in sts:
            lin = 'present'
        elif all(s == 'lost' for s in sts):
            lin = 'lost'
        elif 'weak' in sts:
            lin = 'weak'
        else:
            lin = 'no_homology'
        cand = st['Bcor'][0] == 'present' and lin == 'lost'
        rows.append(dict(tag=tag, motif=h['motif'], tf=h['tf'],
                         seg_a=h['a'], seg_b=h['b'], strand=h['strand'],
                         p=f'{h["p"]:.3g}', q=f'{h["q"]:.3g}',
                         **{sp: st[sp][0] for sp in OTH},
                         specialist_lineage_status=lin,
                         candidate='Y' if cand else 'N',
                         cols=cols_of('Bdor', h['a'], h['b'])))
    return rows, len([x for x in H1['Bdor'] if cols_of('Bdor', x['a'], x['b']) & core_cols])

rows_e, n_e = classify_region('E')
rows_c, n_c = classify_region('control')
nc_e = sum(1 for r in rows_e if r['candidate'] == 'Y')
nc_c = sum(1 for r in rows_c if r['candidate'] == 'Y')

flds = ['tag', 'motif', 'tf', 'seg_a', 'seg_b', 'strand', 'p', 'q'] + OTH + \
       ['specialist_lineage_status', 'candidate']
with open(f'{OUT}/Gr28b_E_TFBS_peakanchored.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=flds, delimiter='\t', extrasaction='ignore', lineterminator='\n')
    w.writeheader()
    for r in rows_e + rows_c: w.writerow(r)

with open(f'{OUT}/Gr28b_E_TFBS_counts_peakanchored.tsv', 'w') as f:
    f.write('region\tbdor_hits_p1e-4\tcandidates(Bcor+&specialist-lineage-lost)\n')
    f.write(f'E\t{n_e}\t{nc_e}\ncontrol\t{n_c}\t{nc_c}\n')

# sites on alignment for the figure
with open(f'{OUT}/sites_on_aln_peakanchored.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['motif', 'tf', 'species', 'col_a', 'col_b', 'p', 'status'])
    aln = load_aln(f'{OUT}/E_homologs.mafft.fa')
    for r in rows_e:
        for sp in ['Bdor'] + [s for s in OTH if r[s] not in ('no_homolog', 'no_homology')]:
            if sp not in aln: continue
            for h in fimo_hits('E', sp, '1e-3'):
                if h['motif'] != r['motif']: continue
                p2c_, kk = {}, 0
                for c, ch in enumerate(aln[sp]):
                    if ch != '-': kk += 1; p2c_[kk] = c
                cs = [p2c_[p] for p in range(h['a'], h['b'] + 1) if p in p2c_]
                ov = len(set(cs) & r['cols']) / len(r['cols'])
                if ov > 0:
                    st = 'Bdor_hit' if sp == 'Bdor' else r[sp]
                    w.writerow([r['motif'], r['tf'], sp, min(cs), max(cs), h['p'], st])
print('done: homologs, alignments, ATAC, FIMO, tables ->', OUT)
