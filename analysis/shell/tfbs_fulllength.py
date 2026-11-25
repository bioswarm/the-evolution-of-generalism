#!/usr/bin/env python3
"""Chained full-length identity for peak-anchored Gr28b locus homologs.

For every B. dorsalis query (the 10 mature-ovipositor peaks + control window,
each searched as +/-100 bp) and each of the five other species' Gr28b loci, chain
the E < 1e-3 HSPs and measure the fraction of the *peak interval itself* (not the
flanks) that aligns identically. Also measure the locus background by tiling the
B. dorsalis locus into 400-bp windows.

Writes into <out>: fulllength_identity.tsv, locus_background.tsv,
                   blast/qseq_<sp>.tsv (raw HSPs with aligned sequences).
Usage: tfbs_fulllength.py <out-dir>
"""
import sys, os, csv, subprocess, math

OUT = sys.argv[1]
TDIR = os.path.dirname(OUT)
HIC = '/path/to/data/Hic'
BLAST = '/opt/ncbi-blast-2.10.0+/bin'
GTF = f'{HIC}/Bdor/results_hicexplorer/Bdor.gtf'

# ---------- inputs ----------
sp_locus = {}
for l in open(f'{OUT}/species_locus.tsv'):
    v = l.rstrip('\n').split('\t')
    if not v or not v[0]:
        continue
    sp_locus[v[0]] = dict(fa=v[1], chr=v[2], lo=int(v[3]), hi=int(v[4]),
                          strand=v[5], n1=int(v[6]))
OTH = [s for s in sp_locus if s != 'Bdor']

# flank100 bed: chrom, flank_s, flank_e, name, peak_id, mlog10q, core_s, core_e
queries = []
for l in open(f'{OUT}/queries_flank100.bed'):
    v = l.rstrip('\n').split('\t')
    queries.append(dict(s=int(v[1]), e=int(v[2]), name=v[3],
                        core_s=int(v[6]), core_e=int(v[7])))

# query sequences (flank-extended), in the same order
qseqs = {}
name = None
for l in open(f'{OUT}/queries_flank100.fa'):
    l = l.rstrip('\n')
    if l.startswith('>'):
        name = l[1:].split('::')[0]
        qseqs[name] = ''
    else:
        qseqs[name] += l.strip().upper()

# ---------- blastn with aligned sequences ----------
def run_blast(query_fa, subject, out):
    subprocess.run(
        [f'{BLAST}/blastn', '-task', 'blastn', '-word_size', '7',
         '-evalue', '1e-2', '-dust', 'no', '-soft_masking', 'false',
         '-query', query_fa, '-subject', subject, '-out', out,
         '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend '
                    'sstart send evalue bitscore qseq sseq'],
        check=True, capture_output=True)

def read_hsps(f):
    out = []
    for l in open(f):
        v = l.rstrip('\n').split('\t')
        if len(v) < 14:
            continue
        out.append(dict(q=v[0].split('::')[0], s=v[1], pident=float(v[2]), alen=int(v[3]),
                        qlen=int(v[4]), qs=int(v[6]), qe=int(v[7]),
                        ss=int(v[8]), se=int(v[9]), ev=float(v[10]),
                        bs=float(v[11]), qseq=v[12], sseq=v[13]))
    return out

# ---------- chaining ----------
def chain(hsps):
    """Greedy-by-bitscore collinear chain on the best HSP's subject strand.
    Returns (ordered_hsps, {id: [clip_lo, clip_hi]}, minus_strand)."""
    if not hsps:
        return [], {}, False, None, None
    best = max(hsps, key=lambda h: (h['bs'], -h['ev']))
    minus = best['ss'] > best['se']
    def olo(h):
        lo, hi = min(h['ss'], h['se']), max(h['ss'], h['se'])
        return -hi if minus else lo
    def ohi(h):
        lo, hi = min(h['ss'], h['se']), max(h['ss'], h['se'])
        return -lo if minus else hi
    acc = [best]
    for h in sorted((x for x in hsps if x is not best),
                    key=lambda h: (-h['bs'], h['ev'])):
        if (h['ss'] > h['se']) != minus:
            continue
        # collinearity: query order and oriented-subject order must agree
        if any((h['qs'] - a['qs']) * (olo(h) - olo(a)) < 0 for a in acc):
            continue
        trial = sorted(acc + [h], key=lambda x: x['qs'])
        i = trial.index(h)
        ok = True
        if i > 0:
            p = trial[i - 1]
            if p['qe'] - h['qs'] + 1 > 20 or olo(h) - ohi(p) > 2000:
                ok = False
        if ok and i + 1 < len(trial):
            n = trial[i + 1]
            if h['qe'] - n['qs'] + 1 > 20 or olo(n) - ohi(h) > 2000:
                ok = False
        if ok:
            acc = trial
    acc = sorted(acc, key=lambda x: x['qs'])
    clip = {id(h): [h['qs'], h['qe']] for h in acc}
    for p, n in zip(acc, acc[1:]):
        if p['qe'] >= n['qs']:
            if p['bs'] >= n['bs']:
                clip[id(n)][0] = p['qe'] + 1
            else:
                clip[id(p)][1] = n['qs'] - 1
    return acc, clip, minus, olo, ohi

def measure(ch, clip, core_lo, core_hi, peak_len):
    """identical/covered aligned query bases inside [core_lo,core_hi]."""
    cov = idn = 0
    slo = shi = None
    for h in ch:
        lo, hi = clip[id(h)]
        qpos = h['qs'] - 1
        for qb, sb in zip(h['qseq'], h['sseq']):
            if qb != '-':
                qpos += 1
            else:
                continue
            if lo <= qpos <= hi and core_lo <= qpos <= core_hi:
                cov += 1
                if sb != '-' and sb.upper() == qb:
                    idn += 1
        hlo, hhi = min(h['ss'], h['se']), max(h['ss'], h['se'])
        slo = hlo if slo is None else min(slo, hlo)
        shi = hhi if shi is None else max(shi, hhi)
    return cov, idn, slo, shi

# ---------- peak queries ----------
os.makedirs(f'{OUT}/blast', exist_ok=True)
for sp in OTH:
    run_blast(f'{OUT}/queries_flank100.fa', f'{OUT}/locus_{sp}.fa',
              f'{OUT}/blast/qseq_{sp}.tsv')
qhsps = {sp: read_hsps(f'{OUT}/blast/qseq_{sp}.tsv') for sp in OTH}

# old single-HSP columns (pident/qcov of the best HSP) for recordkeeping
peak_rows = []          # (query, species, fl_cov or '') for the low-coverage check
with open(f'{OUT}/fulllength_identity.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['query', 'species', 'fl_ident', 'fl_cov', 'chain_hsps',
                's_start', 's_end', 'pident', 'qcov', 'q_lo', 'q_hi', 'strand',
                'chain_qints'])
    for q in queries:
        cl = q['core_s'] - q['s'] + 1
        ch = q['core_e'] - q['s']
        plen = q['core_e'] - q['core_s']
        for sp in OTH:
            sl = sp_locus[sp]
            hs = [h for h in qhsps[sp] if h['q'] == q['name'] and h['ev'] < 1e-3]
            chn, clip, minus, olo, ohi = chain(hs)
            if not chn:
                w.writerow([q['name'], sp, '', '', 0, '', '', '', '', '', '', '', ''])
                peak_rows.append((q['name'], sp, '', '', '', ''))
                continue
            best = max(chn, key=lambda h: (h['bs'], -h['ev']))
            cov, idn, slo, shi = measure(chn, clip, cl, ch, plen)
            ints = sorted(clip.values())
            qints = ';'.join(f'{a}-{b}' for a, b in ints)
            fident, fcov = f'{idn / plen * 100:.2f}', f'{cov / plen * 100:.2f}'
            w.writerow([q['name'], sp, fident, fcov,
                        len(chn), sl['lo'] - 1 + slo, sl['lo'] - 1 + shi,
                        f"{best['pident']:.1f}", f"{best['alen']/best['qlen']:.3f}",
                        ints[0][0], ints[-1][1],
                        '-' if minus else '+', qints])
            peak_rows.append((q['name'], sp, fident, fcov,
                              sl['lo'] - 1 + slo, sl['lo'] - 1 + shi))
print('fulllength_identity.tsv written')

# ---------- locus background ----------
locus_fa = None
for l in open(f'{OUT}/locus_Bdor.fa'):
    if l.startswith('>'):
        hdr = l[1:].strip().split(':')[-1]
        lo0 = int(hdr.split('-')[0])
        locus_fa = ''
    else:
        locus_fa += l.strip().upper()
L = len(locus_fa)

# exclusion intervals (locus-fasta 0-based): the 10 peaks + Bdor exons/CDS
excl = [(q['core_s'], q['core_e']) for q in queries if q['name'].startswith('peak_')]
with open(GTF) as g:
    for l in g:
        if not l.startswith('NC_064303.1\t'):
            if l[0] == '#':
                continue
            v = l.split('\t')
            if v[0] != 'NC_064303.1':
                continue
        else:
            v = l.split('\t')
        if v[2] not in ('exon', 'CDS'):
            continue
        s, e = int(v[3]) - 1 - lo0, int(v[4]) - lo0
        if s < L and e > 0:
            excl.append((max(0, s), min(L, e)))

tiles = []
for s in range(0, L - 399, 200):
    seq = locus_fa[s:s + 400]
    if seq.upper().count('N') > 40:
        continue
    if any(a < s + 400 and b > s for a, b in excl):
        continue
    tiles.append((lo0 + s, lo0 + s + 400, seq))

tfa = f'{OUT}/_bg_tiles.fa'
with open(tfa, 'w') as f:
    for i, (s, e, seq) in enumerate(tiles):
        f.write(f'>tile_{i}::{s}-{e}\n')
        for j in range(0, 400, 60):
            f.write(seq[j:j + 60] + '\n')

bg = {sp: [] for sp in OTH}
bg_med = {}
for sp in OTH:
    run_blast(tfa, f'{OUT}/locus_{sp}.fa', f'{OUT}/blast/bg_{sp}.tsv')
    thsps = read_hsps(f'{OUT}/blast/bg_{sp}.tsv')
    per = {}
    for h in thsps:
        per.setdefault(h['q'], []).append(h)
    for i, (s, e, seq) in enumerate(tiles):
        hs = [h for h in per.get(f'tile_{i}', []) if h['ev'] < 1e-3]
        chn, clip, _m, _o, _h = chain(hs)
        if not chn:
            bg[sp].append((f'tile_{i}', s, e, '', '', 0, '', ''))
            continue
        cov, idn, slo, shi = measure(chn, clip, 1, 400, 400)
        bg[sp].append((f'tile_{i}', s, e, f'{idn / 4:.2f}', f'{cov / 4:.2f}',
                       len(chn), slo, shi))
    vals = sorted(float(r[3]) if r[3] != '' else 0.0 for r in bg[sp])
    if vals:
        n = len(vals)
        med = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2
        q1 = vals[n // 4]; q3 = vals[(3 * n) // 4]
        hits = sum(1 for r in bg[sp] if r[3] != '')
        print(f'{sp}: {n} tiles ({hits} aligned), median fl-ident {med:.2f}, '
              f'IQR {q1:.2f}-{q3:.2f}')
        bg_med[sp] = (med, q1, q3)

with open(f'{OUT}/locus_background.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['species', 'tile', 'bd_start', 'bd_end', 'fl_ident', 'fl_cov',
                'chain_hsps', 's_lo', 's_hi'])
    for sp in OTH:
        for r in bg[sp]:
            w.writerow([sp] + list(r))
with open(f'{OUT}/locus_background_summary.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['species', 'n_tiles', 'n_aligned', 'median', 'q1', 'q3'])
    for sp in OTH:
        vals = sorted(float(r[3]) if r[3] != '' else 0.0 for r in bg[sp])
        n = len(vals)
        med = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2
        w.writerow([sp, n, sum(1 for r in bg[sp] if r[3] != ''),
                    f'{med:.2f}', f'{vals[n // 4]:.2f}', f'{vals[(3 * n) // 4]:.2f}'])
os.unlink(tfa)
print('locus_background.tsv written')

# ---------- whole-genome check for low-coverage query x species cells ----------
# For every query x species whose locus chain covers < 50% of the peak, search
# the query against the whole genome and report whether a better collinear hit
# exists on the locus chromosome within 100 kb of the locus (window-artefact check).
need = {(q, sp) for q, sp, _i, c, _s, _e in peak_rows
        if c == '' or float(c) < 50}
if need:
    lowq = sorted({q for q, _s in need})
    qfa = f'{OUT}/_lowq.fa'
    with open(qfa, 'w') as f:
        for qn in lowq:
            f.write(f'>{qn}\n{qseqs[qn]}\n')
    gw = {}
    for sp in sorted({s for _q, s in need}):
        out = f'{OUT}/blast/genome_{sp}.tsv'
        run_blast(qfa, sp_locus[sp]['fa'], out)
        gw[sp] = read_hsps(out)
    os.unlink(qfa)
    qlen = {q['name']: q['e'] - q['s'] for q in queries}
    with open(f'{OUT}/lowcov_genome_check.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['query', 'species', 'locus_fl_cov',
                    'gw_best_chr', 'gw_best_range', 'gw_best_pident',
                    'gw_best_alen', 'gw_best_evalue',
                    'near_locus_hit', 'near_range', 'near_strand_ok',
                    'near_fl_cov', 'near_fl_ident', 'note'])
        for q, sp, fid, cov, _s, _e in peak_rows:
            if (q, sp) not in need:
                continue
            sl = sp_locus[sp]
            qq = next(x for x in queries if x['name'] == q)
            hs = [h for h in gw[sp] if h['q'] == q and h['ev'] < 1e-3]
            best = max(hs, key=lambda h: (h['bs'], -h['ev'])) if hs else None
            near = [h for h in hs if h['s'] == sl['chr'] and
                    min(h['ss'], h['se']) <= sl['hi'] + 100000 and
                    max(h['ss'], h['se']) >= sl['lo'] - 100000]
            nrow = [q, sp, cov or '0']
            if best:
                nrow += [best['s'],
                         f"{min(best['ss'], best['se'])}-{max(best['ss'], best['se'])}",
                         f"{best['pident']:.1f}", best['alen'], f"{best['ev']:.2g}"]
            else:
                nrow += ['', '', '', '', '']
            if near:
                chn, clip, minus, _o, _h = chain(near)
                cl = qq['core_s'] - qq['s'] + 1
                chh = qq['core_e'] - qq['s']
                plen = qq['core_e'] - qq['core_s']
                cov2, idn2, slo2, shi2 = measure(chn, clip, cl, chh, plen)
                strand_ok = (not minus) == (sl['strand'] == '+')
                inwin = slo2 >= sl['lo'] and shi2 <= sl['hi']
                # hits in the 20 kb added to the 3' margin would have been
                # missed by the previous (10 kb margin) pipeline
                was_out = (shi2 > sl['hi'] - 20000) if sl['strand'] == '+' \
                    else (slo2 < sl['lo'] + 20000)
                note = ('inside new window' if inwin else
                        'near locus but outside window') + \
                    ('; missed by old 10 kb margin' if was_out else '')
                nrow += ['Y', f"{slo2}-{shi2}", 'Y' if strand_ok else 'N',
                         f'{cov2 / plen * 100:.2f}', f'{idn2 / plen * 100:.2f}', note]
            else:
                nrow += ['N', '', '', '', '', 'no hit on locus chr +/-100 kb']
            w.writerow(nrow)
    print('lowcov_genome_check.tsv written')
