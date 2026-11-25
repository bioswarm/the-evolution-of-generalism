#!/usr/bin/env python3
"""Supplementary Data 7 TFBS section from the peak-anchored analysis.
Reads <out>/Gr28b_E_TFBS_peakanchored.tsv, identities_E.tsv,
results/locus/tfbs/accessibility.tsv; writes
<out>/SuppData7_section_TFBS.tsv  (Bdor E hits deduplicated + footer).
Expression columns are joined exactly as in tfbs_classify.py.
"""
import sys, os, re, csv

OUT = sys.argv[1]
TDIR = os.path.dirname(OUT)

rows = [r for r in csv.DictReader(open(f'{OUT}/Gr28b_E_TFBS_peakanchored.tsv'), delimiter='\t')
        if r['tag'] == 'E']
counts = {l.split('\t')[0]: l.split('\t')[1:] for l in open(f'{OUT}/Gr28b_E_TFBS_counts_peakanchored.tsv')
          if not l.startswith('region')}
ide = {}
for l in open(f'{OUT}/identities_E.tsv'):
    v = l.rstrip('\n').split('\t')
    if v[0] != 'species': ide[v[0]] = (v[3], v[4])      # id_homologous, ncols
acc = {}
for l in open(f'{TDIR}/accessibility.tsv'):
    v = l.rstrip('\n').split('\t')
    if v[0] != 'track':
        acc[v[0] + ('@E' if ':108181027' in v[1] else '@CT')] = v[2]

# ---- expression join (identical to tfbs_classify.py) ----
TPM = '/path/to/data/genesets/bactrocera/orthologue/new/gene_expression_wide_output.txt'
hdr = open(TPM).readline().rstrip('\n').split('\t'); cols = {c: i for i, c in enumerate(hdr)}
want = [f'Bd-{s}-O{k}' for s in ('FN', 'FM') for k in (1, 2, 3)] + [f'Bm-FM-O{k}' for k in (1, 2, 3)]
tp = {}
for i, l in enumerate(open(TPM)):
    if i == 0: continue
    v = l.rstrip('\n').split('\t')
    tp[v[0].split('_')[0]] = {w: float(v[cols[w]]) for w in want if w in cols}
fly = {}
for l in open('/path/to/data/genesets/bactrocera/kallisto/results/mashR/Bdor-fly.genenames.txt'):
    v = l.rstrip('\n').split('\t')
    if len(v) >= 2: fly[v[1].lower()] = v[0]
orth = {}
for i, l in enumerate(open('/path/to/data/genesets/bactrocera/kallisto/results/mashR/bactrocera_orthologues.txt')):
    if i == 0: continue
    v = l.rstrip('\n').split('\t')
    orth[v[0]] = dict(zip(['Bcor', 'Zcuc', 'Bmin', 'Ztau', 'Btsu'], v[1:6]))
SYNONYM = {'vfl': 'zld'}
def bd_gene(tf):
    t = tf.lower()
    keys = (t, t.split('::')[0], re.sub(r'[^a-z0-9].*$', '', t))
    for cand in keys + tuple(SYNONYM.get(k, '') for k in keys):
        if cand and cand in fly: return fly[cand]
    return ''
def tpm(g, pre):
    return round(sum(tp[g][k] for k in want if k.startswith(pre)) / 3, 3) if g in tp else ''

# ---- dedupe strand duplicates, genomic coords ----
GS0 = 108180926                     # E query segment start - 1 (peak 108181027 - 100 flank, 0-based bed start)
grp = {}
for r in rows:
    a, b = GS0 + int(r['seg_a']), GS0 + int(r['seg_b'])
    k = (r['motif'], a, b)
    g = grp.setdefault(k, {'row': r, 'str': set(), 'a': a, 'b': b})
    g['str'].add(r['strand'])
    if float(r['p']) < float(g['row']['p']): g['row'] = r
EX1 = 108181492
nmotif = sum(1 for l in open(f'{TDIR}/jaspar_insects.meme') if l.startswith('MOTIF'))
with open(f'{OUT}/SuppData7_section_TFBS.tsv', 'w') as f:
    f.write('jaspar_id\ttf\tbd_start\tbd_end\tstrand\tp\tq\tdist_to_exon1_bp\t'
            'Bcor\tZcuc\tZtau\tspecialist_lineage\tcandidate\tbd_fn_tpm\tbd_fm_tpm\tbm_fm_tpm\n')
    for k in sorted(grp, key=lambda x: x[1]):
        r, st = grp[k]['row'], grp[k]['str']
        g = bd_gene(r['tf'])
        bm = orth.get(g, {}).get('Bmin', '')
        f.write('\t'.join(map(str, [r['motif'], r['tf'], grp[k]['a'], grp[k]['b'],
            '+/-' if len(st) > 1 else r['strand'], r['p'], r['q'], EX1 - grp[k]['b'],
            r['Bcor'], r['Zcuc'], r['Ztau'], r['specialist_lineage_status'],
            r['candidate'], tpm(g, 'Bd-FN'), tpm(g, 'Bd-FM'),
            tpm(g, 'Bm-FM') if bm else ''])) + '\n')
    n_e, nc_e = counts['E']; n_c, nc_c = counts['control']
    n_usite = len({(r['motif'], r['seg_a']) for r in rows if r['candidate'] == 'Y'})
    f.write(f'# motifs scanned {nmotif}; Bdor E hits (p<1e-4) {n_e}; '
            f'E candidates {nc_e} ({n_usite} unique sites); '
            f'control window hits {n_c}, candidates {nc_c} '
            f'(control has no homolog outside the dorsalis lineage, so specialist status is '
            f'no_homology and no site qualifies as a specialist loss)\n')
    f.write(f'# mean fold-enrichment: E MO {acc["Bdor_MO@E"]} / NO {acc["Bdor_NO@E"]}; '
            f'control MO {acc["Bdor_MO@CT"]} / NO {acc["Bdor_NO@CT"]}\n')
    f.write('# identity to B. dorsalis over the BLAST-homologous columns of the E alignment '
            '(peak-anchored L-INS-i, restricted to each species\' chained HSPs): ' +
            '; '.join(f'{sp} {ide.get(sp, ("NA",))[0]} over {ide.get(sp, ("", "NA"))[1]} cols'
                      for sp in ('Bcor', 'Bmin', 'Btsu', 'Zcuc', 'Ztau')) +
            '; specialist homologs cover only the exon-proximal ~270 bp\n')
    f.write('# status "no_homology": the site\'s alignment columns lie outside the species\' '
            'chained BLASTN HSPs, so there is no detectable similarity and motif loss cannot '
            'be assessed; candidates require present in B. correcta AND lost in both '
            'B. minax and B. tsuneonis\n')
print('wrote', f'{OUT}/SuppData7_section_TFBS.tsv')
