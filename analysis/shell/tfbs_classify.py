#!/usr/bin/env python3
"""TFBS classification for the Gr28b C-clade upstream element E.
Usage: tfbs_classify.py <out_dir> <E_region> <CT_region> <tmp_dir>
  <out_dir>   = results/locus/tfbs (contains species.tsv, upstream_Cexon.mafft.fa,
                fimo/<sp>_{1e-4,1e-3}/fimo.tsv, jaspar_insects.meme)
  E_region / CT_region = 1-based inclusive ranges inside the 1600-bp extracted
  sequence, e.g. 1036-1463 / 608-1035 (Bdor coordinates).
Writes: identities.tsv, Gr28b_Cexon_upstream_TFBS.tsv, sites_on_aln.tsv,
        Gr28b_Cexon_upstream_TFBS_counts.tsv, tf_names.txt
If snrna_detection.tsv exists, its pct column is merged into the main table.
"""
import sys, os, re, csv, glob

OUT, E_RG, CT_RG, TMP = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
E_A, E_B = map(int, E_RG.split('-')); C_A, C_B = map(int, CT_RG.split('-'))
SPP = [l.split('\t')[0] for l in open(f'{OUT}/species.tsv') if l.strip()]
OTH = [s for s in SPP if s != 'Bdor']

# ---------- alignment ----------
aln = {}
name = None
for l in open(f'{OUT}/upstream_Cexon.mafft.fa'):
    if l.startswith('>'): name = l[1:].split()[0]; aln[name] = ''
    else: aln[name] += l.strip()
L = len(aln['Bdor'])
# position (1-based ungapped) -> column and column -> position per species
pos2col, col2pos = {}, {}
for sp, s in aln.items():
    p2c, c2p, k = {}, {}, 0
    for c, ch in enumerate(s):
        if ch != '-':
            k += 1; p2c[k] = c; c2p[c] = k
    pos2col[sp], col2pos[sp] = p2c, c2p
E_cols  = {pos2col['Bdor'][p] for p in range(E_A, E_B + 1)}
CT_cols = {pos2col['Bdor'][p] for p in range(C_A, C_B + 1)}

def ident(sp, cols):
    bd = aln['Bdor']; s = aln[sp]
    m = sum(1 for c in cols if s[c] == bd[c])
    ng = sum(1 for c in cols if s[c] != '-')
    return m / len(cols), ng

# ---------- synteny segments ----------
syn_beds = {'Bcor': 'Bdor_vs_Bcor', 'Bmin': 'Bdor_vs_Bmin', 'Zcuc': 'Bdor_vs_Bcuc',
            'Ztau': 'Bdor_vs_Btau', 'Btsu': None}
SYND = '/path/to/data/Hic/Bdor/results_hicexplorer'
def in_syn(sp, a, b):
    if syn_beds.get(sp) is None: return 'NA'
    for i, l in enumerate(open(f'{SYND}/{syn_beds[sp]}.gui-gapped1000.genomic.bed')):
        if i == 0 or not l.strip(): continue
        v = l.split('\t')
        if v[0] == 'NC_064303.1' and int(v[1]) < b and int(v[2]) > a:
            return 'Y'
    return 'N'
# genomic E / CT spans (Bdor + strand; region file started at 108179992)
GS0 = 108179992
E_G, CT_G = (GS0 + E_A - 1, GS0 + E_B - 1), (GS0 + C_A - 1, GS0 + C_B - 1)

with open(f'{OUT}/identities.tsv', 'w') as f:
    f.write('species\tid_E\tngap_E\tid_CT\tngap_CT\tE_in_synteny\tCT_in_synteny\n')
    for sp in SPP:
        ie, ne = ident(sp, E_cols); ic, nc = ident(sp, CT_cols)
        f.write(f'{sp}\t{ie:.4f}\t{ne}\t{ic:.4f}\t{nc}\t'
                f'{in_syn(sp, *E_G)}\t{in_syn(sp, *CT_G)}\n')

# ---------- motif names ----------
motif_name = {}
for m in re.finditer(r'MOTIF (\S+)(?: +(\S+))?', open(f'{OUT}/jaspar_insects.meme').read()):
    motif_name[m.group(1)] = m.group(2) or m.group(1)

# ---------- FIMO hits ----------
def hits(sp, thr):
    f = f'{OUT}/fimo/{sp}_{thr}/fimo.tsv'
    h = []
    if not os.path.exists(f): return h
    for l in open(f):
        if l.startswith('#') or l.startswith('motif_id'): continue
        v = l.rstrip('\n').split('\t')
        if len(v) < 9: continue
        h.append(dict(motif=v[0], tf=motif_name.get(v[0], v[1]), a=int(v[3]), b=int(v[4]),
                      strand=v[5], p=float(v[7]), q=float(v[8]) if v[8] not in ('', '-') else 1.0))
    return h
H1, H3 = {sp: hits(sp, '1e-4') for sp in SPP}, {sp: hits(sp, '1e-3') for sp in SPP}

def cols_of(sp, a, b):
    p2c = pos2col[sp]
    return {p2c[p] for p in range(a, b + 1) if p in p2c}

def region_hits(hl, cols):
    out = []
    for h in hl:
        cs = cols_of('Bdor', h['a'], h['b'])
        ov = len(cs & cols) / len(cols)
        if ov > 0:
            h = dict(h); h['cols'] = cs; h['ov_region'] = ov; out.append(h)
    out.sort(key=lambda h: h['p'])
    return out

def status(sp, motif, S):
    """present / lost / weak for one other species at Bdor site columns S"""
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

# ---------- expression data ----------
def rd_tsv(f):
    return list(csv.DictReader(open(f), delimiter='\t'))
TPM = f'/path/to/data/genesets/bactrocera/orthologue/new/gene_expression_wide_output.txt'
tp = {}
hdr = open(TPM).readline().rstrip('\n').split('\t')
cols = {c: i for i, c in enumerate(hdr)}
want = [f'Bd-{s}-O{k}' for s in ('FN', 'FM') for k in (1, 2, 3)] + [f'Bm-FM-O{k}' for k in (1, 2, 3)]
for i, l in enumerate(open(TPM)):
    if i == 0: continue
    v = l.rstrip('\n').split('\t')
    g = v[0].split('_')[0]
    tp[g] = {w: float(v[cols[w]]) for w in want if w in cols}
fly = {}
for l in open('/path/to/data/genesets/bactrocera/kallisto/results/mashR/Bdor-fly.genenames.txt'):
    v = l.rstrip('\n').split('\t')
    if len(v) >= 2: fly[v[1].lower()] = v[0]
orth = {}
for i, l in enumerate(open('/path/to/data/genesets/bactrocera/kallisto/results/mashR/bactrocera_orthologues.txt')):
    if i == 0: continue
    v = l.rstrip('\n').split('\t')
    orth[v[0]] = dict(zip(['Bcor', 'Zcuc', 'Bmin', 'Ztau', 'Btsu'], v[1:6]))
pat = {r['gene']: float(r['Bd']) for r in rd_tsv(f'{OUT}/../../mashr/maturation/patterns_ovipositor.tsv')}
lfs = {r['gene']: float(r['Bd']) for r in rd_tsv(f'{OUT}/../../mashr/maturation/lfsr_ovipositor.tsv')}

SYNONYM = {'vfl': 'zld'}          # JASPAR legacy TF name -> current Dmel gene name

def bd_gene(tf):
    t = tf.lower()
    keys = (t, t.split('::')[0], re.sub(r'[^a-z0-9].*$', '', t))
    for cand in keys + tuple(SYNONYM.get(k, '') for k in keys):
        if cand and cand in fly: return fly[cand]
    return ''

def tpm(g, pre):
    return round(sum(tp[g][k] for k in want if k.startswith(pre)) / 3, 3) if g in tp else ''

sn = {}
sd = f'{OUT}/snrna_detection.tsv'
if os.path.exists(sd):
    sn = {r['tf']: r['pct_Gr28bplus'] for r in rd_tsv(sd)}

# ---------- classify ----------
def classify(cols, tag):
    bh = region_hits(H1['Bdor'], cols)          # Bdor hits p<1e-4 in the region
    rows = []
    for h in bh:
        S = h['cols']
        st = {sp: status(sp, h['motif'], S) for sp in OTH}
        cand = st['Bcor'][0] == 'present' and st['Bmin'][0] == 'lost' and st['Btsu'][0] == 'lost'
        g = bd_gene(h['tf'])
        bm = orth.get(g, {}).get('Bmin', '')
        rows.append(dict(
            tag=tag, motif=h['motif'], tf=h['tf'],
            bd_start=GS0 + h['a'] - 1, bd_end=GS0 + h['b'] - 1, strand=h['strand'],
            p=f'{h["p"]:.3g}', q=f'{h["q"]:.3g}', ov_region=f'{h["ov_region"]:.2f}',
            **{sp: st[sp][0] for sp in OTH},
            **{f'{sp}_p': f'{st[sp][1]:.2g}' if st[sp][1] else '' for sp in OTH},
            bd_gene=g, bd_fn_tpm=tpm(g, 'Bd-FN'), bd_fm_tpm=tpm(g, 'Bd-FM'),
            mashr_log2fc=f'{pat[g]:.3f}' if g in pat else '',
            mashr_lfsr=f'{lfs[g]:.2e}' if g in lfs else '',
            bm_ortho=bm, bm_fm_tpm=tpm(g, 'Bm-FM') if bm else '',   # orthology-keyed table: Bm cols live on the Bd gene row
            snrna_pct=sn.get(h['tf'], ''), candidate='Y' if cand else ''))
    nc = sum(1 for r in rows if r['candidate'] == 'Y')
    return rows, len(bh), nc

rows_e, n_e, nc_e = classify(E_cols, 'E')
rows_c, n_c, nc_c = classify(CT_cols, 'control')
rows = rows_e + rows_c

with open(f'{OUT}/Gr28b_Cexon_upstream_TFBS.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', extrasaction='ignore', lineterminator='\n')
    w.writeheader(); w.writerows(rows)

# sites table for the figure: every same-motif hit (p<1e-3) overlapping a Bdor E-site
sitecols = {}
for r in rows_e:
    sitecols.setdefault(r['motif'], []).append((r['bd_start'], r['bd_end']))
with open(f'{OUT}/sites_on_aln.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['motif', 'tf', 'species', 'col_a', 'col_b', 'p', 'strand', 'status'])
    for r in rows_e:
        S = cols_of('Bdor', r['bd_start'] - GS0 + 1, r['bd_end'] - GS0 + 1)
        for sp in SPP:
            sts = 'Bdor_hit' if sp == 'Bdor' else status(sp, r['motif'], S)[0]
            for h in H3[sp]:
                if h['motif'] != r['motif']: continue
                cs = cols_of(sp, h['a'], h['b'])
                if len(cs & S) / len(S) > 0 or (sp == 'Bdor' and cs == S):
                    w.writerow([r['motif'], r['tf'], sp, min(cs), max(cs), f'{h["p"]:.3g}',
                                h['strand'], sts])
w_counts = open(f'{OUT}/Gr28b_Cexon_upstream_TFBS_counts.tsv', 'w')
w_counts.write('region\tn_bdor_hits_p1e-4\tn_candidates\n')
w_counts.write(f'E\t{n_e}\t{nc_e}\ncontrol\t{n_c}\t{nc_c}\n')
w_counts.write(f'total_motifs_scanned\t{len(motif_name)}\t\n')
w_counts.close()

# ---------- Supplementary Data 7 section: Bdor E hits, deduplicated ----------
EX1 = 108181492                     # Bd1G02958.1 exon-1 start (all sites are upstream)
grp = {}
for r in rows_e:
    k = (r['motif'], int(r['bd_start']), int(r['bd_end']))
    g = grp.setdefault(k, {'row': r, 'str': set()})
    g['str'].add(r['strand'])
    if float(r['p']) < float(g['row']['p']): g['row'] = r
acc = {}
for l in open(f'{OUT}/accessibility.tsv'):
    v = l.rstrip('\n').split('\t')
    if v[0] != 'track':
        acc[v[0] + ('@E' if ':108181027' in v[1] else '@CT')] = v[2]
ide = {}
for l in open(f'{OUT}/identities.tsv'):
    v = l.rstrip('\n').split('\t')
    if v[0] != 'species':
        ide[v[0]] = (v[1], v[3])
with open(f'{OUT}/SuppData7_section_TFBS.tsv', 'w') as f:
    f.write('jaspar_id\ttf\tbd_start\tbd_end\tstrand\tp\tq\tdist_to_exon1_bp\t'
            'Bcor\tZcuc\tZtau\tBmin\tBtsu\tcandidate\tbd_fn_tpm\tbd_fm_tpm\tbm_fm_tpm\n')
    for k in sorted(grp, key=lambda x: x[1]):
        r, st = grp[k]['row'], grp[k]['str']
        f.write('\t'.join(map(str, [r['motif'], r['tf'], r['bd_start'], r['bd_end'],
            '+/-' if len(st) > 1 else r['strand'], r['p'], r['q'], EX1 - int(r['bd_end']),
            r['Bcor'], r['Zcuc'], r['Ztau'], r['Bmin'], r['Btsu'],
            'Y' if r['candidate'] == 'Y' else 'N',
            r['bd_fn_tpm'], r['bd_fm_tpm'], r['bm_fm_tpm']])) + '\n')
    n_usite = len({r['motif'] for r in rows_e if r['candidate'] == 'Y'})
    f.write(f'# motifs scanned {len(motif_name)}; Bdor E hits (p<1e-4) {n_e}; '
            f'E candidates {nc_e} ({n_usite} unique sites); '
            f'control window hits {n_c}, candidates {nc_c}\n')
    f.write(f'# mean fold-enrichment: E MO {acc["Bdor_MO@E"]} / NO {acc["Bdor_NO@E"]}; '
            f'control MO {acc["Bdor_MO@CT"]} / NO {acc["Bdor_NO@CT"]}\n')
    f.write('# identity to Bdor over E columns (control columns): ' +
            '; '.join(f'{sp} {ide[sp][0]} ({ide[sp][1]})' for sp in ('Bcor', 'Zcuc', 'Ztau', 'Bmin', 'Btsu')) + '\n')

# tf name list for the snRNA detection step
with open(f'{OUT}/tf_names.txt', 'w') as f:
    f.write('\n'.join(sorted({r['tf'] for r in rows})))

print(f'Bdor E hits (p<1e-4): {n_e}; candidates: {nc_e}; control: {n_c} hits / {nc_c} candidates')
print('identities:', [l.strip() for l in open(f"{OUT}/identities.tsv").readlines()[1:]])
