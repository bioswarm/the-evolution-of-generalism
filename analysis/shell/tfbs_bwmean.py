#!/usr/bin/env python
"""Mean fold-enrichment of a bigWig over one or more regions.
Usage: tfbs_bwmean.py <bigwig> <chr>:<start>-<end> [<chr>:<start>-<end> ...]
Prints: region <TAB> mean <TAB> n_bases_with_data
Runs under the chipseq env (pyBigWig).
"""
import sys, pyBigWig

bw = pyBigWig.open(sys.argv[1])
for reg in sys.argv[2:]:
    c, se = reg.split(':')
    s, e = (int(x) for x in se.replace(',', '').split('-'))
    v = bw.values(c, s, e, numpy=True)
    import numpy as np
    ok = np.isfinite(v)
    print(f'{c}:{s}-{e}\t{np.nanmean(v):.4f}\t{int(ok.sum())}')
bw.close()
