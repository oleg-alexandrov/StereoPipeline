#!/usr/bin/env python3
"""
wv3_green_combine.py - combine several acquisitions' green->PAN per-column
disparities into one shipped CCD correction for a single (satellite, TDI,
scan-direction) cell, for the WorldView-3 green band (band 3). This is the
"combine + detrend" step; see README_MULTISPECTRAL for the full recipe.

INPUT : one <avg-dx.txt> per acquisition (the column-averaged green->PAN
        disparity from the measure step). The matching <avg-dy.txt> is loaded
        automatically (same name with dx -> dy). All files must have the same
        length (the green image width, e.g. 10651).
OUTPUT: <out-prefix>-dx.txt and <out-prefix>-dy.txt, the combined, detrended,
        edge-clamped per-column correction, ready for form_corrections_image.py.

Refinements over the older ccd_process.py (2026):
  * nan-safe MEAN across acquisitions (a whole invalid column is left NaN);
  * a MOVING-AVERAGE high-pass detrend (window ~1.8x the CCD sub-array width):
    it removes the broad parallax / terrain / boresight trend while KEEPING the
    per-column CCD steps, centered on zero;
  * EDGE CLAMPING: the outer columns at each end are off-border correlation =
    garbage; they are forced flat at a clean interior value, otherwise wv_correct
    would shift the first/last image columns by the garbage amount.
Before combining, the script prints each acquisition's interior std - scan it and
drop any obvious outlier (a scene whose correlation failed can inflate the mean).

Example:
  python wv3_green_combine.py cell14fwd_corr \\
    run1/avg-dx.txt run2/avg-dx.txt run3/avg-dx.txt --win 1600
"""
import sys, argparse
import numpy as np
from scipy.ndimage import uniform_filter1d

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("out_prefix", help="output prefix; writes <prefix>-dx.txt and -dy.txt")
    p.add_argument("dx_files", nargs="+", help="per-acquisition avg-dx.txt files (dy inferred)")
    p.add_argument("--win", type=int, default=1600,
                   help="moving-average detrend window in columns (default 1600, ~1.8x the CCD sub-array)")
    p.add_argument("--edge", type=int, default=150,
                   help="columns held flat past the last valid column at each end (default 150)")
    p.add_argument("--interior", type=int, default=700,
                   help="columns trimmed at each end when reporting the interior std (default 700)")
    return p.parse_args()

def garbage_mask(dx):
    # Off-border correlation gives wild dx; flag columns far from the robust median.
    med = np.nanmedian(dx)
    mad = np.nanmedian(np.abs(dx - med)) + 1e-9
    return ~np.isfinite(dx) | (np.abs(dx - med) > 8 * mad)

def movavg_detrend(a, bad, win):
    v = a.copy()
    v[bad] = np.nanmedian(a[~bad])           # fill garbage with the median so the filter is stable
    trend = uniform_filter1d(v, size=win, mode="nearest")
    return a - trend

def clamp_edges(a, bad, edge):
    valid = np.where(~bad)[0]
    if valid.size == 0:
        return a
    lo, hi = valid[0] + edge, valid[-1] - edge + 1
    out = a.copy()
    out[:lo] = np.nanmedian(a[lo:lo + 100])  # flat at a clean interior value, not the edge garbage
    out[hi:] = np.nanmedian(a[hi - 100:hi])
    return out

def combine(files, win, edge, interior):
    arrs = [np.loadtxt(f) for f in files]
    n = arrs[0].size
    for f, a in zip(files, arrs):
        if a.size != n:
            sys.exit(f"ERROR: {f} has length {a.size}, expected {n}")
        s = np.nanstd((a - np.nanmedian(a))[interior:-interior])
        print(f"  {f}: interior std = {s:.4f}   (scan for outliers)")
    grand = np.nanmean(np.vstack(arrs), axis=0)
    bad = garbage_mask(grand)
    corr = movavg_detrend(grand, bad, win)
    corr = clamp_edges(corr, bad, edge)
    return corr

def main():
    o = parse_args()
    dy_files = [f.replace("dx", "dy") for f in o.dx_files]
    print(f"Combining {len(o.dx_files)} acquisitions (win={o.win}):")
    print(" dx:")
    cx = combine(o.dx_files, o.win, o.edge, o.interior)
    print(" dy:")
    cy = combine(dy_files, o.win, o.edge, o.interior)
    np.savetxt(o.out_prefix + "-dx.txt", cx)
    np.savetxt(o.out_prefix + "-dy.txt", cy)
    i = o.interior
    print(f"Wrote {o.out_prefix}-dx.txt (interior std {np.nanstd(cx[i:-i]):.4f}) "
          f"and {o.out_prefix}-dy.txt (interior std {np.nanstd(cy[i:-i]):.4f})")
    print(f"Edges: dx[0]={cx[0]:.4f} dx[-1]={cx[-1]:.4f}  dy[0]={cy[0]:.4f} dy[-1]={cy[-1]:.4f} "
          f"(must be near 0, not near the boresight offset)")

if __name__ == "__main__":
    main()
