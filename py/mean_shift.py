#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import numpy as np
from math import ceil,floor
from statistics import mean

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-c",
                        "--raw-coverage",
                        type=str,
                        required=True,
                        help="Path to raw coverage TXT file")
    parser.add_argument("-bs",
                        "--histogram-bin-size",
                        type=float,
                        default=0.025,
                        help="Bin size for the histogram used for coverage curve fitting")
    parser.add_argument("-w",
                        "--window-size",
                        type=int,
                        default=1,
                        help="Non-overlapping window size")
    args = parser.parse_args()
    return args

def gaussian(x, a, x0, sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def ms_gradient(i, j, r_i, r_j, H_b, H_r,):
    e_b = gaussian(x=i, a=2, x0=j, sigma=H_b)
    e_r = gaussian(x=r_i, a=2, x0=r_j, sigma=H_r)
    return (j-i)*e_b*e_r

def get_neg_pos_breaks(l):
    result = list()
    for i in range(1, len(l)):
        if l[i] > 0 and not l[i-1] > 0:
            result.append(i)
    return result

def plot_coverage(coverage, ticks=list(), predicted_junctions=list(), outpath='coverage_plot.pdf'):
    plt.figure(figsize=(30,4))
    plt.plot(range(len(coverage)), coverage)
    for tick in ticks:
        plt.plot([tick,tick], [0,1], 'k--', alpha=0.2)
    for tick in predicted_junctions:
        plt.plot([tick,tick], [0,1], 'g', alpha=0.2)
    plt.tight_layout()
    plt.savefig(outpath)

def get_segments(coverage_l=0, ticks=set()):
    segs = list()
    segs.append((
        0,
        coverage_l,
    ))
    for i in range(coverage_l):
        if i in ticks:
            segs[-1] = (
                segs[-1][0],
                i,
            )
            segs.append((
                i,
                coverage_l,
            ))
    return segs

def get_variation(coverage, step):
    x = ar(list(np.arange(step/2, 1-step/2, step)))
    y = ar(
        np.histogram(
            [c for c in coverage if float(c) > 0.05],
            bins=np.arange(0,1,step)
        )[0]
    )
    popt,pcov = curve_fit(gaussian,x,y,maxfev=50000)
    return popt[2]**2

def get_hb_ticks(coverage, coverage_var, window_size, Hb_list, segmentations_count=3, lim=200):
    ticks = {x for x in range(window_size, len(coverage), window_size)}
    Hb_ticks = list()
    for H_b in Hb_list:
        segments = get_segments(len(coverage), ticks)
        for round_count in range(segmentations_count):
            gradients = list()
            for i,seg in enumerate(segments):
                i_pos = floor((seg[0]+seg[1])/2)
                i_rd = mean(coverage[seg[0]:seg[1]])
    #             print(i, seg, i_pos, i_rd)
                start = max(0,i-lim)
                end = min(i+lim, len(segments)-1)
                S = 0
                for j in range(start, end+1):
                    j_pos = floor((segments[j][0]+segments[j][1])/2)
                    j_rd = mean(coverage[segments[j][0]:segments[j][1]])
                    S += ms_gradient(i=i_pos, j=j_pos, r_i=i_rd, r_j=j_rd, H_b=H_b, H_r=coverage_var)
                gradients.append(S)
            new_segments = list()
            [(0,len(coverage))]

            ticks = {
                segments[seg_idx][0] for seg_idx in get_neg_pos_breaks(gradients)
            }
            segments = get_segments(len(coverage), ticks)
            print('Hb = {} Round = {} |ticks| = {}'.format(H_b, round_count, len(ticks)))
            print(ticks)
            Hb_ticks.append((
                ticks,
                'Hb = {}; Round = {}'.format(H_b, round_count)
            ))
    return Hb_ticks

def plot_coverage(coverage, Hb_ticks, out_path):
    subplot_count = len(Hb_ticks)
    f, ax = plt.subplots(subplot_count, 1, sharex=True, figsize=(30,5*subplot_count))
    for idx,(ticks,text) in enumerate(Hb_ticks):
        ax[idx].plot(range(len(coverage)), coverage)
        ax[idx].title.set_text(text)
        for tick in ticks:
            ax[idx].plot([tick,tick], [0,1], 'g', alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path)

def main():
    args = parse_args()

    coverage = [float(c) for c in open(args.raw_coverage).readlines()]
    coverage_var = get_variation(coverage=coverage, step=args.histogram_bin_size)

    Hb_list = [2,4,8,16,32,64,128]
    segmentations_count=3
    lim=200
    Hb_ticks = get_hb_ticks(coverage=coverage, coverage_var=coverage_var, window_size=args.window_size, Hb_list=Hb_list, segmentations_count=segmentations_count, lim=lim)

    outpath = '{}.meanshift_coverage.pdf'.format(args.raw_coverage[0:args.raw_coverage.rfind('.')])
    plot_coverage(coverage=coverage, Hb_ticks=Hb_ticks, out_path=outpath)

if __name__ == "__main__":
    main()
