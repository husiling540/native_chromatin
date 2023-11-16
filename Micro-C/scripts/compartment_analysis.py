#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
import cooler
import cooltools.lib.plotting
import sys
import cooltools
import warnings
from cytoolz import merge
from packaging import version

if version.parse(cooltools.__version__) < version.parse('0.5.4'):
    raise AssertionError("tutorials rely on cooltools version 0.5.4 or higher,"+
                         "please check your cooltools version and update to the latest")

sample_name = sys.argv[1]

cool_file_path = "/mnt/disk3/husiling/private/data/native_chromatin_Micro-C/20230630_K562_Micro-C_JIACE_BJ_data_analysis/4hicpro/K562_MicroC_allvalidpairs/cooler_res100000"
clr = cooler.Cooler(f'{cool_file_path}/K562_MicroC_{sample_name}.mcool::resolutions/100000')

if not os.path.isfile('/mnt/disk3/husiling/private/pipeline/Hi-C/saddle_plot/hg38.fa'):
    subprocess.call('wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O /mnt/disk3/husiling/private/pipeline/Hi-C/saddle_plot/hg38.fa.gz', shell=True)
    subprocess.call('gunzip /mnt/disk3/husiling/private/pipeline/Hi-C/saddle_plot/hg38.fa.gz', shell=True)

import bioframe
bins = clr.bins()[:]
hg38_genome = bioframe.load_fasta('/mnt/disk3/husiling/private/pipeline/Hi-C/saddle_plot/hg38.fa');

gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)
gc_cov.to_csv('hg38_gc_cov_100kb.tsv',index=False,sep='\t')

view_df = pd.DataFrame({'chrom': clr.chromnames,
                        'start': 0,
                        'end': clr.chromsizes.values,
                        'name': clr.chromnames}
                      )

cis_eigs = cooltools.eigs_cis(
                        clr,
                        gc_cov,
                        view_df=view_df,
                        n_eigs=3,
                        )

eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
eigenvector_track.to_csv(f'Compartment_PC1_{sample_name}.tsv',index=False,sep='\t')

if os.path.exists(f'Expected_cis_100kb_{sample_name}.csv'):
    cvd = pd.read_csv(f'Expected_cis_100kb_{sample_name}.csv')
else:
    cvd = cooltools.expected_cis(clr=clr, view_df=view_df)
    cvd.to_csv(f'Expected_cis_100kb_{sample_name}.csv')

Q_LO = 0.02
Q_HI = 0.98
N_GROUPS = 48


interaction_sum, interaction_count =  cooltools.saddle(
        clr,
        cvd,
        eigenvector_track,
        'cis',
        n_bins=N_GROUPS,
        qrange=(Q_LO,Q_HI),
        view_df=view_df
)

def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    strength_label=None,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    """
    Generate a saddle plot.
    Parameters
    ----------
    track : pd.DataFrame
        See cooltools.digitize() for details.
    saddledata : 2D array-like
        Saddle matrix produced by `make_saddle`. It will include 2 flanking
        rows/columns for outlier signal values, thus the shape should be
        `(n+2, n+2)`.
    cmap : str or matplotlib colormap
        Colormap to use for plotting the saddle heatmap
    scale : str
        Color scaling to use for plotting the saddle heatmap: log or linear
    vmin, vmax : float
        Value limits for coloring the saddle heatmap
    color : matplotlib color value
        Face color for margin bar plots
    fig : matplotlib Figure, optional
        Specified figure to plot on. A new figure is created if none is
        provided.
    fig_kws : dict, optional
        Passed on to `plt.Figure()`
    heatmap_kws : dict, optional
        Passed on to `ax.imshow()`
    margin_kws : dict, optional
        Passed on to `ax.bar()` and `ax.barh()`
    cbar_kws : dict, optional
        Passed on to `plt.colorbar()`
    subplot_spec : GridSpec object
        Specify a subregion of a figure to using a GridSpec.
    Returns
    -------
    Dictionary of axes objects.
    """

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib import ticker
    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]
    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)
        
    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    if subplot_spec is not None:
        GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    grid = {}
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    plt.rcParams["font.family"] = "Arial"

    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    grid["ax_heatmap"].set_title(title)

    if strength_label is not None:
        ax.text(
            0.5, 0.5, strength_label,
            transform=ax.transAxes,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=16,
            color='black'
        )

    ax.text(0.2, -0.07, "inactive (B)", transform=ax.transAxes, ha='center', va='center', fontsize=15, color='black')
    ax.text(0.8, -0.07, "active (A)", transform=ax.transAxes, ha='center', va='center', fontsize=15, color='black')

    ax.text(-0.07, 0.2, "active (A)", rotation='vertical', transform=ax.transAxes, ha='center', va='center', fontsize=15, color='black')
    ax.text(-0.07, 0.8, "inactive (B)", rotation='vertical', transform=ax.transAxes, ha='center', va='center', fontsize=15, color='black')

    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})

    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        print('cbar')

        cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid['ax_heatmap'].grid(False)
    grid["ax_heatmap"].xaxis.set_visible(False)
    
    if title is not None:
        grid["ax_heatmap"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    plt.savefig(f'Saddle_plot_heatmap_{sample_name}.pdf',dpi=600)
    plt.savefig(f'Saddle_plot_heatmap_{sample_name}.png',dpi=600)

n=50
k=10
S=interaction_sum
C=interaction_count

intra_sum = np.nansum(S[0:k, 0:k]) + np.nansum(S[n - k : n, n - k : n])
intra_count = np.nansum(C[0:k, 0:k]) + np.nansum(C[n - k : n, n - k : n])
intra = intra_sum / intra_count

inter_sum = np.nansum(S[0:k, n - k : n]) + np.nansum(S[n - k : n, 0:k])
inter_count = np.nansum(C[0:k, n - k : n]) + np.nansum(C[n - k : n, 0:k])
inter = inter_sum / inter_count

strength_label=round(intra / inter,1)

saddleplot(eigenvector_track,
           interaction_sum/interaction_count,
           N_GROUPS,vmin=0.5, vmax=2,
           qrange=(Q_LO,Q_HI),title=sample_name,
           cbar_kws={'label':'Observed / Expected'},strength_label=strength_label
          );

print(sample_name+" finished!")

