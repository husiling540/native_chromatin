#!/usr/bin/env python
# coding: utf-8

import warnings

warnings.filterwarnings("ignore")
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib import colors

plt.style.use('seaborn-poster')
import numpy as np
import pandas as pd
import bioframe
import cooler
import cooltools

# Sample 1
resolution = 10000
clr = cooler.Cooler('/mnt/disk3/K562_MicroC_allvalidpairs/cooler_res10000/K562_MicroC_DMSO.cool')
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

cvd_smooth_agg = cooltools.expected_cis(
    clr=clr,
    view_df=hg38_arms,
    smooth=True,
    aggregate_smoothed=True,
    nproc=80
)

cvd_smooth_agg['s_bp'] = cvd_smooth_agg['dist'] * resolution
cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan

cvd_merged_dmso = cvd_smooth_agg.drop_duplicates(subset=['dist'])[['s_bp', 'balanced.avg.smoothed.agg']]
der_dmso = np.gradient(np.log(cvd_merged_dmso['balanced.avg.smoothed.agg']), np.log(cvd_merged_dmso['s_bp']))

# Sample 2
resolution = 10000
clr = cooler.Cooler('/mnt/disk3/K562_MicroC_allvalidpairs/cooler_res10000/K562_MicroC_TPA.cool')
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

cvd_smooth_agg = cooltools.expected_cis(
    clr=clr,
    view_df=hg38_arms,
    smooth=True,
    aggregate_smoothed=True,
    nproc=40
)

cvd_smooth_agg['s_bp'] = cvd_smooth_agg['dist'] * resolution
cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan

cvd_merged_tpa = cvd_smooth_agg.drop_duplicates(subset=['dist'])[['s_bp', 'balanced.avg.smoothed.agg']]
der_tpa = np.gradient(np.log(cvd_merged_tpa['balanced.avg.smoothed.agg']), np.log(cvd_merged_tpa['s_bp']))

# Plotting
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


# Define plotting function
def plot_dual_graph(cvd_merged, axs, der, label, color):
    ax = axs[0]
    ax.loglog(
        cvd_merged['s_bp'],
        cvd_merged['balanced.avg.smoothed.agg'],
        '-',
        markersize=5,
        label=label,
        color=color,
        alpha=1
    )
    ax.legend(loc="upper right", frameon=False, prop={'size': 15})
    ax.set(
        # xlabel='Genomic distance (Mb)',
        ylabel='Contact probability',
        xlim=(1e4, 2e8),
        ylim=(10e-7, 0.1)
    )
    ax.set_aspect(0.75)

    ax = axs[1]
    ax.semilogx(
        cvd_merged['s_bp'],
        der,
        color=color,
        alpha=1
    )

    ax.set(
        xlabel='Genomic distance (Mb)',
        ylabel='slope'
    )


# Create a 2x1 subplot
fig, axs = plt.subplots(
    figsize=(6, 7),
    nrows=2,
    gridspec_kw={'height_ratios': [2, 1.2]},
    sharex=True
)

# Plot DMSO curve
plot_dual_graph(cvd_merged_dmso, axs, der_dmso, 'DMSO', "#264992")

# Plot TPA curve
plot_dual_graph(cvd_merged_tpa, axs, der_tpa, 'TPA', "#C2212B")

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.001)

# Set aspect ratio of second subplot
axs[1].set_aspect(1)

# Customize tick labels
axs[0].xaxis.set_major_locator(plt.FixedLocator([1e5, 1e6, 1e7, 1e8]))
axs[0].xaxis.set_major_formatter(
    ticker.FuncFormatter(lambda x, _: f'{x / 1e6:.1f}' if x % 1e6 != 0 else f'{int(x / 1e6):d}'))
axs[0].xaxis.set_tick_params(labelsize=14)

axs[1].xaxis.set_major_locator(plt.FixedLocator([1e5, 1e6, 1e7, 1e8]))
axs[1].xaxis.set_major_formatter(
    ticker.FuncFormatter(lambda x, _: f'{x / 1e6:.1f}' if x % 1e6 != 0 else f'{int(x / 1e6):d}'))
axs[1].xaxis.set_tick_params(labelsize=14)

fig.set_size_inches(4, 6)
plt.rcParams.update({'font.size': 14, 'font.family': 'Arial'})

# Save the figure
plt.savefig("figS3.cooler_IC_contact_frequency_DMSO_TPA.pdf", format="pdf", bbox_inches='tight')
