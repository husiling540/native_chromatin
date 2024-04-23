#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cooler
import bioframe
import cooltools
import cooltools.lib.plotting
from packaging import version
from matplotlib.font_manager import FontProperties
import matplotlib as mpl

if version.parse(cooltools.__version__) < version.parse('0.5.2'):
    raise AssertionError("tutorials rely on cooltools version 0.5.2 or higher," +
                         "please check your cooltools version and update to the latest")

data_dir = 'K562_MicroC_allvalidpairs/cooler_res5000/'
# Open cool file with Micro-C data:
clr_dmso = cooler.Cooler(data_dir + '/K562_MicroC_DMSO.cool')
clr_tpa = cooler.Cooler(data_dir + '/K562_MicroC_TPA.cool')
# Set up selected data resolution:
resolution = clr_dmso.binsize

# Use bioframe to fetch the genomic features from the UCSC.
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
hg38_arms = hg38_arms.set_index("chrom").loc[clr_dmso.chromnames].reset_index()

loop_prefix_list = ["S1", "S2", "S3", "S4", "S5", "S6"]

for loop_prefix in loop_prefix_list:
    file_name = f'K562_MicroC_mustache_DMSO_loop.{loop_prefix}.bedpe'
    column_names = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', "strand"]
    interaction_file = bioframe.read_table(file_name, names=column_names)
    paired_sites = interaction_file.query(f'chrom1 in {clr_dmso.chromnames}')
    columns_to_convert = ['start1', 'end1', 'start2', 'end2']

    for col in columns_to_convert:
        paired_sites[col] = paired_sites[col].astype(int)

    paired_sites.loc[:, 'mid1'] = (paired_sites['start1'] + paired_sites['end1']) // 2
    paired_sites.loc[:, 'mid2'] = (paired_sites['start2'] + paired_sites['end2']) // 2

    globals()[f'loops_{loop_prefix}'] = paired_sites

expected_dmso = cooltools.expected_cis(clr_dmso, view_df=hg38_arms, nproc=80, chunksize=10_000)
expected_tpa = cooltools.expected_cis(clr_tpa, view_df=hg38_arms, nproc=80, chunksize=10_000)

loop_prefix_list = ["S1", "S2", "S3", "S4", "S5", "S6"]

for loop_prefix in loop_prefix_list:
    for condition in ["dmso", "tpa"]:
        stack = cooltools.pileup(eval(f'clr_{condition}'), eval(f'loops_{loop_prefix}'), view_df=hg38_arms,
                                 expected_df=eval(f'expected_{condition}'), flank=100_000)
        sites = eval(f'loops_{loop_prefix}')
        mask = np.array(sites.strand == '-', dtype=bool)
        stack[:, :, mask] = stack[::-1, ::-1, mask]

        globals()[f'mtx_{condition}_{loop_prefix}'] = np.nanmean(stack, axis=2)

        stack_shape = stack.shape
        num_layers = stack_shape[2]
        averages = []

        for layer in range(num_layers):
            current_layer = stack[:, :, layer]
            center_region = current_layer[19:22, 19:22]
            center_region[np.isnan(center_region)] = 0
            average = np.nanmean(center_region)
            averages.append(average)

        averages = np.array(averages)
        np.savetxt(f'fig7.mtx_{condition}_{loop_prefix}_center_signal_every_loop.tsv', averages, delimiter='\t')

        print(f'mtx_{condition}_{loop_prefix} finish!')

for loop_prefix in loop_prefix_list:
    for condition in ["dmso", "tpa"]:
        globals()[f'mtx_{condition}_{loop_prefix}_p'] = np.where(eval(f'mtx_{condition}_{loop_prefix}') == 0, np.nan,
                                                                 eval(f'mtx_{condition}_{loop_prefix}'))
        np.savetxt(f'fig7.APA_plot_mtx_{condition}_{loop_prefix}.tsv', eval(f'mtx_{condition}_{loop_prefix}_p'),
                   delimiter='\t')

    globals()[f'mtx_fc_{loop_prefix}_p'] = eval(f'mtx_tpa_{loop_prefix}_p') / eval(f'mtx_dmso_{loop_prefix}_p')
    np.savetxt(f'fig7.APA_plot_mtx_fc_{loop_prefix}.tsv', eval(f'mtx_fc_{loop_prefix}_p'), delimiter='\t')


def get_enrichment(amap, n):
    """Get values from the center of a pileup for a square with side *n*

    Parameters
    ----------
    amap : 2D array
        Pileup.
    n : int
        Side of the central square to use.

    Returns
    -------
    enrichment : float
        Mean of the pixels in the central square.

    """
    c = amap.shape[0] // 2
    if c < n:
        raise ValueError(f"Central pixel value {n} is too large, can be maximum {c}")
    return np.nanmean(amap[c - n // 2: c + n // 2 + 1, c - n // 2: c + n // 2 + 1])


mpl.rcParams['pdf.fonttype'] = 42

font_properties = FontProperties(family='Arial')
conditions = ['dmso', 'tpa', "fc"]
loop_prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]

flank = 100_000

cmap = plt.cm.get_cmap('coolwarm')
cmap.set_bad(color='#E6E2DF')

fig, axes = plt.subplots(nrows=len(conditions), ncols=len(loop_prefixes), figsize=(4, 4), sharex=True, sharey=True)

for i, condition in enumerate(conditions):
    for j, loop_prefix in enumerate(loop_prefixes):
        mtx_p = eval(f'mtx_{condition}_{loop_prefix}_p')
        ax = axes[i, j]
        score = round(get_enrichment(mtx_p, 3), 2)
        ax.axis('off')
        if condition == "fc":
            im = ax.imshow(mtx_p, vmax=1.6, vmin=0.4, cmap=cmap, interpolation='none')
        else:
            im = ax.imshow(mtx_p, vmax=4, vmin=0, cmap=cmap, interpolation='none')
            normal_im = im

        ticks_pixels = np.linspace(0, flank * 2 // resolution, 5)
        ticks_kbp = ((ticks_pixels - ticks_pixels[-1] / 2) * resolution // 1000).astype(int)
        ax.set_xticks(ticks_pixels)
        ax.set_yticks([])
        ax.set_xticklabels(ticks_kbp, fontsize=16)
        ax.set_yticklabels([])

        if i == len(conditions) - 1 and j == len(loop_prefixes) - 1:
            ax.set_xlabel('Relative Position (kbp)', fontsize=16)
        if j == 0:
            if condition == "fc":
                ax.text(-0.1, 0.5, "TPA / DMSO", transform=ax.transAxes, va='center', ha='right', rotation='vertical',
                        fontsize=16)
            else:
                ax.text(-0.1, 0.5, condition.upper(), transform=ax.transAxes, va='center', ha='right',
                        rotation='vertical', fontsize=16)
        if i == 0:
            ax.text(0.5, 1.05, loop_prefix, transform=ax.transAxes, va='bottom', ha='center', fontsize=16)

        ax.text(0.05, 0.9, score, transform=ax.transAxes, va='top', ha='left', fontsize=15)

plt.subplots_adjust(wspace=0.05, hspace=0.05)
cax = plt.axes([0.92, 0.4, 0.02, 0.45])
colorbar = plt.colorbar(normal_im, cax=cax, ticks=[0, 2, 4], label='Observed / Expected')
cax2 = plt.axes([0.92, 0.12, 0.02, 0.22])
colorbar2 = plt.colorbar(im, cax=cax2, ticks=[0.4, 1, 1.6], label='TPA / DMSO')
fig.set_size_inches(10, 5)
plt.rcParams.update({'font.size': 16, 'font.family': 'Arial'})
plt.savefig('fig7.APA_plot.pdf', format="pdf", bbox_inches='tight')