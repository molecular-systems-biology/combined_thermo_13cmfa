import sys
import os
import argparse
import warnings

sys.path.append(os.path.abspath('../3_src/'))

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib import ticker
import analysis.io_func as io_func
from analysis import config_13c as config

warnings.simplefilter(action='ignore', category=FutureWarning)


def parse_args():
    parser = argparse.ArgumentParser(description="Plot KDE distributions for flux samples.")

    parser.add_argument("condition", type=str)
    parser.add_argument("rxn_names", type=str,
                        help="Comma-separated reaction names")
    parser.add_argument("output_folder", type=str)

    parser.add_argument(
        "--plot_all",
        action="store_true",
        help="Plot KDE for all samples"
    )
    # multiple percentiles
    parser.add_argument(
        "--nperc",
        nargs="+",
        type=int,
        default=[],
        help="Additional percentiles to plot (e.g. --nperc 5 10 30)"
    )

    return parser.parse_args()


args = parse_args()

condition = args.condition
rxn_names = args.rxn_names
plot_all = args.plot_all
OUTPUT_FOLDER = args.output_folder
npercs = args.nperc

print(vars(args))
print(f'plotting distrib for {rxn_names} of {condition} strain')

# import vars
model = config.model
net_flux_df = config.full_flux_df(condition)
obj_vals = config.obj_vals(condition)
bootstrap_flux_df = config.bootstrap_flux_df(condition)
print('loaded, plotting')

palette, cfg = io_func.import_plotconfigs()
palette_cond=[palette[0],palette[4],palette[1]] if condition == 'WT' else [palette[2],palette[5],palette[3]]

def plot_onefig(rxn_name):

    if rxn_name not in net_flux_df.index:
        raise ValueError(f'{rxn_name} not known')

    fig, ax = plt.subplots(1, 1, figsize=(1, 0.5), tight_layout=False)
    io_func.scale_rcparams(fig, scalingf=0.3)

    samp_samesize = np.random.choice(net_flux_df.columns,
                                     200,
                                     replace=False)


    if plot_all:
        sns.kdeplot(
            net_flux_df.loc[rxn_name, samp_samesize],
            ax=ax,
            color='black',
            label='All samples',
            linewidth=1
        )

    # additional percentiles
    for i, perc in enumerate(npercs):
        accepted_additional = obj_vals[
            obj_vals <= np.percentile(obj_vals, perc)
        ].index

        color_idx = i

        sns.kdeplot(
            bootstrap_flux_df.loc[rxn_name],
            ax=ax,
            color=palette_cond[color_idx],
            label=f'Top {perc}% samples',
            linewidth=1
        )

    # style
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_yticks([])
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=10)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.locator_params(axis='both', nbins=3)
    ax.spines[['right', 'top', 'left']].set_visible(False)
    ax.set_title(rxn_name, fontsize=10, loc='right')
    ax.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

    OUTPUT_name = f'{OUTPUT_FOLDER}/13C_smallkde_{condition}_{rxn_name}_plotall{plot_all}'

    plt.savefig(f'{OUTPUT_name}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{OUTPUT_name}.svg', dpi=300, bbox_inches='tight')

    print(f'done for {rxn_name}, saved in {OUTPUT_FOLDER}')


if __name__ == "__main__":

    rxn_nameslist = rxn_names.split(',')

    for i, r_name in enumerate(rxn_nameslist):
        plot_onefig(r_name)

        # apply rcparams twice as in your original script
        if i == 0:
            plot_onefig(r_name)
