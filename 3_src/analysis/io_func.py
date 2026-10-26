import yaml
import seaborn as sns
import numpy as np

from cmap import Colormap

def import_plotconfigs(cpalette=False):
    '''function to import plot settings from cfg 
    return the palette objcet and cfg dict'''
    # Import general configs:
    with open("/home/nathan/PhD/04.Computational/config.yml", "r") as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        
    sns.set(**cfg["general_plotting"])

    # Create custom palettes
    paired = sns.color_palette(cfg["colors"]["palette"])
    for color in cfg["colors"]["extraColors"]:
        paired.append(eval(color))

    ##set batlopex cmap
    cm = Colormap('crameri:batlow')  
    cm_colors=cm.as_dict()['value']
    batlopex = [ cm_colors[0][1],  cm_colors[32][1], cm_colors[64][1], cm_colors[96][1], cm_colors[128][1], cm_colors[160][1],
                 cm_colors[192][1], cm_colors[224][1], 
    ]
    palette=([batlopex[2]]+ [sns.light_palette(batlopex[2],n_colors=4)[1]]  
    + [batlopex[4] ]+ [sns.light_palette(batlopex[4],n_colors=4)[1]]
    )
    extra_color_tm6=sns.light_palette(palette[2],n_colors=8)[4]
    extra_color_wt=sns.light_palette(palette[0],n_colors=8)[4]
    palette.append(extra_color_wt)
    palette.append(extra_color_tm6)
    sns.set_palette(palette)
    if cpalette:
        return palette,cm.to_matplotlib(),cfg
    else :
        return palette,cfg

def scale_rcparams(fig,scalingf=None):
    default_figsize = (5, 5)
    current_figsize = fig.get_size_inches()
    scale_x = current_figsize[0] / default_figsize[0]
    scale_y = current_figsize[1] / default_figsize[1]
    scale = min(scale_x, scale_y) if scalingf is None else scalingf
    
    # Update rcParams with scaled values
    plt.rcParams['axes.linewidth'] = 2.05 * scale
    plt.rcParams['lines.linewidth'] = 1.5 * scale
    plt.rcParams['font.size'] = 12.5# 25 * scale #fixed size of 12.5 for all
    plt.rcParams['xtick.major.size'] = 9.84 * scale
    plt.rcParams['ytick.major.size'] = 9.84 * scale
    plt.rcParams['xtick.minor.size'] = 6.56 * scale
    plt.rcParams['ytick.minor.size'] = 6.56 * scale
    #width
    # 'ytick.minor.width': 1.0
    plt.rcParams['xtick.major.width'] = 2.05 * scale
    plt.rcParams['ytick.major.width'] = 2.05 * scale
    plt.rcParams['xtick.minor.width'] = 1.64 * scale
    plt.rcParams['ytick.minor.width'] = 1.64 * scale
    # plt.rcParams['lines.linewidth'] = 1.5 * scale
    plt.rcParams['legend.fontsize'] = 12.5#5 * scale
    plt.rcParams['xtick.labelsize'] =12.5 #25 * scale
    plt.rcParams['ytick.labelsize'] = 12.5#25 * scale
    plt.rcParams['axes.labelsize'] =12.5# 25 * scale
    
    return plt.rcParams

def rescale_axlim(ax,how='both',nbins=4):
    '''## just great on linear scale '''
    def rescale_ylim(ax,nbins):
        yticks = ax.get_yticks()
        ax.set_yticks(np.linspace(yticks[0], yticks[-1], nbins))
        return ax
    def rescale_xlim(ax,nbins):
        xticks = ax.get_xticks()
        ax.set_xticks(np.linspace(xticks[0], xticks[-1], nbins))
        return ax
    if how=='both':
        ax=rescale_ylim(rescale_xlim(ax,nbins),nbins)
    elif how=='x':
        ax=rescale_xlim(ax,nbins)
    elif how=='y':
        ax=rescale_ylim(ax,nbins)
    return ax

def plot_diag(ax,cfg):
    ax.plot([0,1],[0,1],transform=ax.transAxes,**cfg['settings_diagonal'])
    return ax


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator
### from stackoverflow
class MinorSymLogLocator(Locator):
    ###to use : ax.xaxis.set_minor_locator(io_func.MinorSymLogLocator(1e-1)) and can use ax.set_xticks([-10,-0.1,0,0.1,10])
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()
        # add temporary major tick locs at either end of the current range
        # to fill in minor tick gaps
        dmlower = majorlocs[1] - majorlocs[0]    # major tick difference at lower end
        dmupper = majorlocs[-1] - majorlocs[-2]  # major tick difference at upper end

        # add temporary major tick location at the lower end
        if majorlocs[0] != 0. and ((majorlocs[0] != self.linthresh and dmlower > self.linthresh) or (dmlower == self.linthresh and majorlocs[0] < 0)):
            majorlocs = np.insert(majorlocs, 0, majorlocs[0]*10.)
        else:
            majorlocs = np.insert(majorlocs, 0, majorlocs[0]-self.linthresh)

        # add temporary major tick location at the upper end
        if majorlocs[-1] != 0. and ((np.abs(majorlocs[-1]) != self.linthresh and dmupper > self.linthresh) or (dmupper == self.linthresh and majorlocs[-1] > 0)):
            majorlocs = np.append(majorlocs, majorlocs[-1]*10.)
        else:
            majorlocs = np.append(majorlocs, majorlocs[-1]+self.linthresh)
            
        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))
    

def symlog_distance(wt, tm6, linthresh):
    ##for global flux plot WT vs TM6
    """
    Compute the perpendicular distance to the y=x line in symlog space.
    
    Parameters:
    - wt: WT values
    - tm6: TM6 values
    - linthresh: Linear threshold for symlog scaling
    
    Returns:
    - Distances in symlog space
    """
    # Transform WT and TM6 to symlog space
    wt_symlog = np.sign(wt) * np.log1p(np.abs(wt / linthresh))
    tm6_symlog = np.sign(tm6) * np.log1p(np.abs(tm6 / linthresh))
    
    # Compute the distance in symlog space
    distance = np.abs(wt_symlog - tm6_symlog) / np.sqrt(2)
    return distance

def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)