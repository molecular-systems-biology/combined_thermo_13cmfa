'''script to compute the shape similarity between two distributions. 
With 2 distributions, get both normalized distribs and align them on the same mean. then compute overlap area
Overlap area :area under the curve of the minimum of the two distributions.

Input : net flux samples dataframe (all samples : first distribution)
list of accepted samples (considering only those : second distribution)

Output :dict of overlap areas for each reaction'''
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde,wasserstein_distance
from scipy.integrate import simpson

class DistribComparer:
    def __init__(self, dist1, dist2,  num_points=10000):
        """
        net_flux_df: DataFrame with reactions as index and samples as columns.
        accepted_samples: list of column names to define second distribution.
        dist2: Optional DataFrame of precomputed accepted flux distributions.
        num_points: Resolution of KDE evaluation.
        """

        self.dist1 = dist1
        self.dist2 = dist2
        self.index= dist1.index.intersection(dist2.index)
        self.num_points = num_points
        # self.align=align
    
    @staticmethod
    def _align_means(dist1, dist2):
        mean1 = np.mean(dist1)
        mean2 = np.mean(dist2)
        return dist2 + (mean1 - mean2)
    
    @staticmethod
    def _get_kde(dist):
        return gaussian_kde(dist)

    def _compute_overlap_area(self, kde1, kde2):
        x_min = min(kde1.dataset.min(), kde2.dataset.min()) - 1
        x_max = max(kde1.dataset.max(), kde2.dataset.max()) + 1
        x = np.linspace(x_min, x_max, self.num_points)
        y1 = kde1(x)
        y2 = kde2(x)
        min_y = np.minimum(y1, y2)

        return simpson(min_y, x=x)

    def compute_overlap_one(self, idx,align=True):
        try:
            pre_dist1= self.dist1.loc[idx] #should be calle dist1_rxn_i or something
            pre_dist2 = self.dist2.loc[idx]
            dist2_aligned = self._align_means(pre_dist1, pre_dist2) if align else pre_dist2
            kde1 = self._get_kde(pre_dist1) ## could mutualize this with a dict of kdes1 and kdes2 by rxn
            kde2 = self._get_kde(dist2_aligned)
            return self._compute_overlap_area(kde1, kde2)
        except Exception as e:
            print(f"null area {idx}",e)
            return None

    def compute_overlap_all(self,align=True):
        overlaps = {}
        for idx in self.index:
            area = self.compute_overlap_one(idx,align=align)
            if area is not None :
                overlaps[idx] = area if area >= 1e-4 else 0
        return overlaps
    
    ###EMD instead of overlap 

    def compute_EMD_one(self,idx):
        # pre_dist1= self.dist1.loc[idx] 
        # pre_dist2 = self.dist2.loc[idx]
        # return wasserstein_distance(pre_dist1,pre_dist2)
        pre_dist1 = self.dist1.loc[idx].values
        pre_dist2 = self.dist2.loc[idx].values

        # joint scaling range
        xmin = min(pre_dist1.min(), pre_dist2.min())
        xmax = max(pre_dist1.max(), pre_dist2.max())

        # create normalized copies
        if xmax > xmin:
            norm_dist1 = 2 * (pre_dist1 - xmin) / (xmax - xmin) - 1
            norm_dist2 = 2 * (pre_dist2 - xmin) / (xmax - xmin) - 1
        else:
            norm_dist1 = pre_dist1.copy()
            norm_dist2 = pre_dist2.copy()

        return wasserstein_distance(norm_dist1, norm_dist2)


    def compute_EMD_all(self):
        overlaps = {}
        for idx in self.index:
            area = self.compute_EMD_one(idx)
            if area is not None :
                overlaps[idx] = area if area >= 1e-4 else 0
        return overlaps

    def plot_reaction(self, ax, reaction,align=True,colors=None):
        pre_dist1 = self.dist1.loc[reaction]
        pre_dist2 = self.dist2.loc[reaction]


        dist2_aligned = self._align_means(pre_dist1, pre_dist2) if align else pre_dist2

        kde1 = self._get_kde(pre_dist1)
        kde2 = self._get_kde(dist2_aligned)

        x_min = min(pre_dist1.values.min(), dist2_aligned.min()) - 1
        x_max = max(pre_dist1.values.max(), dist2_aligned.max()) + 1
        x = np.linspace(x_min, x_max, self.num_points)

        y1 = kde1(x)
        y2 = kde2(x)
        min_y = np.minimum(y1, y2)
        overlap_area = simpson(min_y, x=x)

        ##cond spec colors can be given as a list of 2 colors
        c1 = None if colors is None else colors[0]
        c2 = None if colors is None else colors[1]
        #plot
        ax.plot(x, y1, label='All samples', lw=2,color=c1)
        ax.plot(x, y2, label='Best scoring samples', lw=2,color=c2)

        ax.axvline(np.mean(pre_dist1), color='black', linestyle='--', lw=2, label='Mean (all samples)')
        ax.fill_between(x, min_y, alpha=0.4, color='gray', label=f'Overlap area = {overlap_area:.3f}')
        return ax

    def get_deltamean(self):
        ### get delta_mean
        delta_means = abs(self.dist2.mean(axis=1) - self.dist1.mean(axis=1))
        self.deltameans= pd.Series(delta_means, index=self.index)
        return self.deltameans  
