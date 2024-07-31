# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 13:00:06 2022

@author: Angel.BAUDON
"""

import pandas as pd, matplotlib.pyplot as plt, scipy, seaborn as sns, glob
from sklearn.cluster import KMeans


folder = r"C:\Angel.BAUDON\Biblio\Publi\MoprhAstro in prep\v1 - Copy\Data\Figure 5\EPhy_Puff KCL biased agonists\PCA"
df = pd.read_excel(glob.glob(rf'{folder}\*.xlsx')[0])
show_fig = False
targets = ('Threshold', 'Rm')

for target in targets:
    # =====================================
    # K-Means Clustering
    # =====================================
    kmeans = KMeans(n_clusters=2)
    df['Resp'] = kmeans.fit_predict(df[['ΔVm']])
    centroids = kmeans.cluster_centers_
    
    
    Resp = df.loc[df['Resp'] == 1][target].values
    Non_resp = df.loc[df['Resp'] == 0][target].values
    
    stat, pval = scipy.stats.ttest_ind(Resp, Non_resp)
    print(stat, pval, '\n\n')
    
    plt.figure(), plt.title(f'{target} comparison   pval = {pval}')
    sns.barplot(data=df, x='Resp', y=target)
    sns.swarmplot(data=df, x='Resp', y=target, size=10,
                  color='w', edgecolor='k', linewidth=1)
    plt.tight_layout()
    plt.savefig(rf'{folder}\{target}.pdf')
    plt.close()
