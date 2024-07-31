# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 17:17:02 2022

@author: Angel.BAUDON
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt, glob
from mlxtend.plotting import plot_pca_correlation_graph
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


folder = r"C:\Angel.BAUDON\Biblio\Publi\MoprhAstro in prep\v1 - Copy\Data\Figure 5\EPhy_Puff KCL biased agonists\PCA"
df = pd.read_excel(glob.glob(rf'{folder}\*.xlsx')[0])
show_fig = False

# =====================================
# K-Means Clustering
# =====================================
kmeans = KMeans(n_clusters=2)
df['Resp'] = kmeans.fit_predict(df[['ΔVm']])
centroids = kmeans.cluster_centers_


# =====================================
# PCA
# =====================================
features = ['Threshold', '0 = ES / 1 = LS', 'NA = 0 / A = 1', 'Cm', 'Rm',
            'Rise', 'Decay', 'Vm', 'ΔVm', 'Δspike', 'Left = 0 / Right = 1']
target, col = 'Resp', ('PC 1', 'PC 2')
x = StandardScaler().fit_transform(df.loc[:, features].values)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns=col)

Df = pd.concat([principalDf, df[[target]]], axis = 1)
per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)

plt.figure(), plt.title('2 component PCA'), plt.grid()
plt.xlabel(f'PC1 - {per_var[0]}%'), plt.ylabel(f'PC2 - {per_var[1]}%')
for target_label, color in zip((1, 0), ('r', 'g')):
    ToKeep = Df[target] == target_label
    plt.scatter(Df.loc[ToKeep, 'PC 1'], Df.loc[ToKeep, 'PC 2'], c=color)
plt.legend(('Responsive', 'Non responsive')), plt.tight_layout()
plt.savefig(rf'{folder}\PCA.pdf')
if not show_fig: plt.close()
print('\n\n\nExplained Variance = ', per_var, '\n\n\n')


# =====================================
# PCA Correlations Plot
# =====================================
figure, correlation_matrix = plot_pca_correlation_graph(x, features)
plt.tight_layout()
plt.savefig(rf'{folder}\PCA Correlation plot.pdf')
if not show_fig: plt.close()

# # =====================================
# # Random Forest Classifier
# # =====================================
# # https://www.lovelyanalytics.com/2020/06/08/random-forest-tutoriel-python/
# from sklearn.model_selection import train_test_split
# from sklearn.ensemble import RandomForestRegressor

# df_rf = df.copy(deep=False)
# pd.DataFrame.set_index(df_rf, 'Name', inplace=True)
# cible = np.array(df_rf['ΔVm'])
# df_rf = df_rf.drop('Resp', axis = 1)
# df_rf = df_rf.drop('ΔVm', axis = 1)
# x_train, x_test, y_train, y_test = train_test_split(df_rf, cible, test_size=0.25)


# # =================================================
# #              Without R & D              
# # =================================================

# x_train, x_test = x_train.drop('Rise', axis = 1), x_test.drop('Rise', axis = 1)
# x_train, x_test = x_train.drop('Decay', axis = 1), x_test.drop('Decay', axis = 1)
# x_train, x_test = x_train.drop('Δspike', axis = 1), x_test.drop('Δspike', axis = 1)

# liste_variables = list(x_train.columns)

# stock_mape, stock_mae, stock_features = [], [], []
# for i in range(10):
  
#     #Train 
#     rf = RandomForestRegressor(n_estimators = 10000)
#     rf.fit(x_train, y_train)
    
#     #Test
#     predictions = rf.predict(x_test)
    
#     # MAE : Mean Asolute Error
#     erreurs = abs(predictions - y_test)
#     mae = round(np.mean(erreurs), 2)

#     # MAPE : Mean Absolute Percentage Error
#     mape = round(np.mean(100*erreurs/y_test), 2)
    
#     importances = rf.feature_importances_
#     # indices = np.argsort(importances)

#     stock_mape.append(mape), stock_mae.append(mae), stock_features.append(importances)
    


# x = np.asarray(stock_features)
# means, sems = np.nanmean(x, axis=0), stat.sem(x, axis=0)
# y = np.arange(len(liste_variables))
# indices = np.argsort(means)

# means, sems = means[indices], sems[indices]
# yticks = [liste_variables[i] for i in indices]

# plt.figure(), plt.title('Feature Importances')
# plt.barh(y, means, color='b', align='center'), plt.yticks(y, yticks)
# plt.errorbar(means, y, xerr=sems, c='k', ls='', capsize=1)
# plt.xlabel('Relative Importance'), plt.tight_layout(), plt.grid()
# plt.savefig(rf'{folder}\Features dvm without R&D.pdf')
# plt.close()

# writer = pd.ExcelWriter(f'{folder}/Random Forest.xlsx')
# pd.DataFrame(stock_mae).to_excel(writer, sheet_name = 'MAE')
# pd.DataFrame(stock_mape).to_excel(writer, sheet_name = 'MAPE')
# pd.DataFrame(stock_features).to_excel(writer, sheet_name = 'Features')
# writer.save()

# print(np.nanmean(stock_mape))

