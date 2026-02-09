import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
import seaborn as sns
import palettable


def loadtxtmethod(filename):
    data = np.loadtxt(filename,dtype=np.float32,delimiter=',')
    return data

if __name__=="_main_":
    data=loadtxtmethod()
    print(data)


def read_tablemethod(filename):
    data = pd.read_table(filename, header=None, delim_whitespace=True)
    return data

if __name__ == "__main__":
    data = read_tablemethod()
    print(data)

NT_RNAmatrix_rho
print(NT_RNAmatrix_rho)

seaborn.clustermap(NT_RNAmatrix_rho, pivot_kws=None, method='average', metric='euclidean', z_score=None, standard_scale=None,
                   figsize=(10, 10), cbar_kws=None, row_cluster=True, col_cluster=True, row_linkage=None, col_linkage=None,
                   row_colors=None, col_colors=None, mask=None, dendrogram_ratio=0.2, colors_ratio=0.03,
                   cbar_pos=(0.02, 0.8, 0.05, 0.18), tree_kws=None, **kwargs)


row_c = dict(zip(NT_05_RNAmatrix_top50_rho['class'].unique(), 'terrain'))

NT=sns.clustermap(NT_RNAmatrix_rho, row_cluster=True, col_cluster=True, method='ward', metric='euclidean',
                        #row_colors=pd_iris['class'].map(row_c),
                        cmap="seismic", vmax=0.85, vmin=-0.85,
                        dendrogram_ratio=(.1, .1),
                        cbar = False,
                        )

# 添加标题, fontweight='bold'
#plt.title('AAA', fontsize=35) # 'WS(m/s)'  u'T(°C)'  'P(mbar)'
#plt.xticks(fontsize=28)
#plt.yticks(fontsize=28)
#plt.xlabel('longitude', fontsize=34)  # 经度
#plt.ylabel('latitude', fontsize=34)   # 纬度
#plt.rcParams['savefig.dpi'] = 600
#plt.tight_layout()
#plt.savefig()

# 显示图形
#plt.show()


NT_ind = NT.dendrogram_col.reordered_ind
file = open('NT_ind.txt','w')
file.write(str(NT_ind))
file.close()



