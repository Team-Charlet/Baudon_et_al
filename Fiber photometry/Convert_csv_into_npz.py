import numpy as np
import pandas as pd

##
#import data
Import_data_list = ("20240507_RM8_Ret28_taillift_3", "20240506_RM8_Ret28_3")

for Import_data_name in Import_data_list:
    #convert csv into npz
    df = pd.read_csv(f'{Import_data_name}.csv', delimiter=',', header=1, usecols=['Time(s)', 'AIn-1 - Dem (AOut-1)', 'AIn-2 - Dem (AOut-1)', 'AIn-3 - Dem (AOut-2)', 'AIn-4 - Dem (AOut-2)'])
    df= df.T
    df = df.to_numpy()
    # print(df)

    #separate column into individual arrays
    (T, C1, C2, C3, C4) = (df[0], df[1], df[2], df[3], df[4])

    np.savez(f'{Import_data_name}.npz', Time=T, C1=C1, C2=C2, C3=C3, C4=C4)



