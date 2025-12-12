import os
import numpy as np, flexpart as fp

folder = r'C:\Projects\bioice\data\FLEXWRF'
fnamelist = os.listdir(folder)
print(np.sort(fnamelist)[0])

'''
partout columns are structured the following way:
(1) release; (2) timepoint; (3) longitude; (4) latitude; (5) altitude; (6) topography; (7) potential vorticity; (8) water mixing ratio; (9) air density;
(10) PBL height, m agl; (11) tropopause height, m agl; (12) temperature, K; (13) mass for each particle
'''

fname = os.path.join(folder,fnamelist[0])
partout, time, numpart = fp.read_partposition(fname)
mixingratio
