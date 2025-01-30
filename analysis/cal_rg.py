"""
calculate the radius of gyration of one single RNA chain
usage: python cal_rg.py psf_file trajectory_file
output: data file name as "rg.dat"
"""

import numpy as np
import MDAnalysis as mda
from sys import argv

psf = argv[1]
dcd = argv[2]

u = mda.Universe(psf, dcd)
sel = u.select_atoms('not name MG')
rgs = []
times = []
for ts in u.trajectory:
    rg = sel.radius_of_gyration()
    rgs.append(rg)
    times.append(ts.time)

with open('rg.dat', 'w') as f:
    for i in range(len(rgs)):
        print(times[i], rgs[i], file=f)

print('Rg: ', np.mean(rgs), 'err: ', np.std(rgs))
