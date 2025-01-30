"""
calculate the end-to-end distance of one single RNA chain
usage: python cal_re.py psf_file trajectory_file
output: data file name as "re.dat"
"""


import numpy as np
from sys import argv
import mdtraj

psf = argv[1]
dcd = argv[2]

traj = mdtraj.load(dcd, top=psf)
top = traj.topology
start = top.select('resid 0 and name P')    # change the resid as the first residue 
end = top.select('resid 32 and name C1')    # change the resid as the last residue
res = mdtraj.compute_distances(traj, [[start[0], end[0]],])
all_re = []
with open('re.dat', 'w') as f:
    for i in range(len(res)):
        print(i, res[i][0], file=f)
        all_re.append(res[i][0])

print('Re (nm): ', np.mean(all_re), 'err: ', np.std(all_re))
