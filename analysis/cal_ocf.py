import numpy as np
import MDAnalysis as mda
from sys import argv

psf = argv[1]

u = mda.Universe(psf, 'system.dcd')
sel = u.select_atoms('name P')
n_atom = sel.n_atoms
ocf1 = [0]*(n_atom-1)
nums = [0]*(n_atom-1)
for ts in u.trajectory[1000:]:
    rxys = []
    bds = []
    for atom in sel.atoms:
        rxys.append(np.array(atom.position))
    for i in range(n_atom-1):
        bd = rxys[i+1] - rxys[i]
        bds.append(bd/np.linalg.norm(bd))
    for i in range(len(bds)-1):
        for j in range(i+1, len(bds)):
            ij = j - i
            ocf = np.dot(bds[i], bds[j])
            ocf1[ij] += ocf
            nums[ij] += 1
for i in range(1, n_atom-1):
    ocf1[i] = ocf1[i]/nums[i]

with open('ocf.dat', 'w') as f:
    for i in range(1, n_atom-1):
        print(i, ocf1[i], abs(ocf1[i]/10), file=f)

print('Done!')
