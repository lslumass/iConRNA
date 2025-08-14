import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import polymer
from sys import argv
import matplotlib.pyplot as plt
import mdtraj as md

psf = 'conf.psf'
dcd = 'system.xtc'

u = mda.Universe(psf, dcd)
chains = u.atoms.fragments
bbs = [ch.select_atoms('name P') for ch in chains]

plen = polymer.PersistenceLength(bbs)
plen.run()

print('lp', plen.lp, 'A')

plen.plot()
plt.show()
