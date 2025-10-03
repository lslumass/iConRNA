#!/usr/bin/python

"""
this script is used for the analysis of general hydrogen bond, meaning donor-acceptor system, based on MDAnalysis
authour: Shanlong Li@UMass
date: 09-04-2024
"""

from __future__ import division
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from sys import argv

dcd = argv[1]
out = argv[2]

u = mda.Universe('conf.psf', dcd)
nframes = u.trajectory.n_frames
rna = u.select_atoms('all')
dim = np.array([100, 100, 100, 90, 90, 90])
#for ts in u.trajectory:
#    ts.dimensions=dim
#    center = rna.center_of_geometry()
#    box_center = np.sum(dim, axis=0)/2
#    u.atoms.translate(box_center - center)


CGs = HydrogenBondAnalysis(
    universe=u,
    donors_sel="resname GUA and name NB",
    hydrogens_sel="resname GUA and name NC",
    acceptors_sel="resname CYT and name NB",
    d_a_cutoff=10.0,
    d_h_cutoff=3.5,
    d_h_a_angle_cutoff=150, 
    update_selections=False
)
CGs.run(step=1, verbose=True)

## remove the unreasonable hbonds
CGs_count = [[] for _ in range(nframes)]
for hbond in CGs.results.hbonds:
    frame, donor = hbond[:2].astype(int)
    CGs_count[frame].append(donor)

for i in range(nframes):
    CGs_count[i] = list(set(CGs_count[i]))

AUs = HydrogenBondAnalysis(
    universe=u,
    donors_sel="resname ADE and name NB",
    hydrogens_sel="resname ADE and name NC",
    acceptors_sel="resname URA and name NB",
    d_a_cutoff=10.0,
    d_h_cutoff=3.5,
    d_h_a_angle_cutoff=150, 
    update_selections=False
)
AUs.run(step=1, verbose=True)

## remove the unreasonable hbonds
AUs_count = [[] for _ in range(nframes)]
for hbond in AUs.results.hbonds:
    frame, donor = hbond[:2].astype(int)
    AUs_count[frame].append(donor)

for i in range(nframes):
    AUs_count[i] = list(set(AUs_count[i]))

with open(out, 'w') as f:
    for frame, num1, num2 in zip(CGs.frames, CGs_count, AUs_count):
        print(frame+1, len(num1)+len(num2), file=f)

print('Finished!')

 
