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

psf = argv[1]
dcd = argv[2]

u = mda.Universe(psf, dcd)
pairs = HydrogenBondAnalysis(
    universe=u,
    donors_sel="resname GUA ADE and name NB",
    hydrogens_sel="resname GUA ADE and name NC",
    acceptors_sel="resname CYT URA and name NB",
    d_a_cutoff=10.0,
    d_h_cutoff=3.5,
    d_h_a_angle_cutoff=150, 
    update_selections=False
)

pairs.run(step=1, verbose=True)

counts = pairs.count_by_time()
with open('pair_number.dat', 'w') as f:
    for frame, num in zip(pairs.frames, counts):
        print((frame+1)*u.trajectory.dt, num, file=f)

print('Finished!')

 
