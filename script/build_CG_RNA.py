import sys
import numpy as np
from HyresBuilder import RNAbuilder

pdb = sys.argv[1]
seq = sys.argv[2]
RNAbuilder.build(seq, pdb)
