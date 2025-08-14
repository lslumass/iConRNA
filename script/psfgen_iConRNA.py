import sys
from psfgen import PsfGen
from HyresBuilder import utils
import argparse
import os

# Global varibale
RNA_topology = "./top_iConRNA.inp"

def main():
    parser = argparse.ArgumentParser(description="generate PSF for Hyres systems",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_pdb_files",  help="Hyres PDB file(s), it should be the pdb of monomer", required=True, nargs="+")
    parser.add_argument("-o", "--output_psf_file", help="output name/path for Hyres PSF", required=True)
    parser.add_argument("-n", "--num_of_chains", help="Number of copies for each pdb; it should have the same length as the given pdb list specified in the '-i' argument", default=[1,], nargs="+")
    parser.add_argument("-t", "--ter", choices=['neutral', 'charged'], help="Terminal charged status (choose from ['neutral', 'charged'])", default='neutral')
    args = parser.parse_args()
   
    pdb_list = args.input_pdb_files
    outpsf = args.output_psf_file
    num_list = [1,]*len(pdb_list) if len(pdb_list) > 1 and args.num_of_chains == [1,] else [int(i) for i in args.num_of_chains] 
    assert len(pdb_list) == len(num_list), "pdb file list must have the same length as the chain number list (specified by the '-n' argument)"
    ter = args.ter

    gen = PsfGen()
    gen.read_topology(RNA_topology)

    for idx, (pdb, num) in enumerate(zip(pdb_list, num_list), 1):
    # loop through each pdb and make copies
        if pdb in ['mg.pdb', 'Mg.pdb']:
            seg = 'Mg'
        else:
            seg = f"R{chr(64+idx)}"
        for i in range(num):
            segid = f"{seg}{i}" 
            gen.add_segment(segid=segid, pdbfile=pdb, auto_angles=False, auto_dihedrals=False)
        gen.write_psf(filename=outpsf)

if __name__ == '__main__':
    main()
