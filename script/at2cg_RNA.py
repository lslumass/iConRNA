import sys
import numpy as np
import MDAnalysis as mda


'''
Usage: python at2cg_RNA.py input_AA_pdb_file output_CG_pdb_file
'''

inp = sys.argv[1]   #input atomistic pdb file
out = sys.argv[2]   #output hyres pdb file
u = mda.Universe(inp)

# output in pdb-format
def printcg(atom, f):
    f.write("%4s  %5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n" % (atom[0],int(atom[1]),atom[2],atom[3][:3],atom[4], int(atom[5]),float(atom[6]),float(atom[7]),float(atom[8]),float(atom[9]),float(atom[10]), atom[11]))

def aa2cg(sel, resid, segid, cg_bead):
    atom_grp = u.select_atoms(sel)
    com = atom_grp.center_of_mass()
    atom = ['ATOM', idx, cg_bead, res.resname, 'X', resid, com[0], com[1], com[2], 1.00, 0.00, segid]
    printcg(atom, f)

idx = 1
with open(out, 'w') as f:
    for segment in u.segments:
        segid = segment.segid
        for res in segment.residues:
            resname = res.resname
            resid = res.resid
            sel = f"(name P O1P O2P O5' and resid {resid} and segid {segid}) or (name O3' and resid {resid-1} and segid {segid})"
            aa2cg(sel, resid, segid, 'P')
            idx += 1
            sel = f"name C4' and resid {resid} and segid {segid}"
            aa2cg(sel, resid, segid, 'C1')
            idx += 1
            sel = f"name C1' and resid {resid} and segid {segid}"
            aa2cg(sel, resid, segid, 'C2')
            idx += 1

            if resname == 'ADE':
                sel = f"name N9 C4 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NA')
                idx += 1
                sel = f"name C8 H8 N7 C5 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NB')
                idx += 1
                sel = f"name C6 N1 N6 H61 H62 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NC')
                idx += 1
                sel = f"name C2 H2 N3 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'ND')
                idx += 1
            elif resname == 'GUA':
                sel = f"name N9 C4 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NA')
                idx += 1
                sel = f"name C8 H8 N7 C5 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NB')
                idx += 1
                sel = f"name C6 N1 N6 H1 O6 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NC')
                idx += 1
                sel = f"name C2 N2 H21 H22 N3 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'ND')
                idx += 1
            elif resname == 'CYT':
                sel = f"name N1 C5 H5 C6 H6 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NA')
                idx += 1
                sel = f"name C4 N4 H41 H42 N3 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NB')
                idx += 1
                sel = f"name C2 O2 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NC')
                idx += 1
            elif resname == 'URA':
                sel = f"name N1 C5 H5 C6 H6 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NA')
                idx += 1
                sel = f"name C4 O4 H3 H3 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NB')
                idx += 1
                sel = f"name C2 O2 and resid {resid} and segid {segid}"
                aa2cg(sel, resid, segid, 'NC')
                idx += 1
    print('END', file=f)
print('Finished!')
quit()
