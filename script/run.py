from __future__ import division, print_function
import argparse
# OpenMM Imports
from openmm.unit import *
from openmm.app import *
from openmm import *
import numpy as np


# define the iConRNA force field
###### for RNA System with A-U/G-C pairs ######
def SimpleRNASystem(psf, system, ffs):
    top = psf.topology
    # 2) constructe the force field
    print('\n################# constructe the HyRes force field ####################')
    # get nonbonded force
    for force_index, force in enumerate(system.getForces()):
        if force.getName() == "NonbondedForce":
            nbforce = force
            nbforce_index = force_index
        elif force.getName() == "HarmonicAngleForce":
            hmangle = force
            hmangle_index = force_index
    print('\n# get the NonBondedForce and HarmonicAngleForce:', nbforce.getName(), hmangle.getName())
    
    print('\n# get bondlist')
    # get bondlist
    bondlist = []
    for bond in top.bonds():
        bondlist.append([bond[0].index, bond[1].index])
    #get all atom name
    atoms = []
    for atom in psf.topology.atoms():
        atoms.append(atom.name)
    
    print('\n# replace HarmonicAngle with Restricted Bending (ReB) potential')
    # Custom Angle Force
    ReB = CustomAngleForce("kt*(theta-theta0)^2/(sin(theta)^2);")
    ReB.setName('ReBAngleForce')
    ReB.addPerAngleParameter("theta0")
    ReB.addPerAngleParameter("kt")
    for angle_idx in range(hmangle.getNumAngles()):
        ang = hmangle.getAngleParameters(angle_idx)
        ReB.addAngle(ang[0], ang[1], ang[2], [ang[3], ang[4]])
    system.addForce(ReB)
    
    print('\n# add custom nonbondedforce')
    dh = ffs['dh']
    lmd = ffs['lmd']
    er = ffs['er']
    # add custom nonbondedforce: CNBForce, here only charge-charge interactions
    formula = f"""138.935456/er*charge1*charge2/r*exp(-r/dh)*kpmg;
                dh={dh.value_in_unit(unit.nanometer)}; er={er}; kpmg=select(lb1+lb2,1,lmd); lmd={lmd}
              """
    CNBForce = CustomNonbondedForce(formula)
    CNBForce.setName("LJ_ElecForce")
    CNBForce.setNonbondedMethod(nbforce.getNonbondedMethod())
    CNBForce.setUseSwitchingFunction(use=True)
    #CNBForce.setUseLongRangeCorrection(use=True)
    CNBForce.setCutoffDistance(1.8*unit.nanometers)
    CNBForce.setSwitchingDistance(1.6*unit.nanometers)
    CNBForce.addPerParticleParameter('charge')
    CNBForce.addPerParticleParameter('lb')
    
    for idx in range(nbforce.getNumParticles()):
        particle = nbforce.getParticleParameters(idx)
        if atoms[idx] == 'P':
            lb = 1
        elif atoms[idx] == 'MG':
            lb = -1
        else:
            lb = 2
        perP = [particle[0], lb]
        CNBForce.addParticle(perP)
    
    CNBForce.createExclusionsFromBonds(bondlist, 2)
    system.addForce(CNBForce)
    
    
    print('\n# add base stacking force')
    # base stakcing and paring
    # define relative strength of base pairing and stacking
    eps_base = ffs['eps_base']
    scales = {'AA':1.0, 'AG':1.0, 'AC':0.8, 'AU':0.8, 'GA':1.0, 'GG':1.0, 'GC':1.0, 'GU':1.0,
              'CA':0.4, 'CG':0.5, 'CC':0.5, 'CU':0.3, 'UA':0.3, 'UG':0.3, 'UC':0.2, 'UU':0.0,
              'A-U':0.395, 'C-G':0.545}
    # get all the groups of bases
    grps = []
    for atom in psf.topology.atoms():
        if atom.name == "NA":
            if atom.residue.name in ['A', 'G']:
                grps.append([atom.residue.name, [atom.index, atom.index+1]])
                grps.append([atom.residue.name, [atom.index+2, atom.index+3]])
            elif atom.residue.name in ['C', 'U']:
                grps.append([atom.residue.name, [atom.index, atom.index+1]])
                grps.append([atom.residue.name, [atom.index+1, atom.index+2]])
    # base stacking
    fstack = CustomCentroidBondForce(2, "eps_stack*(5*(r0/r)^10-6.0*(r0/r)^6); r=distance(g1,g2);")
    fstack.setName('StackingForce')
    fstack.addPerBondParameter('eps_stack')
    fstack.addGlobalParameter('r0', 0.34*unit.nanometers)
    # add all group
    for grp in grps:
        fstack.addGroup(grp[1])
    # get the stacking pairs
    sps = []
    for i in range(0,len(grps)-2,2):
        grp = grps[i]
        pij = grps[i][0] + grps[i+2][0]
        sps.append([[i+1, i+2], scales[pij]*eps_base])
    for sp in sps:
        fstack.addBond(sp[0], [sp[1]])
    print('    add ', fstack.getNumBonds(), 'stacking pairs')
    system.addForce(fstack)
    
    # base pairing
    print('\n# add base pair force')
    a_b, a_c, a_d = [], [], []
    g_b, g_c, g_d = [], [], []
    c_a, c_b, c_c, u_a, u_b, u_c = [], [], [], [], [], []
    a_p, g_p, c_p, u_p = [], [], [], []
    num_A, num_G, num_C, num_U = 0, 0, 0, 0
    for atom in psf.topology.atoms():
        if atom.residue.name == 'A':
            num_A += 1
            if atom.name == 'NC':
                a_c.append(int(atom.index))
            elif atom.name == 'NB':
                a_b.append(int(atom.index))
            elif atom.name == 'ND':
                a_d.append(int(atom.index))
            elif atom.name == 'P':
                a_p.append(int(atom.index))
        elif atom.residue.name == 'G':
            num_G += 1
            if atom.name == 'NC':
                g_c.append(int(atom.index))
            elif atom.name == 'NB':
                g_b.append(int(atom.index))
            elif atom.name == 'ND':
                g_d.append(int(atom.index))
            elif atom.name == 'P':
                g_p.append(int(atom.index))
        elif atom.residue.name == 'U':
            num_U += 1
            if atom.name == 'NA':
                u_a.append(int(atom.index))
            elif atom.name == 'NB':
                u_b.append(int(atom.index))
            elif atom.name == 'NC':
                u_c.append(int(atom.index))
            elif atom.name == 'P':
                u_p.append(int(atom.index))
        elif atom.residue.name == 'C':
            num_C += 1
            if atom.name == 'NA':
                c_a.append(int(atom.index))
            elif atom.name == 'NB':
                c_b.append(int(atom.index))
            elif atom.name == 'NC':
                c_c.append(int(atom.index))
            elif atom.name == 'P':
                c_p.append(int(atom.index))
    # add A-U pair through CustomHbondForce
    eps_AU = eps_base*scales['A-U']
    r_au = 0.304*unit.nanometer
    r_au2 = 0.37*unit.nanometer
    
    if num_A != 0 and num_U != 0:
        formula = f"""eps_AU*(5.0*(r_au/r)^10-6.0*(r_au/r)^6 + 5*(r_au2/r2)^10-6.0*(r_au2/r2)^6)*step(cos3)*cos3;
                  r=distance(a1,d1); r2=distance(a3,d2); cos3=-2*cos(phi)^3; phi=angle(d1,a1,a2);
                  eps_AU={eps_AU.value_in_unit(unit.kilojoule_per_mole)};
                  r_au={r_au.value_in_unit(unit.nanometer)}; r_au2={r_au2.value_in_unit(unit.nanometer)}
                  """
        pairAU = CustomHbondForce(formula)
        pairAU.setName('AUpairForce')
        pairAU.setNonbondedMethod(nbforce.getNonbondedMethod())
        pairAU.setCutoffDistance(0.65*unit.nanometer)
        for idx in range(len(a_c)):
            pairAU.addAcceptor(a_c[idx], a_b[idx], a_d[idx])
        for idx in range(len(u_b)):
            pairAU.addDonor(u_b[idx], u_c[idx], -1)
        system.addForce(pairAU)
        print(pairAU.getNumAcceptors(), pairAU.getNumDonors(), 'AU')
        
    # add C-G pair through CustomHbondForce
    eps_CG = eps_base*scales['C-G']
    r_cg = 0.304*unit.nanometer
    r_cg2 = 0.35*unit.nanometer
    
    if num_C != 0 and num_G != 0:
        formula = f"""eps_CG*(5.0*(r_cg/r)^10-6.0*(r_cg/r)^6 + 5*(r_cg2/r2)^10-6.0*(r_cg2/r2)^6)*step(cos3)*cos3;
                  r=distance(a1,d1); r2=distance(a3,d2); cos3=-2*cos(phi)^3; phi=angle(d1,a1,a2); psi=dihedral(a3,a1,d1,d2);
                  eps_CG={eps_CG.value_in_unit(unit.kilojoule_per_mole)};
                  r_cg={r_cg.value_in_unit(unit.nanometer)}; r_cg2={r_cg2.value_in_unit(unit.nanometer)}
                  """
        pairCG = CustomHbondForce(formula)
        pairCG.setName('CGpairForce')
        pairCG.setNonbondedMethod(nbforce.getNonbondedMethod())
        pairCG.setCutoffDistance(0.65*unit.nanometer)
        for idx in range(len(g_c)):
            pairCG.addAcceptor(g_c[idx], g_b[idx], g_d[idx])
        for idx in range(len(c_b)):
            pairCG.addDonor(c_b[idx], c_c[idx], -1)
        system.addForce(pairCG)
        print(pairCG.getNumAcceptors(), pairCG.getNumDonors(), 'CG')
    
    # delete the NonbondedForce and HarmonicAngleForce
    system.removeForce(nbforce_index)
    system.removeForce(hmangle_index)
    return system


# 0) set variables in the simulation
gpu_id = "0"
top_inp, param_inp = 'param_iConRNA.inp', 'top_iConRNA.inp'

# input parameters
parser = argparse.ArgumentParser()
parser.add_argument('-c', "--pdb", default='conf.pdb', help="pdb file, default is conf.pdb")
parser.add_argument('-p', "--psf", default='conf.psf', help="psf file, default is conf.psf")
parser.add_argument('-t', "--temp", default=303, type=float, help="system temperature, default is 303 K")
parser.add_argument('-b', "--box", nargs='+', type=float, help="box dimensions in nanometer, e.g., '50 50 50' ")
parser.add_argument('-s', "--salt", default=150.0, type=float, help="salt concentration in mM, default is 0.0 mM")
parser.add_argument('-e', "--ens", default='NVT', type=str, help="simulation ensemble, NPT, NVT, or non, non is for non-periodic system")
parser.add_argument('-m', "--Mg", default=0.0, type=float, help="Mg2+ concentration in mM")

args = parser.parse_args()
pdb_file = args.pdb
psf_file = args.psf
T = args.temp
c_ion = args.salt/1000.0                                   # concentration of ions in M
c_Mg = args.Mg                                           # concentration of Mg in mM
ensemble = args.ens

## set pbc and box vector
if ensemble == 'non' and c_Mg != 0.0:
    print("Error: Mg ion cannot be usde in non-periodic system.")
    exit(1)
if ensemble in ['NPT', 'NVT']:
    # pbc box length
    if len(args.box) == 1:
        lx, ly, lz = args.box[0], args.box[0], args.box[0]
    elif len(args.box) == 3:
        lx = args.box[0]
        ly = args.box[1]
        lz = args.box[2]
    else:
        print("Error: You must provide either one or three values for box.")
        exit(1)
    a = Vec3(lx, 0.0, 0.0)
    b = Vec3(0.0, ly, 0.0)
    c = Vec3(0.0, 0.0, lz)
elif ensemble not in ['NPT', 'NVT', 'non']:
    print("Error: The ensemble must be NPT, NVT or non. The input value is {}.".format(ensemble))
    exit(1)

# simulation parameters
dt_equil = 0.001*unit.picoseconds		                # equilibration time step
dt_prod = 0.008*unit.picoseconds				# production time step
total_step = 250000000                                             # total step
equil_step = 10000
temperture = T*unit.kelvin                                      # temperature
log_freq = 1000                                                 # frequency of log file
dcd_freq = 10000                                               # frequency of dcd file
pdb_freq = 12500000                                             # frequence of dpd_traj file
pressure = 1*unit.atmosphere                                    # pressure in NPT
friction = 0.1/unit.picosecond                                  # friction coefficient in Langevin

# force field parameters
Td = T-273
er_t = 87.74-0.4008*Td+9.398*10**(-4)*Td**2-1.41*10**(-6)*Td**3
print('relative electric constant: ', er_t*20.3/77.6)                        
dh = 0.304/(np.sqrt(c_ion))
print('Debye-Huckel screening length: ', dh)
if c_Mg == 0:
    nMg = 0
else:
    #### modify the nMg calculation for different RNA
    ### for CAG/CUG
    nMg = 0.46*(c_Mg/0.42)**(0.408)/(1+(c_Mg/0.42)**(0.408)) + 0.0012*(Td-30)    
    ### poly(rA)
    # nMg = 0.54*(c_Mg/0.94)**(0.59)/(1+(c_Mg/0.94)**(0.59)) + 0.0012*(Td-30)
    ### poly(rU)
    # nMg = 0.48*(c_Mg/1.31)**(0.85)/(1+(c_Mg/1.31)**(0.85)) + 0.0012*(Td-30)  
                                                 
lmd0 = 1.255*(nMg/0.172)**0.625/(1+(nMg/0.172)**0.625)
print('lmd: ', lmd0)
ffs = {
    'temp': T,                                                  # Temperature
    'lmd': lmd0,                                                # Charge scaling factor of P-
    'dh': dh*unit.nanometer,                                  # Debye Huckel screening length
    'ke': 138.935456,                                           # Coulomb constant, ONE_4PI_EPS0
    'er': er_t*20.3/77.6,                                         # relative dielectric constant
    'eps_hb': 2.0*unit.kilocalorie_per_mole,                    # hydrogen bond strength
    'sigma_hb': 0.29*unit.nanometer,                            # sigma of hydrogen bond
    'eps_base': 2.05*unit.kilocalorie_per_mole,                 # base stacking strength
}

# 1) import coordinates and topology form charmm pdb and psf
print('\n################## load coordinates, topology and parameters ###################')
pdb = PDBFile(pdb_file)
psf = CharmmPsfFile(psf_file)
top = psf.topology
params = CharmmParameterSet(top_inp, param_inp)
if ensemble == 'non':
    system = psf.createSystem(params, nonbondedMethod=CutoffNonPeriodic, constraints=HBonds)
else:
    psf.setBox(lx, ly, lz)
    top.setPeriodicBoxVectors((a, b, c))
    top.setUnitCellDimensions((lx, ly,lz))
    system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic, constraints=HBonds)
    system.setDefaultPeriodicBoxVectors(a, b, c)

system = SimpleRNASystem(psf, system, ffs)

with open('system.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))

print('\n# Now, the system has:')
for force in system.getForces():
    print('      ', force.getName())

# simulation
print('\n################### prepare simulation system####################')
if ensemble == 'NPT':
    print('This is a NPT system')
    system.addForce(MonteCarloBarostat(pressure, temperture, 25))
elif ensemble == 'NVT':
    print('This is a NVT system')
elif ensemble == 'non':
    print('This is a non-periodic system')
else:
    print('Only NPT, NVT, and non-periodic system are supported!')
integrator = LangevinMiddleIntegrator(temperture, friction, dt_equil)
plat = Platform.getPlatformByName('CUDA')
prop = {'Precision': 'mixed', 'DeviceIndex': gpu_id}
simulation = Simulation(top, system, integrator, plat, prop)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperture)
print(f'Langevin, CUDA, {temperture}')

print('\n################### Minimization, Equilibriation, Production simulation ####################')
print('# minimizeEnergy:')
print('before: ', simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy(maxIterations=500000, tolerance=0.01)
print('after: ', simulation.context.getState(getEnergy=True).getPotentialEnergy())

print('\n# Equilibriation running:')
simulation.step(equil_step)

## save a pdb traj using large step, dcd traj using small step, and log file
simulation.reporters.append(PDBReporter('system.pdb', pdb_freq))
simulation.reporters.append(XTCReporter('system.xtc', dcd_freq))
#simulation.reporters.append(DCDReporter('system.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter('system.log', log_freq, progress=True, totalSteps=total_step, temperature=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True
))
#simulation.reporters.append(CheckpointReporter('system.chk', dcd_freq*10))

simulation.integrator.setStepSize(dt_prod)
print('\n# Production simulation running:')
simulation.step(total_step)

simulation.saveCheckpoint('system.chk')
print('\n# Finished!')
