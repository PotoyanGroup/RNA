import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import unit
import open3SPN2
import pandas 
import numpy as np
import openmm.unit as u
import parmed
import warnings


##angle
#PSP = np.deg2rad(82.735)
PSP = np.deg2rad(82.735)
SPS = np.deg2rad(87.410)
PSA  = np.deg2rad(97.569)
PSU = np.deg2rad(90.155)
PSG = np.deg2rad(101.356)
PSC = np.deg2rad(90.545)
ASP = np.deg2rad(110.346)
USP = np.deg2rad(112.661)
GSP = np.deg2rad(109.721)
CSP = np.deg2rad(112.615)

##Bond
SP = 3.8157
PS = 4.6010
SA = 4.8515
SU = 4.2733
SG = 4.9659
SC = 4.2738

##dihedral angle
PSPS = np.deg2rad(-148.215)
SPSP = np.deg2rad(175.975)


###H-bond distance
A_U = 5.8815
G_C = 5.6550

###
SA_U = np.deg2rad(156.320)
SU_A = np.deg2rad(143.910)
SG_C = np.deg2rad(161.746)
SC_G = np.deg2rad(142.306)


###H-bond dihedral
SA_US = 71.958
SG_CS = 79.653
PSA_U = 54.689
PSU_A = 67.305
PSG_C = 43.654
PSC_G = 69.752

####stack distance
AC = 3.8260185
AG = 4.4255305
AU  = 3.8260185
CA  = 4.7010580
CC  = 4.2500910
CG  = 4.9790760
CU  = 4.2273615
GA  = 4.0128560
GC  = 3.6784360
GG  = 4.2427250
GU  = 3.6616930
UA  = 4.7010580
UC  = 4.2679180
UG  = 4.9977560
UU  = 4.2453650

###   h       s      Tm
s_AA = [4.348,  -0.319,   298.9]
s_AC = [4.311,  -0.319,   298.9]
s_AG = [5.116,   5.301,   341.2]
s_AU  =[4.311,  -0.319,   298.9]
s_CA  = [4.287,  -0.319,   298.9]
s_CC  = [4.015,  -1.567,   285.8]
s_CG  = [4.602,   0.774,   315.5]
s_CU  = [3.995,  -1.567,   285.8]
s_GA  = [5.079,   5.301,   341.2]
s_GC  = [5.075,   4.370,   343.2]
s_GG  = [5.555,   7.346,   366.3]
s_GU  = [4.977,   2.924 ,  338.2]
s_UA  = [4.287,  -0.319,   298.9]
s_UC  = [3.992,  -1.567,   285.8]
s_UG  = [5.032,   2.924,   338.2]
s_UU  = [3.370,  -3.563,   251.6]

kB=1.380649*10**(-23)
def s_u0(x,T):
    return -x[0] + kB*(T - x[2])* x[1]


ss_gg = s_u0(s_GG,300)
ss_gc = s_u0(s_GC,300)
ss_cg = s_u0(s_CG,300)
ss_cc = s_u0(s_CC,300)

all_ss = {'GG':(GG,ss_gg),'GC':(GC,ss_gc),'CC':(CC,ss_cc),'CG':(CG,ss_cg)} 


bonds_parm = {'SrP':3.8157,'PSr':4.6010,'SrA' : 4.8515,
        'SrU' : 4.2733,'SrG': 4.9659,'SrC' :4.2738,
        'ASr' : 4.8515, 'USr' : 4.2733,'GSr': 4.9659,'CSr' : 4.2738}


angles_parm = {'PSrP':np.deg2rad(82.735),'SrPSr':np.deg2rad(87.410),'PSrA' : np.deg2rad(97.569),
        'PSrU' : np.deg2rad(90.155),'PSrG': np.deg2rad(101.356),'PSrC' :np.deg2rad(90.545),
        'ASrP' : np.deg2rad(110.346), 'USrP' : np.deg2rad(112.661),'GSrP': np.deg2rad(109.721),'CSrP' : np.deg2rad(112.615)}


import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import unit
import open3SPN2
import pandas 
import numpy as np
import openmm.unit as u
import parmed

##angle
#PSP = np.deg2rad(82.735)
PSP = np.deg2rad(82.735)
SPS = np.deg2rad(87.410)
PSA  = np.deg2rad(97.569)
PSU = np.deg2rad(90.155)
PSG = np.deg2rad(101.356)
PSC = np.deg2rad(90.545)
ASP = np.deg2rad(110.346)
USP = np.deg2rad(112.661)
GSP = np.deg2rad(109.721)
CSP = np.deg2rad(112.615)

##Bond
SP = 3.8157
PS = 4.6010
SA = 4.8515
SU = 4.2733
SG = 4.9659
SC = 4.2738

##dihedral angle
PSPS = np.deg2rad(-148.215)
SPSP = np.deg2rad(175.975)


###H-bond distance
A_U = 5.8815
G_C = 5.6550

###
SA_U = np.deg2rad(156.320)
SU_A = np.deg2rad(143.910)
SG_C = np.deg2rad(161.746)
SC_G = np.deg2rad(142.306)


###H-bond dihedral
SA_US = 71.958
SG_CS = 79.653
PSA_U = 54.689
PSU_A = 67.305
PSG_C = 43.654
PSC_G = 69.752

####stack distance
AC = 3.8260185
AG = 4.4255305
AU  = 3.8260185
CA  = 4.7010580
CC  = 4.2500910
CG  = 4.9790760
CU  = 4.2273615
GA  = 4.0128560
GC  = 3.6784360
GG  = 4.2427250
GU  = 3.6616930
UA  = 4.7010580
UC  = 4.2679180
UG  = 4.9977560
UU  = 4.2453650

###   h       s      Tm
s_AA = [4.348,  -0.319,   298.9]
s_AC = [4.311,  -0.319,   298.9]
s_AG = [5.116,   5.301,   341.2]
s_AU  =[4.311,  -0.319,   298.9]
s_CA  = [4.287,  -0.319,   298.9]
s_CC  = [4.015,  -1.567,   285.8]
s_CG  = [4.602,   0.774,   315.5]
s_CU  = [3.995,  -1.567,   285.8]
s_GA  = [5.079,   5.301,   341.2]
s_GC  = [5.075,   4.370,   343.2]
s_GG  = [5.555,   7.346,   366.3]
s_GU  = [4.977,   2.924 ,  338.2]
s_UA  = [4.287,  -0.319,   298.9]
s_UC  = [3.992,  -1.567,   285.8]
s_UG  = [5.032,   2.924,   338.2]
s_UU  = [3.370,  -3.563,   251.6]

kB=1.380649*10**(-23)
def s_u0(x,T):
    return -x[0] + kB*(T - x[2])* x[1]


ss_gg = s_u0(s_GG,300)
ss_gc = s_u0(s_GC,300)
ss_cg = s_u0(s_CG,300)
ss_cc = s_u0(s_CC,300)

all_ss = {'GG':(GG,ss_gg),'GC':(GC,ss_gc),'CC':(CC,ss_cc),'CG':(CG,ss_cg)} 


import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import numpy as np
import MDAnalysis as mda
import parmed
import simtk
import sys


class CustomRNAForceSystem:
    def __init__(self, pdb_file, psf_file, forcefield_file):
        self.pdb = app.PDBFile(pdb_file)
        self.coord=self.pdb.positions
        
        self.forcefield = app.ForceField(forcefield_file)
        self.structure = parmed.load_file(psf_file)

        self.topology = self.pdb.topology
        self.topology.setPeriodicBoxVectors([[10.0*unit.nanometers,0.0*unit.nanometers,0.0*unit.nanometers], 
                                [0.0*unit.nanometers,10.0*unit.nanometers,0.0*unit.nanometers], 
                                [0.0*unit.nanometers,0.0*unit.nanometers,10.0*unit.nanometers]])
        
        self.system = self._create_system()
        
    def _create_system(self):
        return self.forcefield.createSystem(self.topology)

    def remove_nonbonded_force(self):
        for i,force in enumerate(self.system.getForces()):
            self.system.removeForce(i)

        
    def add_bonds(self, bonds_parm):
        bond_force = mm.HarmonicBondForce()
        bond_force.setUsesPeriodicBoundaryConditions(True)
        
        for bond in self.structure.bonds:
            bond_length = bonds_parm.get(bond.atom1.name + bond.atom2.name)
            if bond_length:
                bond_force.addBond(
                    bond.atom1.idx, bond.atom2.idx,
                    bond_length * unit.angstrom,
                    20.0 * unit.kilocalorie_per_mole / (unit.angstrom**2)
                )
        bond_force.setForceGroup(0)
        self.system.addForce(bond_force)
        
    def add_angles(self, angles_parm):
        angle_force = mm.HarmonicAngleForce()
        angle_force.setUsesPeriodicBoundaryConditions(True)
        
        for angle in self.structure.angles:
            angle_val = angles_parm.get(
                angle.atom1.name + angle.atom2.name + angle.atom3.name
            )
            if angle_val:
                angle_force.addAngle(
                    angle.atom1.idx, angle.atom2.idx, angle.atom3.idx,
                    angle_val * unit.radian,
                    5.0 * unit.kilocalorie_per_mole / (unit.radian**2)
                )
        angle_force.setForceGroup(1)
        self.system.addForce(angle_force)

    def add_excluded_volume(self):
        
        ex_vol_force = mm.CustomNonbondedForce(
            'eps*((val/(r - val))^12 - 2.0*(val/(r - val))^6 + 1.0)'
        )
        ex_vol_force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        ex_vol_force.addGlobalParameter('eps', 0.25 * unit.nanometer)
        ex_vol_force.addGlobalParameter("val", 0.16 * unit.nanometer)
        ex_vol_force.setCutoffDistance(0.32 * unit.nanometer)
        
        for atom in self.topology.atoms():
            ex_vol_force.addParticle([])
        
        for bond in self.structure.bonds:
            ex_vol_force.addExclusion(bond.atom1.idx, bond.atom2.idx)
            
        ex_vol_force.setForceGroup(2)
        self.system.addForce(ex_vol_force)

    def Electrostatics(self,T,C):
        from openmm import unit
        e = 249.4 - 0.788 * (T / unit.kelvin) + 7.2E-4 * (T / unit.kelvin) ** 2
        a = 1 - 0.2551 * (C / unit.molar) + 5.151E-2 * (C / unit.molar) ** 2 - 6.889E-3 * (C / unit.molar) ** 3
        #print(e, a)
        dielectric = e * a
        # Debye length
        kb = unit.BOLTZMANN_CONSTANT_kB  # Bolztmann constant
        Na = unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * unit.coulomb  # proton charge
        pv = 8.8541878176E-12 * unit.farad / unit.meter  # dielectric permittivity of vacuum

        ldby = np.sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        ldby = ldby.in_units_of(unit.nanometer)
        denominator = 4 * np.pi * pv * dielectric / (Na * ec ** 2)
        denominator = denominator.in_units_of(unit.kilocalorie_per_mole**-1 * unit.nanometer**-1)
        #print(ldby, denominator)

        electrostaticForce = openmm.CustomNonbondedForce("""energy;
                                                                energy=q1*q2*exp(-r/dh_length)/denominator/r;""")
        electrostaticForce.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('dh_length', ldby)
        electrostaticForce.addGlobalParameter('denominator', denominator)
        electrostaticForce.setCutoffDistance(1* unit.nanometer)
        
        for i in self.structure.atoms:
            
            if i.name =='P':
                electrostaticForce.addParticle([1])
                
            else:
                electrostaticForce.addParticle([0])
            
        for bond in self.structure.bonds:
            electrostaticForce.addExclusion(bond.atom1.idx, bond.atom2.idx)
                
        electrostaticForce.setForceGroup(4)
        self.system.addForce(electrostaticForce)
        
    def add_hbonds_gc(self,pdb_file):
       
        list_G,list_C,list_A,list_U = self.find_hlist(pdb_file)
        
        Hbondingforce_gc = mm.CustomHbondForce(self._hydrogen_bond_energy_expression())
        Hbondingforce_gc.setCutoffDistance(1.8*nanometers)  # Paper
        Hbondingforce_gc.addGlobalParameter('phi0',np.deg2rad(42.98)*radians)
        Hbondingforce_gc.addGlobalParameter('sigma',0.533*nanometers)
        Hbondingforce_gc.addGlobalParameter('t01',np.deg2rad(42.98)*radians)
        Hbondingforce_gc.addGlobalParameter('t02',np.deg2rad(42.98)*radians)
        Hbondingforce_gc.addGlobalParameter('rng',12)
        Hbondingforce_gc.addGlobalParameter('epsilon',21*kilojoule_per_mole)
        Hbondingforce_gc.addGlobalParameter('alpha',20/nanometers)
        Hbondingforce_gc.addGlobalParameter('pi', np.pi)
        Hbondingforce_gc.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)

        for donor in list_G:
            Hbondingforce_gc.addDonor(donor[0],donor[1],donor[2])
            
        for acceptor in list_C:
            Hbondingforce_gc.addAcceptor(acceptor[0],acceptor[1],acceptor[2])

        ####  Hydrogen bond between G-C
        for ind1,donor in enumerate(list_G):
            for ind2,acceptor in enumerate(list_C):
                res_no_1  = donor[3]
                res_no_2  = acceptor[3]
                chain_1  = donor[4]
                chain_2  = acceptor[4]        
                if abs(res_no_1-res_no_2)<3 and chain_1 == chain_2:   ### to avoid neighbour one
                    Hbondingforce_gc.addExclusion(ind1,ind2)


        Hbondingforce_gc.setForceGroup(3)         
        self.system.addForce(Hbondingforce_gc)

    def _hydrogen_bond_energy_expression(self):
               
        H_energy = '''energy;
                        energy=rep+1/2*(1+cos(dphi))*fdt1*fdt2*attr;
                        rep  = epsilon*(1-exp(-alpha*dr))^2*(1-step(dr));
                        attr = epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                        fdt1 = max(f1*pair0t1,pair1t1);
                        fdt2 = max(f2*pair0t2,pair1t2);
                        pair1t1 = step(pi/2+dt1)*step(pi/2-dt1);
                        pair1t2 = step(pi/2+dt2)*step(pi/2-dt2);
                        pair0t1 = step(pi+dt1)*step(pi-dt1);
                        pair0t2 = step(pi+dt2)*step(pi-dt2);
                        f1 = 1-cos(dt1)^2;
                        f2 = 1-cos(dt2)^2;
                        dphi = dihedral(d2,d1,a1,a2)-phi0;
                        dr    = distance(d1,a1)-sigma;
                        dt1   = rng*(angle(d2,d1,a1)-t01);
                        dt2   = rng*(angle(a2,a1,d1)-t02);'''
        # Long string representing the hydrogen bond energy expression
        return H_energy
        
    def add_gc_ang(self,pdb_file):
        
        list_G,list_C,list_A,list_U = self.find_hlist(pdb_file)
        
        Hbondingforce_gc = mm.CustomHbondForce(self._hbond_mini())
        #Hbondingforce_gc.addGlobalParameter('pi_cut',np.pi*radians)
        Hbondingforce_gc.addGlobalParameter('scale1',1.0)
        Hbondingforce_gc.addGlobalParameter('kdist', 1.5/angstroms**2)
        Hbondingforce_gc.addGlobalParameter('kang', 1.5/radians**2)
        Hbondingforce_gc.addGlobalParameter('kdih', 0.15/radians**2)

        Hbondingforce_gc.addGlobalParameter("UHyd",-5)
        Hbondingforce_gc.addGlobalParameter("dist_0",G_C)
        Hbondingforce_gc.addGlobalParameter("angle_01",SG_C)
        Hbondingforce_gc.addGlobalParameter("angle_02",SC_G)
        Hbondingforce_gc.addGlobalParameter("dihedral_01",SG_CS)
        Hbondingforce_gc.addGlobalParameter('pi', np.pi)

        cutoff = 15*unit.angstroms
        Hbondingforce_gc.setCutoffDistance(cutoff)
        Hbondingforce_gc.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)

        for acceptor in list_G:
            Hbondingforce_gc.addAcceptor(acceptor[0],acceptor[1],acceptor[2])

        for donor in list_C:
            Hbondingforce_gc.addDonor(donor[0],donor[1],donor[2])

        ####  Hydrogen bond between G-C
        for ind1,acceptor in enumerate(list_G):
            for ind2,donor in enumerate(list_C):
                res_no_1  = donor[3]
                res_no_2  = acceptor[3]
                chain_1  = donor[4]
                chain_2  = acceptor[4]        
                if abs(res_no_1-res_no_2)<3 and chain_1 == chain_2:   ### to avoid neighbour one
                    Hbondingforce_gc.addExclusion(ind2,ind1)


        Hbondingforce_gc.setForceGroup(3)         
        self.system.addForce(Hbondingforce_gc)
        print('added gc base pair force')
        
    def add_gc_ang2(self,pdb_file):
        
        list_G,list_C,list_A,list_U = self.find_hlist(pdb_file)
        
        Hbondingforce_gc = mm.CustomHbondForce(self._hbond_ang())
        #Hbondingforce_gc.addGlobalParameter('pi_cut',np.pi*radians)
        Hbondingforce_gc.addGlobalParameter('scale1',1.0)
        Hbondingforce_gc.addGlobalParameter('kdist', 5.0/angstroms**2)
        Hbondingforce_gc.addGlobalParameter('kang', 1.5/radians**2)

        Hbondingforce_gc.addGlobalParameter("UHyd",5.0 * unit.kilocalories_per_mole)
        Hbondingforce_gc.addGlobalParameter("dist_0",G_C)
        Hbondingforce_gc.addGlobalParameter("angle_01",SG_C)
        Hbondingforce_gc.addGlobalParameter("angle_02",SC_G)
        Hbondingforce_gc.addGlobalParameter('pi', np.pi)
        cutoff = 15*unit.angstroms
        Hbondingforce_gc.setCutoffDistance(cutoff)
        Hbondingforce_gc.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)

        for acceptor in list_G:
            Hbondingforce_gc.addAcceptor(acceptor[0],acceptor[1],acceptor[2])

        for donor in list_C:
            Hbondingforce_gc.addDonor(donor[0],donor[1],donor[2])

        ####  Hydrogen bond between G-C
        for ind1,acceptor in enumerate(list_G):
            for ind2,donor in enumerate(list_C):
                res_no_1  = donor[3]
                res_no_2  = acceptor[3]
                chain_1  = donor[4]
                chain_2  = acceptor[4]        
                if abs(res_no_1-res_no_2)<3 and chain_1 == chain_2:   ### to avoid neighbour one
                    Hbondingforce_gc.addExclusion(ind2,ind1)


        Hbondingforce_gc.setForceGroup(3)         
        self.system.addForce(Hbondingforce_gc)
        print('added gc base pair force')

    def add_hbonds_gc_lj_ang(self,pdb_file):
        
        list_G,list_C,list_A,list_U = self.find_hlist(pdb_file)
        
        Hbondingforce_gc = mm.CustomHbondForce(self._hbond_lj_ang2())
        Hbondingforce_gc.addGlobalParameter('r_HB',10.0 * unit.angstroms )
        Hbondingforce_gc.addGlobalParameter('delta_HB', 1.5 * unit.kilocalories_per_mole)
        Hbondingforce_gc.addGlobalParameter('sigma_HB', 5.6 * unit.angstroms)
        
        cutoff = 15*unit.angstroms
        Hbondingforce_gc.setCutoffDistance(cutoff)
        Hbondingforce_gc.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)

        for acceptor in list_G:
            Hbondingforce_gc.addAcceptor(acceptor[0],acceptor[1],acceptor[2])

        for donor in list_C:
            Hbondingforce_gc.addDonor(donor[0],donor[1],donor[2])

        ####  Hydrogen bond between G-C
        for ind1,acceptor in enumerate(list_G):
            for ind2,donor in enumerate(list_C):
                res_no_1  = donor[3]
                res_no_2  = acceptor[3]
                chain_1  = donor[4]
                chain_2  = acceptor[4]        
                if abs(res_no_1-res_no_2)<5 and chain_1 == chain_2:   ### to avoid neighbour one
                    Hbondingforce_gc.addExclusion(ind2,ind1)


        Hbondingforce_gc.setForceGroup(3)         
        self.system.addForce(Hbondingforce_gc)
        print('added lj with ANG force')  
        
    def add_hbonds_gc_lj(self,pdb_file):
        
        list_G,list_C,list_A,list_U = self.find_hlist(pdb_file)
        
        Hbondingforce_gc = mm.CustomHbondForce(self._hbond_lj())
        Hbondingforce_gc.addGlobalParameter('r_HB',7.0 * unit.angstroms )
        Hbondingforce_gc.addGlobalParameter('delta_HB', 2.0 * unit.kilocalories_per_mole)
        Hbondingforce_gc.addGlobalParameter('sigma_HB', 5.6 * unit.angstroms)
        
        cutoff = 15*unit.angstroms
        Hbondingforce_gc.setCutoffDistance(cutoff)
        Hbondingforce_gc.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)

        for acceptor in list_G:
            Hbondingforce_gc.addAcceptor(acceptor[0],acceptor[1],acceptor[2])

        for donor in list_C:
            Hbondingforce_gc.addDonor(donor[0],donor[1],donor[2])

        ####  Hydrogen bond between G-C
        for ind1,acceptor in enumerate(list_G):
            for ind2,donor in enumerate(list_C):
                res_no_1  = donor[3]
                res_no_2  = acceptor[3]
                chain_1  = donor[4]
                chain_2  = acceptor[4]        
                if abs(res_no_1-res_no_2)<3 and chain_1 == chain_2:   ### to avoid neighbour one
                    Hbondingforce_gc.addExclusion(ind2,ind1)


        Hbondingforce_gc.setForceGroup(3)         
        self.system.addForce(Hbondingforce_gc)
        print('added lj force')
        
    def add_hbonds_gc_lj2(self,pdb_file):
        
        list_G,list_C,list_A,list_U = self.find_hlist(pdb_file)
        
        Hbondingforce_gc = mm.CustomHbondForce(self._hbond_lj_ang())
        Hbondingforce_gc.addGlobalParameter('r_HB',7.0 * unit.angstroms )
        Hbondingforce_gc.addGlobalParameter('delta_HB', 2.0 * unit.kilocalories_per_mole)
        Hbondingforce_gc.addGlobalParameter('sigma_HB', 5.6 * unit.angstroms)
        Hbondingforce_gc.addGlobalParameter('ang1_cut',np.deg2rad(161.746-30)*radians )
        Hbondingforce_gc.addGlobalParameter('ang2_cut',np.deg2rad(142.306-30)*radians )
        
        cutoff = 15*unit.angstroms
        Hbondingforce_gc.setCutoffDistance(cutoff)
        Hbondingforce_gc.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)

        for acceptor in list_G:
            Hbondingforce_gc.addAcceptor(acceptor[0],acceptor[1],acceptor[2])

        for donor in list_C:
            Hbondingforce_gc.addDonor(donor[0],donor[1],donor[2])

        ####  Hydrogen bond between G-C
        for ind1,acceptor in enumerate(list_G):
            for ind2,donor in enumerate(list_C):
                res_no_1  = donor[3]
                res_no_2  = acceptor[3]
                chain_1  = donor[4]
                chain_2  = acceptor[4]        
                if abs(res_no_1-res_no_2)<3 and chain_1 == chain_2:   ### to avoid neighbour one
                    Hbondingforce_gc.addExclusion(ind2,ind1)


        Hbondingforce_gc.setForceGroup(3)         
        self.system.addForce(Hbondingforce_gc)
        print('added lj force')
        
    def add_hbonds_gc_exp(self,pdb_file):
        
        list_G,list_C,list_A,list_U = self.find_hlist(pdb_file)        
        Hbondingforce_gc = mm.CustomHbondForce(self._hydrogen_bond_eq())
        Hbondingforce_gc.addGlobalParameter('Uhb', -10*kilocalorie_per_mole)
        Hbondingforce_gc.addGlobalParameter('r0',G_C)
        Hbondingforce_gc.addGlobalParameter('kr',1./angstrom**2)
        Hbondingforce_gc.addGlobalParameter('kt',0/radian**2)
        Hbondingforce_gc.addGlobalParameter('kp', 0.0)
        Hbondingforce_gc.addGlobalParameter('theta1',SC_G*radians)
        Hbondingforce_gc.addGlobalParameter('theta2',SG_C*radians)
        Hbondingforce_gc.addGlobalParameter('phi1',SG_CS*radians)
                
        cutoff = 15*unit.angstroms
        Hbondingforce_gc.setCutoffDistance(cutoff)
        Hbondingforce_gc.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)

        for donor in list_G:
            Hbondingforce_gc.addDonor(donor[0],donor[1],donor[2])
            
        for acceptor in list_C:
            Hbondingforce_gc.addAcceptor(acceptor[0],acceptor[1],acceptor[2])

        ####  Hydrogen bond between G-C
        for ind1,donor in enumerate(list_G):
            for ind2,acceptor in enumerate(list_C):
                res_no_1  = donor[3]
                res_no_2  = acceptor[3]
                chain_1  = donor[4]
                chain_2  = acceptor[4]        
                if abs(res_no_1-res_no_2)<3 and chain_1 == chain_2:   ### to avoid neighbour one
                    Hbondingforce_gc.addExclusion(ind1,ind2)


        Hbondingforce_gc.setForceGroup(3)         
        self.system.addForce(Hbondingforce_gc)
        print('added gc expe force')      

       
    def _hydrogen_bond_eq(self):
                
        energy_function =  "- kr*(distance(a1, d1) - r0)^2"
        energy_function += "- kt*(angle(a2, a1, d1) - theta1)^2"
        energy_function += "- kt*(angle(d2, d1, a1) - theta1)^2"
        energy_function += "- kp*(1. + cos(dihedral(d2, d1, a1, a2) + phi1))"
        energy_function = "2 * Uhb * exp(" + energy_function + ")"
        
        return energy_function
  

    
 
    def _hbond_lj(self):
        H_energy = "step(r_HB - r) * delta_HB * (5 * (sigma_HB/r)^12 - 6 * (sigma_HB/r)^10); r= distance(a1,d1)"
        
        return H_energy
    
    
    def _hbond_lj_ang(self):
        H_energy = "step(r_HB - r) * delta_HB * (5 * (sigma_HB/r)^12 - 6 * (sigma_HB/r)^10)* "
        H_energy += "step(angle(a2, a1, d1) - ang1_cut) * step(angle(d2, d1, a1) - ang2_cut); r= distance(a1,d1)"
        
        return H_energy
    
    def _hbond_mini(self):  
        
        H_energy = "scale1*UHyd/(1.0 + kdist*DIST*DIST + kang*ANG1*ANG1 + "
        H_energy += "kang*ANG2*ANG2 + kdih*DIHED1*DIHED1); "
        #H_energy += "scale1 = 0.25*(1.0 - ((angle(a2,a1,d1)-pi_cut)/(abs(angle(a2,a1,d1)-pi_cut))))*(1.0 - ((angle(d2,d1,a1)-pi_cut)/(abs(angle(d2,d1,a1)-pi_cut))));"
        H_energy += "DIST = (distance(a1,d1) - dist_0); ANG1 = (angle(a2,a1,d1) - angle_01); ANG2 = (angle(d2,d1,a1) - angle_02);"
        #H_energy += "DIHED1 = dihedral(a2,a1,d1,d2) - dihedral_01 + pi*(((pi - dihedral(a2,a1,d1,d2)"
        #H_energy += "+ dihedral_01)/abs(pi - dihedral(a2,a1,d1,d2) + dihedral_01)) - ((pi + dihedral(a2,a1,d2,d1) - dihedral_01)/abs(pi + dihedral(a2,a1,d1,d2) - dihedral_01))); "
        
        return H_energy

    def _hbond_ang(self):  
        
        H_energy = "scale1*UHyd/(1.0 + kdist * DIST^2 + kang * ANG1^2 + "
        H_energy += "kang *ANG2^2); "
        #H_energy += "scale1 = 0.25*(1.0 - ((angle(a2,a1,d1)-pi_cut)/(abs(angle(a2,a1,d1)-pi_cut))))*(1.0 - ((angle(d2,d1,a1)-pi_cut)/(abs(angle(d2,d1,a1)-pi_cut))));"
        H_energy += "DIST = (distance(a1,d1) - dist_0); ANG1 = (angle(a2,a1,d1) - angle_01); ANG2 = (angle(d2,d1,a1) - angle_02);"
        #H_energy += "DIHED1 = dihedral(a2,a1,d1,d2) - dihedral_01 + pi*(((pi - dihedral(a2,a1,d1,d2)"
        #H_energy += "+ dihedral_01)/abs(pi - dihedral(a2,a1,d1,d2) + dihedral_01)) - ((pi + dihedral(a2,a1,d2,d1) - dihedral_01)/abs(pi + dihedral(a2,a1,d1,d2) - dihedral_01))); "
        
        return H_energy
        # Long string representing the hydrogen bond energy expression



    def _hbond_lj_ang3(self):
        
        H_energy = "step(r_HB - r) * delta_HB * (5 * (sigma_HB/r)^12 - 6 * (sigma_HB/r)^10)"
        H_energy += "/(1+kang * ANG1^2 + kang *ANG2^2);r= distance(a1,d1); ANG1 = (angle(a2,a1,d1) - angle_01); "
        H_energy += "ANG2 = (angle(d2,d1,a1) - angle_02); "
        
        return H_energy
    
    def _hbond_lj_ang2(self):
        
        H_energy = "step(r_HB - r) * delta_HB * (5 * (sigma_HB/r)^12 - 6 * (sigma_HB/r)^10)"
        H_energy += "*step(cosA)*cosA * step(cosB)*cosB;r= distance(a1,d1); cosA=-2*cos(ANG1)^3; "
        H_energy += "ANG1=angle(a2,a1,d1); cosB=-2*cos(ANG2)^3; ANG2=angle(d2,d1,a1); "
        
        return H_energy
    
    def find_hlist(self,update_pdb):
        u = mda.Universe(update_pdb)

        list_G = []  ### Take G as donor 
        list_C = []  ### Take C as acceptor
        list_A = []
        list_U = []

        for i in np.unique(u.atoms.chainIDs):
            for j in np.unique(u.select_atoms(f'chainID {i}').resids):

                if 1<j<len(np.unique(u.select_atoms(f'chainID {i}').resids)):

                    sugar = u.select_atoms(f'chainID {i} and (resid {j} and (name Sr))').indices[0]
                    pho_next = u.select_atoms(f'chainID {i} and (resid {j+1} and (name P))').indices[0]
                    base = u.select_atoms(f'chainID {i} and (resid {j} and (name G or name U or name A or name C))')

                    if base.atoms[0].name =='G':
                        list_G.append((base.indices[0],sugar,pho_next,j,i))

                    if base.atoms[0].name =='C':

                        list_C.append((base.indices[0],sugar,pho_next,j,i))

                    if base.atoms[0].name =='A':
                        list_A.append((base.indices[0],sugar,pho_next,j,i))

                    if base.atoms[0].name =='U':
                        list_U.append((base.indices[0],sugar,pho_next,j,i))
                        
        return list_G,list_C,list_A,list_U
    
    def stack_intra(self,pdb_file):
               
        import openmm.unit as u
        Stackingforce = mm.CustomCompoundBondForce(7, "U0/(1.0 + kbond*(distance(p6,p7) - r0)^2  + kphi1*(dihedral(p1,p2,p3,p4) - phi10 + pi*(((pi - dihedral(p1,p2,p3,p4) + phi10)/abs(pi - dihedral(p1,p2,p3,p4) + phi10)) -  ((pi + dihedral(p1,p2,p3,p4) - phi10)/abs(pi + dihedral(p1,p2,p3,p4) - phi10)))  )^2  +  kphi2*(dihedral(p2,p3,p4,p5) - phi20 + pi*(((pi - dihedral(p2,p3,p4,p5) + phi20)/abs(pi - dihedral(p2,p3,p4,p5) + phi20)) -  ((pi + dihedral(p2,p3,p4,p5) - phi20)/abs(pi + dihedral(p2,p3,p4,p5) - phi20)))  )^2 )");
        Stackingforce.addPerBondParameter("U0");
        Stackingforce.addPerBondParameter("r0");
        Stackingforce.addPerBondParameter("phi10");
        Stackingforce.addPerBondParameter("phi20");
        Stackingforce.setUsesPeriodicBoundaryConditions(True)
        Stackingforce.addGlobalParameter('kbond', 1.45/u.angstroms**2)
        Stackingforce.addGlobalParameter('kphi1', 3.0/u.radians**2)
        Stackingforce.addGlobalParameter('kphi2', 3.0/u.radians**2)
        Stackingforce.addGlobalParameter('pi', np.pi)
        Stackingforce.setUsesPeriodicBoundaryConditions(True)
        
        stacking_all = self.find_stack(pdb_file)
        #print(stacking_all)

        for i in stacking_all:
            for j in i:        
                ### ('GG',base.indices[0],base_next.indices[0],pho,sugar,pho_next,s_next,pho_2next)

                r0 = all_ss[j[0]][0]
                U0 = all_ss[j[0]][1]
                phi10 = PSPS
                phi20 = SPSP
                p1 = j[3]
                p2 = j[4]
                p3 = j[5]
                p4 = j[6]
                p5 = j[7]
                p6 = j[1]
                p7 = j[2]        
                group_add = [p1, p2, p3, p4, p5, p6, p7]
                Stackingforce.addBond(group_add, [U0*u.kilocalories_per_mole,
                                                  r0*u.angstroms, phi10*u.radians, 
                                                  phi20*u.radians])

        Stackingforce.setForceGroup(4)          
        self.system.addForce(Stackingforce)
        
    def _stack_bond_energy_expression(self):
        # Long string representing the hydrogen bond energy expression
        return "U0/(1.0 + kbond*(distance(p6,p7) - r0)^2  "\
    "+ kphi1*(dihedral(p1,p2,p3,p4) - phi10 + "\
    "pi*(((pi - dihedral(p1,p2,p3,p4) + phi10)/abs(pi - dihedral(p1,p2,p3,p4) + phi10)) "\
    "-  ((pi + dihedral(p1,p2,p3,p4) - phi10)/abs(pi + dihedral(p1,p2,p3,p4) - phi10)))  )^2"\
    "+  kphi2*(dihedral(p2,p3,p4,p5) - phi20 + "\
    "pi*(((pi - dihedral(p2,p3,p4,p5) + phi20)/abs(pi - dihedral(p2,p3,p4,p5) + phi20)) "\
    "-  ((pi + dihedral(p2,p3,p4,p5) - phi20)/abs(pi + dihedral(p2,p3,p4,p5) - phi20)))  )^2 )";


    def find_stack(self,pdb_file):                      
        u = mda.Universe(pdb_file)
        s_gg = []
        s_gc = []
        s_cc = []
        s_cg = []
        for i in np.unique(u.atoms.chainIDs):
            for j in np.unique(u.select_atoms(f'chainID {i}').resids):

                if j>1 and j<len(np.unique(u.select_atoms(f'chainID {i}').resids))-1:

                    pho = u.select_atoms(f'chainID {i} and (resid {j} and (name P))').indices[0]
                    sugar = u.select_atoms(f'chainID {i} and (resid {j} and (name Sr))').indices[0]
                    pho_next = u.select_atoms(f'chainID {i} and (resid {j+1} and (name P))').indices[0]
                    s_next = u.select_atoms(f'chainID {i} and (resid {j+1} and (name Sr))').indices[0]
                    pho_2next = u.select_atoms(f'chainID {i} and (resid {j+2} and (name P))').indices[0]

                    base = u.select_atoms(f'chainID {i} and (resid {j} and (name G or name U or name A or name C))')
                    base_next = u.select_atoms(f'chainID {i} and (resid {j+1} and (name G or name U or name A or name C))')

                    if base.atoms[0].name =='G'and base_next.atoms[0].name=='G':
                        s_gg.append(('GG',base.indices[0],base_next.indices[0],pho,sugar,pho_next,s_next,pho_2next))

                    if base.atoms[0].name =='G'and base_next.atoms[0].name=='C':
                        s_gc.append(('GC',base.indices[0],base_next.indices[0],pho,sugar,pho_next,s_next,pho_2next))

                    if base.atoms[0].name =='C'and base_next.atoms[0].name=='C':
                        s_cc.append(('CC',base.indices[0],base_next.indices[0],pho,sugar,pho_next,s_next,pho_2next))

                    if base.atoms[0].name =='C'and base_next.atoms[0].name=='G':
                        s_cg.append(('CG',base.indices[0],base_next.indices[0],pho,sugar,pho_next,s_next,pho_2next))

        stacking_all = [s_gg,s_gc,s_cc,s_cg]
        return stacking_all

    def run_nvt(self,output):

#        try:
#            warnings.warn(f"Using GPU platform")
#            platform_name = mm.Platform.getPlatformByName("CUDA")

#        except mm.OpenMMException as e:
#            warnings.warn(f"Using CPU platform")
#            platform_name = mm.Platform.getPlatformByName("CPU")

        temperature = 250 * simtk.openmm.unit.kelvin
        platform_name='OpenCL'


        integrator = simtk.openmm.LangevinIntegrator(temperature,
                                                     10 / simtk.openmm.unit.picosecond,
                                                     15 * simtk.openmm.unit.femtoseconds)
        
        platform = simtk.openmm.Platform.getPlatformByName(platform_name)
        simulation = simtk.openmm.app.Simulation(self.topology,
                                                 self.system,
                                                 integrator,
                                                 platform)

        simulation.context.setPositions(self.coord)
        energy_unit=simtk.openmm.unit.kilojoule_per_mole
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(energy_unit)
        print(energy)

        dcd_reporter=simtk.openmm.app.DCDReporter(f'{output}.dcd', 10000)

        energy_reporter=simtk.openmm.app.StateDataReporter(sys.stdout, 100000, 
                                                           step=True,time=True,
                                                           potentialEnergy=True, 
                                                           temperature=True)
        simulation.reporters.append(dcd_reporter)
        simulation.reporters.append(energy_reporter)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(2000000000)




