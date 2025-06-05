import MDAnalysis as mda
import parmed
import numpy as np
import pandas as pd
from parmed import Bond, Angle

def process_pdb(input_pdb,no_of_atoms, output_pdb='updated_packed.pdb', output_psf='output.psf'):
    # Load the structure using MDAnalysis and ParmEd
    u = mda.Universe(input_pdb)
    file = parmed.load_file(input_pdb)

    # Create a chain ID column based on atom number, assuming 41 atoms per chain
    df = file.to_dataframe()
    df['chain'] = ((df['number'] - 1) // no_of_atoms) + 1
    
    # Update chain IDs in MDAnalysis
    for atom, chain_id in zip(u.atoms, df['chain']):
        atom.chainID = str(chain_id)

    # Save updated structure to a new PDB file
    u.atoms.write(output_pdb)
    
    # Reload updated structure for angle and bond assignments
    u = mda.Universe(output_pdb)
    
    # Define angle and bond lists
    angle_psg, angle_psp, angle_bsp, angle_sps = [], [], [], []
    bond_sb, bond_sp = [], []
    
    # Process each chain and residue for bond and angle calculation
    for chain_id in np.unique(u.atoms.chainIDs):
        for resid in np.unique(u.select_atoms(f'chainID {chain_id}').resids):
            atoms_in_res = u.select_atoms(f'chainID {chain_id} and resid {resid}')
            atoms_in_next_res = u.select_atoms(f'chainID {chain_id} and resid {resid + 1}')
            
            # Identify indices for each atom type in residue
            sugar = atoms_in_res.select_atoms('name Sr').indices[0]
            base = atoms_in_res.select_atoms('name G or name U or name A or name C').indices[0]
            pho = atoms_in_res.select_atoms('name P').indices[0] if 'P' in atoms_in_res.names else None
            pho_next = atoms_in_next_res.select_atoms('name P').indices[0] if len(atoms_in_next_res) > 0 else None
            sugar_next = atoms_in_next_res.select_atoms('name Sr').indices[0] if len(atoms_in_next_res) > 0 else None
            
            # Append bonds
            bond_sb.append((sugar, base))
            if pho is not None:
                bond_sp.append((pho, sugar))
            if pho_next is not None:
                bond_sp.append((sugar, pho_next))
            
            # Append angles
            if pho and sugar and base:
                angle_psg.append((pho, sugar, base))
            if pho and sugar and pho_next:
                angle_psp.append((pho, sugar, pho_next))
            if sugar and base and pho_next:
                angle_bsp.append((base, sugar, pho_next))
            if sugar and pho_next and sugar_next:
                angle_sps.append((sugar, pho_next, sugar_next))
    
    # Combine bond and angle data
    bond_index = [bond_sb, bond_sp]
    angles = [angle_psg, angle_psp, angle_bsp, angle_sps]
    
    # Add bonds and angles to the ParmEd structure
    for bond_d in bond_index:
        for bond_ind in bond_d:
            file.bonds.append(Bond(file.atoms[bond_ind[0]], file.atoms[bond_ind[1]]))
    
    for ang_d in angles:
        for ind_ang in ang_d:
            file.angles.append(Angle(file.atoms[ind_ang[0]], file.atoms[ind_ang[1]], file.atoms[ind_ang[2]]))
    
    # Save to PSF file
    file.save(output_psf, overwrite=True)
    
    print(f"Processed files saved: {output_pdb} and {output_psf}")

# Usage
###process_pdb('packed.pdb')




