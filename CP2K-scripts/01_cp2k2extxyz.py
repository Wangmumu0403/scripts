#!/usr/bin/env python3
# python cp2k2xyz.py [pos.xyz] [frc.xyz] [cell.cell] [stress.stress] [-shifted yes/no]
# -shifted yes, energy shift; -shifted no, not energy shift; The default is "-shifted no".
# If only entering "python cp2k2xyz.py", read according to the default file name; Otherwise, read according to the input file name.
# Merge the box information, atomic coordinates, atomic forces, and stress tensor outputted by CP2K AIMD into the xyz file
# Shift total energy to ~0 by substracting atomic energies using the least square method.

import numpy as np
import sys
import os
import glob

def SVD_A(A, b):
    """Solve A*x=b using Singular Value Decomposition."""
    U, S, V = np.linalg.svd(A)
    # Handle singular values that are very close to zero for numerical stability
    s_inv = np.zeros_like(S)
    # A small tolerance to avoid division by zero or very small numbers
    tol = np.finfo(S.dtype).eps * S.max()
    s_inv[S > tol] = 1 / S[S > tol]
    
    B = np.matmul(U.T, b)
    # Use element-wise multiplication with s_inv to handle potentially different dimensions
    X = B[:len(S)] * s_inv.reshape(-1, 1) 
    x = np.matmul(V.T, X)
    return x

def extract_xyz_data(filename):
    """Extract data from XYZ file, skipping comment lines."""
    frames_data = []
    with open(filename, 'r') as file:
        while True:
            num_atoms_line = file.readline()
            if not num_atoms_line: # End of file
                break
            
            # Skip comment lines or empty lines
            if num_atoms_line.startswith('#') or not num_atoms_line.strip():
                continue
            
            try:
                num_atoms = int(num_atoms_line.strip().split()[0])
            except ValueError:
                continue # Skip non-integer lines
            
            file.readline() # Read and discard the comment line (e.g., energy=...)
            
            atoms_info = []
            for _ in range(num_atoms):
                atom_line = file.readline().strip().split()
                if not atom_line: # Handle unexpected empty lines in atom data
                    continue
                atoms_info.append(atom_line)
            frames_data.append((num_atoms, atoms_info))
    return frames_data # Return as a list directly

def extract_forces_and_energy(frc_file, num_atoms_list):
    """Extract forces and energy from a force file."""
    energies = []
    forces_list = []
    with open(frc_file, 'r') as ff:
        current_frame_idx = 0 # To match with num_atoms_list for frame consistency
        while True:
            line = ff.readline()
            if not line:
                break
            line = line.strip()
            if "E =" in line:
                try:
                    energy = float(line.split("E =")[-1]) * 27.211386245988  # Convert Hartree to eV
                    energies.append(energy)
                except ValueError:
                    continue # Skip if energy parsing fails
                
                # Ensure we have a corresponding num_atoms for the current energy frame
                if current_frame_idx >= len(num_atoms_list):
                    # This indicates an issue where frc file has more frames than pos file
                    print(f"Warning: More energy/force blocks in {frc_file} than atom counts in pos file. Stopping extraction for forces/energies.")
                    break 
                
                num_atoms_current_frame = num_atoms_list[current_frame_idx]
                forces = []
                for _ in range(num_atoms_current_frame):
                    frc_line = ff.readline().strip().split()
                    if len(frc_line) < 4 or not frc_line[1].replace(".", "", 1).lstrip('-').isdigit():
                        continue # Skip malformed lines or lines without numerical force data
                    try:
                        # Convert forces from Hartree/Bohr to eV/Angstrom
                        force_x = float(frc_line[1]) * 51.42206747632590000
                        force_y = float(frc_line[2]) * 51.42206747632590000
                        force_z = float(frc_line[3]) * 51.42206747632590000
                        forces.append((force_x, force_y, force_z))
                    except ValueError:
                        continue
                forces_list.append(forces)
                current_frame_idx += 1
    return energies, forces_list

def extract_cell_data(cell_file):
    """Extract lattice information from cell file."""
    lattices = []
    with open(cell_file, 'r') as cf:
        cf.readline()  # Skip the header line in the cell file
        while True:
            line = cf.readline().strip()
            if not line:
                break
            cell_line = line.split()
            # Read 9 elements for the 3x3 lattice matrix (Ax Ay Az Bx By Bz Cx Cy Cz)
            # Assuming format: <step_id> <time> Ax Ay Az Bx By Bz Cx Cy Cz
            # So, the 9 values start from index 2
            if len(cell_line) >= 11: # Check if enough columns for time + 9 lattice elements
                lattice = " ".join(cell_line[2:11])
                lattices.append(lattice)
            else:
                print(f"Warning: Skipping malformed cell line: {line} (expected at least 11 columns).")
    return lattices

def extract_stress_data(stress_file):
    """
    Extract stress tensor from stress file.
    Assumes file has a header line, and then data lines with 9 stress components
    from column 3 to 11 (0-indexed: [2:11]), in 'bar' units.
    Converts 'bar' to 'eV/A^3' using the factor 6.2415e-7.
    """
    stresses_list = []
    bar_to_eV_A3 = -6.2415e-7 # 1 bar = 6.2415e-7 eV/A^3

    with open(stress_file, 'r') as sf:
        sf.readline() # Skip the first line (header)
        for line in sf:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            try:
                # Read 9 components from column 3 to 11 (Python indices 2 to 10)
                # Assuming format: [idx] [time] Sxx Sxy Sxz Syx Syy Syz Szx Szy Szz
                # Or similar 9 components.
                if len(parts) < 11: # Need at least 2 columns + 9 stress components = 11 columns
                    print(f"Warning: Skipping malformed stress line: {line} (expected at least 11 columns for 9 stress values).")
                    continue
                
                # Extract the 9 values and convert to float
                stress_components_bar = [float(p) for p in parts[2:11]]
                
                # Convert from bar to eV/A^3
                stress_components_eV_A3 = [s * bar_to_eV_A3 for s in stress_components_bar]
                
                stresses_list.append(stress_components_eV_A3)
            except (ValueError, IndexError) as e:
                print(f"Warning: Failed to parse stress line: {line}. Error: {e}. Skipping.")
                continue
    return stresses_list


def write_xyz(output_file, pos_data, forces_list, energies, lattices, stresses_list=None):
    with open(output_file, 'w') as of:
        for idx, (num_atoms, atoms_info) in enumerate(pos_data):
            # Ensure all lists have data for the current index
            # This check is crucial after data trimming in main script
            if idx >= len(energies) or idx >= len(forces_list) or idx >= len(lattices):
                print(f"Error: Data mismatch at index {idx} in write_xyz. Skipping frame.")
                break # Exit loop if any primary data source runs out
            
            energy = energies[idx]
            forces = forces_list[idx]
            lattice = lattices[idx]
            
            # Prepare stress string if available
            stress_str = ""
            properties_str = "species:S:1:pos:R:3:force:R:3" # Default properties string
            
            if stresses_list is not None and idx < len(stresses_list):
                stress_components = stresses_list[idx]
                # Format the 9 components for the output string
                stress_str = f" stress=\"{' '.join(f'{s:.10f}' for s in stress_components)}\""
                properties_str += ":stress:R:9" # Add stress property for 9 components
            
            of.write(f"{num_atoms}\n")
            of.write(f"energy={energy:.10f} config_type=cp2k2xyz pbc=\"T T T\" Lattice=\"{lattice}\" {stress_str}Properties={properties_str}\n")
            
            for i, atom_info in enumerate(atoms_info):
                if i >= len(forces): # Ensure force data exists for this atom
                    print(f"Warning: Force data missing for atom {i} in frame {idx}. Setting to 0s.")
                    force_x, force_y, force_z = 0.0, 0.0, 0.0
                else:
                    force_x, force_y, force_z = forces[i]

                symbol = atom_info[0]
                x, y, z = map(float, atom_info[1:4])
                
                of.write(f"{symbol:<2} {x:>20.10f} {y:>20.10f} {z:>20.10f} {force_x:>20.10f} {force_y:>20.10f} {force_z:>20.10f}\n")

def find_file(pattern):
    """Find a file with a given pattern. Raise error if not exactly one found."""
    files = glob.glob(pattern)
    if len(files) == 0:
        raise FileNotFoundError(f"No file found matching pattern: {pattern}")
    elif len(files) > 1:
        # If multiple, try to pick one, otherwise raise error.
        # This part might need further refinement based on specific naming conventions.
        print(f"Warning: Multiple files found matching pattern '{pattern}'. Found: {files}.")
        # Attempt to pick the one that looks like a main output if specific enough
        if 'pos-1' in pattern:
            selected = next((f for f in files if '-pos-1.xyz' in f), None)
            if selected: return selected
        if 'frc-1' in pattern:
            selected = next((f for f in files if '-frc-1.xyz' in f), None)
            if selected: return selected
        # If still ambiguous, raise error.
        raise SystemError(f"Ambiguous file selection for pattern: {pattern}. Please specify files explicitly.")
    return files[0]

# --- Main Script ---
if __name__ == "__main__":
    # Parse command line arguments
    args = sys.argv[1:]  # Skip the script name
    shifted = "no"  # Default behavior is no energy shifting
    
    # Check for the -shifted argument and remove it and its value from args
    if "-shifted" in args:
        try:
            shifted_index = args.index("-shifted")
            if shifted_index + 1 < len(args):
                shifted = args[shifted_index + 1].lower()
                del args[shifted_index : shifted_index + 2] 
            else:
                print("Warning: '-shifted' argument found without a value ('yes' or 'no'). Defaulting to 'no'.")
                del args[shifted_index] 
        except ValueError: # Should not happen if "-shifted" is in args
            pass

    # Determine file names based on remaining arguments or default patterns
    pos_file, frc_file, cell_file, stress_file = None, None, None, None

    if len(args) == 4: # python cp2k2xyz.py [pos.xyz] [frc.xyz] [cell.cell] [stress.stress]
        pos_file, frc_file, cell_file, stress_file = args
    elif len(args) == 3: # python cp2k2xyz.py [pos.xyz] [frc.xyz] [cell.cell] (no stress)
        pos_file, frc_file, cell_file = args
    elif len(args) == 0: # Default file names
        try:
            pos_file = find_file("*-pos-1.xyz")
            frc_file = find_file("*-frc-1.xyz")
            cell_file = find_file("*.cell")
            # Attempt to find a default stress file, but it's optional
            try:
                stress_file = find_file("*.stress")
            except (FileNotFoundError, SystemError):
                print("Info: No default '*.stress' file found. Stress information will not be included.")
                stress_file = None 
        except (FileNotFoundError, SystemError) as e:
            print(f"Error: Could not find all required default input files: {e}")
            sys.exit(1)
    else:
        print("Usage: python cp2k2xyz.py [pos.xyz] [frc.xyz] [cell.cell] [stress.stress] [-shifted yes/no]")
        print("   Or: python cp2k2xyz.py                                                 (for default files)")
        sys.exit(1)

    # Check if selected files exist
    required_files = [pos_file, frc_file, cell_file]
    # Only check stress_file existence if it was actually provided/found
    if stress_file: 
        required_files.append(stress_file)

    for f in required_files:
        if not os.path.exists(f):
            print(f"Error: Input file '{f}' does not exist.")
            sys.exit(1)

    # Read data from files
    pos_data = list(extract_xyz_data(pos_file))
    num_atoms_list = [num_atoms for num_atoms, _ in pos_data] # Collect number of atoms for each frame
    
    energies, forces_list = extract_forces_and_energy(frc_file, num_atoms_list)
    lattices = extract_cell_data(cell_file)
    
    stresses_list = None
    if stress_file:
        stresses_list = extract_stress_data(stress_file)

    # --- Data Length Consistency Check ---
    # Find the minimum length among all collected data lists
    min_len = len(pos_data) # Start with pos_data length
    min_len = min(min_len, len(energies))
    min_len = min(min_len, len(forces_list))
    min_len = min(min_len, len(lattices))
    if stresses_list: # Only if stress data was actually loaded
        min_len = min(min_len, len(stresses_list))

    # If any list is longer than min_len, trim all lists to min_len
    if len(pos_data) > min_len or len(energies) > min_len or \
       len(forces_list) > min_len or len(lattices) > min_len or \
       (stresses_list and len(stresses_list) > min_len):
        print(f"Warning: Mismatch in number of frames among input files. Trimming all data to {min_len} frames.")
        pos_data = pos_data[:min_len]
        energies = energies[:min_len]
        forces_list = forces_list[:min_len]
        lattices = lattices[:min_len]
        if stresses_list:
            stresses_list = stresses_list[:min_len]
        num_atoms_list = num_atoms_list[:min_len] # Also trim num_atoms_list

    print(f"Processing {min_len} frames.")

    # --- Energy Shifting Logic ---
    if shifted == "yes":
        print('Normalizing energy....')
        # Determine unique elements
        all_elements = sorted(list(set(atom_info[0] for _, atoms_frame in pos_data for atom_info in atoms_frame)))

        coeff_matrix = np.zeros((len(pos_data), len(all_elements)))
        energy_matrix = np.zeros((len(pos_data), 1))

        for idx, (num_atoms, atoms_info) in enumerate(pos_data):
            for j, element in enumerate(all_elements):
                coeff_matrix[idx][j] = sum(1 for atom in atoms_info if atom[0] == element)
            energy_matrix[idx][0] = energies[idx]

        # Check for underdetermined system and add constraints
        if np.linalg.matrix_rank(coeff_matrix) < len(all_elements):
            print("Warning! The coeff_matrix is underdetermined, adding heuristic constraints (setting adjacent element energy differences to zero).")
            if len(all_elements) > 1:
                # Add N-1 constraints for N elements
                for i in range(len(all_elements) - 1):
                    additional_matrix = np.zeros(len(all_elements))
                    additional_matrix[i], additional_matrix[i+1] = 1, -1 # E_i - E_{i+1} = 0
                    additional_energy = np.zeros(1) # We assume difference is zero
                    coeff_matrix = np.r_[coeff_matrix, [additional_matrix]]
                    energy_matrix = np.r_[energy_matrix, [additional_energy]]
            else:
                print("Cannot add constraints for a single element system.")


        # Solve for atomic energy shifts
        atomic_shifted_energy = SVD_A(coeff_matrix, energy_matrix)
        for i, element in enumerate(all_elements):
            print(f"{element}:{atomic_shifted_energy[i][0]:.10f} eV")

        # Calculate shifted energies
        shifted_energies_array = (energy_matrix - np.matmul(coeff_matrix, atomic_shifted_energy)).flatten()
        print("Averaged shifted energy now: %f eV." % shifted_energies_array.mean())
        print("Absolute maximum shifted energy now: %f eV." % np.max(np.abs(shifted_energies_array)))

        # Write shifted XYZ file
        write_xyz("shifted.xyz", pos_data, forces_list, shifted_energies_array.tolist(), lattices, stresses_list)
        print("Done! 'shifted.xyz' file is generated.")

    # Always write the original XYZ file (or the primary output if no shift)
    write_xyz("original-stress.extxyz", pos_data, forces_list, energies, lattices, stresses_list)
    print("Done! 'original-stress.xyz' file is generated.")
