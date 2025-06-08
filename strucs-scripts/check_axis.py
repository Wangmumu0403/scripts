import numpy as np
import sys

def calculate_lattice_parameters(matrix):
    """
    Calculates crystallographic lattice parameters (a, b, c, alpha, beta, gamma)
    from a 3x3 lattice vector matrix.
    Assumes rows of the matrix are the lattice vectors a_vec, b_vec, c_vec.
    """
    a_vec = matrix[0]
    b_vec = matrix[1]
    c_vec = matrix[2]

    # Calculate lengths (a, b, c)
    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)

    # Calculate angles (alpha, beta, gamma) in radians first
    # Alpha (angle between b and c)
    cos_alpha = np.dot(b_vec, c_vec) / (b * c)
    alpha_rad = np.arccos(np.clip(cos_alpha, -1.0, 1.0)) # Clip to avoid floating point errors

    # Beta (angle between a and c)
    cos_beta = np.dot(a_vec, c_vec) / (a * c)
    beta_rad = np.arccos(np.clip(cos_beta, -1.0, 1.0))

    # Gamma (angle between a and b)
    cos_gamma = np.dot(a_vec, b_vec) / (a * b)
    gamma_rad = np.arccos(np.clip(cos_gamma, -1.0, 1.0))

    # Convert angles to degrees
    alpha_deg = np.degrees(alpha_rad)
    beta_deg = np.degrees(beta_rad)
    gamma_deg = np.degrees(gamma_rad)

    return a, b, c, alpha_deg, beta_deg, gamma_deg

def read_poscar_lattice_vectors(filepath):
    """
    Reads the 3x3 lattice vector matrix from lines 3-5 of a POSCAR file.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Check if there are enough lines
        if len(lines) < 5:
            raise ValueError("POSCAR file must have at least 5 lines for lattice vectors.")

        # Line 2 is the scaling factor
        scaling_factor = float(lines[1].strip())

        # Lines 3, 4, 5 are the lattice vectors
        matrix = []
        for i in range(2, 5): # Lines are 0-indexed, so 2, 3, 4 correspond to lines 3, 4, 5
            parts = lines[i].strip().split()
            # Ensure 3 numerical values per line
            if len(parts) != 3:
                raise ValueError(f"Line {i+1} in POSCAR does not contain 3 numerical values.")
            vector = [float(p) for p in parts]
            matrix.append(vector)

        # Apply scaling factor to the matrix
        scaled_matrix = np.array(matrix) * scaling_factor
        return scaled_matrix

    except FileNotFoundError:
        print(f"Error: POSCAR file not found at '{filepath}'")
        sys.exit(1)
    except ValueError as e:
        print(f"Error reading POSCAR file: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_axis.py <POSCAR_filepath>")
        sys.exit(1)

    poscar_filepath = sys.argv[1]
    
    # Read the lattice vector matrix from POSCAR
    lattice_matrix = read_poscar_lattice_vectors(poscar_filepath)

    # Calculate the lattice parameters
    a, b, c, alpha, beta, gamma = calculate_lattice_parameters(lattice_matrix)

    # Output the results in the desired format
    print("      a         b         c       alpha     beta     gamma")
    print(f"{a:9.5f} {b:9.5f} {c:9.5f} {alpha:9.4f} {beta:9.4f} {gamma:9.4f}")
