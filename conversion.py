"""
conversion.py

Patrick Wang
patrick.wang@chem.ox.ac.uk


CLI script to convert from PySCF configuration interaction 
vector notation of alpha and beta determinant to Molpro-like
representation: i.e.

PySCF                                 Molpro
alpha     beta     coeff              ci vector
[0 1 2]  [0 1 2]   0.99               222000
[0 1 2]  [0 1 3]   0.03               22ab00

To run:
    python conversion.py input.txt output.csv -o 6 -e 6

    -o number of orbitals in active space
    -e number electrons in active space

"""


import argparse
from typing import List, Tuple
import sys
import csv

class CIConverter:
    def __init__(self, n_orbitals: int, n_electrons: int):
        """
        Initialize converter for CI determinant to occupation notation.
        
        Args:
            n_orbitals: Total number of orbitals in the active space
            n_electrons: Total number of electrons in the active space
        """
        self.n_orbitals = n_orbitals
        self.n_electrons = n_electrons
        
    def det_to_occupation(self, det_alpha: List[int], det_beta: List[int]) -> str:
        """
        Convert alpha/beta determinant representation to occupation notation.
        
        Args:
            det_alpha: List of orbital indices occupied by alpha electrons
            det_beta: List of orbital indices occupied by beta electrons
            
        Returns:
            String in occupation notation (2,a,b,0)
        """
        # Initialize all orbitals as empty
        occupation = ['0'] * self.n_orbitals
        
        # Process alpha electrons
        for orbital_idx in det_alpha:
            if orbital_idx >= self.n_orbitals:
                raise ValueError(f"Orbital index {orbital_idx} exceeds number of orbitals {self.n_orbitals}")
            if occupation[orbital_idx] == '0':
                occupation[orbital_idx] = 'a'
            elif occupation[orbital_idx] == 'b':
                occupation[orbital_idx] = '2'
            else:
                raise ValueError(f"Invalid: orbital {orbital_idx} already has an alpha electron")
                
        # Process beta electrons
        for orbital_idx in det_beta:
            if orbital_idx >= self.n_orbitals:
                raise ValueError(f"Orbital index {orbital_idx} exceeds number of orbitals {self.n_orbitals}")
            if occupation[orbital_idx] == '0':
                occupation[orbital_idx] = 'b'
            elif occupation[orbital_idx] == 'a':
                occupation[orbital_idx] = '2'
            else:
                raise ValueError(f"Invalid: orbital {orbital_idx} already has a beta electron")
        
        # Verify electron count
        n_electrons_found = sum(2 if x == '2' else 1 if x in ['a', 'b'] else 0 for x in occupation)
        if n_electrons_found != self.n_electrons:
            raise ValueError(f"Expected {self.n_electrons} electrons but found {n_electrons_found}")
            
        return ''.join(occupation)

def parse_determinant(det_str: str) -> List[int]:
    """Parse a determinant string into a list of integers."""
    try:
        # Remove brackets and split
        det_str = det_str.strip('[]')
        if not det_str:
            return []
        return [int(x) for x in det_str.split()]
    except ValueError as e:
        raise ValueError(f"Error parsing determinant {det_str}: {e}")

def process_input_file(input_path: str, output_path: str, n_orbitals: int, n_electrons: int):
    """Process the input file and write results to output file."""
    converter = CIConverter(n_orbitals, n_electrons)
    
    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            # Write header
            outfile.write("det-alpha,det-beta,CI coefficient,occupation\n")
            
            # Skip header if present
            first_line = infile.readline()
            if not any(c.isdigit() for c in first_line):
                # Line doesn't contain numbers, assume it's a header
                pass
            else:
                # Line contains data, reset to start
                infile.seek(0)
            
            # Process each line
            for line in infile:
                # Skip empty lines
                if not line.strip():
                    continue
                    
                try:
                    # Split line and clean up whitespace
                    parts = [p.strip() for p in line.split() if p.strip()]
                    if len(parts) != 3:
                        print(f"Warning: Skipping malformed line: {line.strip()}", file=sys.stderr)
                        continue
                    
                    # Parse determinants and coefficient
                    det_alpha = parse_determinant(parts[0])
                    det_beta = parse_determinant(parts[1])
                    coeff = float(parts[2])
                    
                    # Convert to occupation notation
                    occupation = converter.det_to_occupation(det_alpha, det_beta)
                    
                    # Write to output file
                    outfile.write(f"{parts[0]},{parts[1]},{coeff},{occupation}\n")
                    
                except (ValueError, IndexError) as e:
                    print(f"Warning: Error processing line '{line.strip()}': {e}", file=sys.stderr)
                    continue
                    
    except FileNotFoundError:
        print(f"Error: Could not find input file: {input_path}", file=sys.stderr)
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied when accessing files", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Convert CI determinants to occupation notation')
    parser.add_argument('input', help='Input file containing CI determinants')
    parser.add_argument('output', help='Output file path')
    parser.add_argument('--orbitals', '-o', type=int, required=True,
                      help='Number of orbitals in active space')
    parser.add_argument('--electrons', '-e', type=int, required=True,
                      help='Number of electrons in active space')
    
    args = parser.parse_args()
    
    try:
        process_input_file(args.input, args.output, args.orbitals, args.electrons)
        print(f"Successfully converted determinants to occupation notation.")
        print(f"Results written to: {args.output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
