
"""
ci_combi.py
A simple script to generate all possible electronic configurations for a given molecule.
The script generates all possible distributions of alpha and beta electrons across the molecular orbitals
for a given set of irreducible representations (irreps).

Patrick Wang
patrick.wang@chem.ox.ac.uk

10/02/2025

"""



from itertools import combinations, product
from typing import List, Dict, Tuple, Set
from collections import defaultdict



class ConfigurationGenerator:
    def __init__(self, total_electrons: int, frozen_core: int, irrep_specs: Dict[str, int]):
        """
        Initialize the configuration generator.
        
        Args:
            total_electrons: Total number of electrons in the system
            frozen_core: Number of frozen core electrons
            irrep_specs: Dictionary mapping irrep names to their count
                        e.g., {'Ag': 5, 'B3u': 2, 'B2u': 2, ...}
        """
        self.total_electrons = total_electrons
        self.active_electrons = total_electrons - frozen_core
        self.frozen_core = frozen_core
        self.irrep_specs = irrep_specs
        self.n_orbitals = sum(irrep_specs.values())
        
        # Keep track of irrep positions for formatting
        self.irrep_positions = []
        pos = 0
        for irrep, count in irrep_specs.items():
            self.irrep_positions.append((pos, pos + count))
            pos += count
    
    def _format_configuration(self, distribution: List[str]) -> str:
        """Format the distribution with spaces between irrep groups."""
        result = []
        for start, end in self.irrep_positions:
            result.append(''.join(distribution[start:end]))
        return ' '.join(result)
    
    def _generate_spin_distributions(self) -> Set[str]:
        """Generate possible distributions of alpha and beta electrons."""
        n_alpha = self.active_electrons // 2
        n_beta = self.active_electrons - n_alpha
        distributions = set()
        
        # Generate all possible alpha electron positions
        for alpha_pos in combinations(range(self.n_orbitals), n_alpha):
            # Then, for each alpha distribution, generate beta positions
            for beta_pos in combinations(range(self.n_orbitals), n_beta):
                # Create the orbital string
                dist = ['0'] * self.n_orbitals
                
                # Place alpha electrons
                for pos in alpha_pos:
                    dist[pos] = 'a'
                
                # Place beta electrons
                valid_config = True
                for pos in beta_pos:
                    if dist[pos] == '0':
                        dist[pos] = 'b'
                    elif dist[pos] == 'a':
                        dist[pos] = '2'  # Doubly occupied orbital
                    else:
                        valid_config = False
                        break
                
                if valid_config:
                    # Format with proper spacing between irreps
                    config = self._format_configuration(dist)
                    distributions.add(config)
        
        return distributions

    def generate_configurations(self) -> Set[str]:
        """Generate all possible electronic configurations."""
        return self._generate_spin_distributions()

def main():
    # Example usage for Ne system
    irrep_specs = {
        'Ag': 5,    # xxxx
        'B3u': 2,   # xx
        'B2u': 2,   # xx
        'B1g': 1,   # x
        'B1u': 2,   # xx
        'B2g': 1,   # x
        'B3g': 1    # x
    }
    
    generator = ConfigurationGenerator(
        total_electrons=10,
        frozen_core=2,
        irrep_specs=irrep_specs
    )
    
    configurations = generator.generate_configurations()
    
    print(f"Generated {len(configurations)} unique configurations:")
    print("\nFirst few configurations:")
    for config in sorted(configurations)[:10]:  # Show first 10 configurations
        print(f"{config}")
    with open('configurations.txt', 'w') as f:
        for config in configurations:
            f.write(f"{config}\n")

if __name__ == "__main__":
    main()