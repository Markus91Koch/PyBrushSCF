# PyBrushSCF

PyBrushSCF is a Python implementation of the Scheutjens-Fleer self-consitent field method [1,2] for the purpose of creating polymer brushes with various special modifications.

Via the self-consitent solution, it is possible to extract the properties of a polymer brush, such as total monomer density, single monomer distributions (typically for end monomers), and free energies. An attractive interaction between the end monomers of the brush chains and the grafting surface can be activated. Additionally, it is possible to embded a test chain (or minority chain) into a given polymer brush, make them adsorption-attractive, and extract characteristic properties.

## Planned Extensions

- Adsorption-attractive copolymers
- Consistency checks
- Error messages
- Bidispersity or general polydispersity of the brush or of the minority chains (Schulz-Zimm distribution)
- Further performance optimization

## References

This code has been used in the following publications:

- M. Koch, "Adsorption Behavior of Chains in Polymer Brushes", Bachelor Thesis, TU Dresden, 2013
- M. Koch, D. Romeis, J.-U. Sommer, *[Macromolecules, 2020, 53, 17, 7356–7368](https://doi.org/10.1021/acs.macromol.0c01094)*

This code is originally based on a Fortran95 Code developed by my coworker Dirk Romeis.

### Additional Sources

- [1] Scheutjens, J. M. H. M., Fleer, G. J., J. Phys. Chem. 1979, 83, 1619−1635 *[J. Phys. Chem. 1979, 83, 1619−1635](https://dx.doi.org/10.1021/j100475a012)*
- [2] Fleer, G. J., Cohen Stuart, M. A., Scheutjens, J. M. H. M., Cosgrove, T., Vincent, B. "Polymers at Interfaces", Chapman and Hall: London, 1993
