# Bader integration using weighted integration
This code takes the grid-based output from [VASP](https://www.vasp.at/) as from CHGCAR or AECCAR files and performs *Bader integration* over the basins of attraction. The call

    weight_int CHGCAR [grid1] [grid2] ...

uses the first grid to define atom-centered basins of attraction, and integrates each grid-defined quantity over those same basins of attraction. At a minimum, one grid file can be used. The grids must all be *compatible*: same grid dimensions defined for each. It does not check that the atomic coordinates are the same, and it *only* reads in the first grid found in the file; thus, it requires some hacking of the files in order to, e.g., integrate the magnetization. There are a few options available:

* **-a** assigns volumes to true basins (maxima) rather than explicitly to atoms.
* **-s** divides (scales) the integral by the total volume of the cell; required to compute a Bader charge.
* **-o base** output the weights on a grid, to be used for visualization; "base" is the name of the weight file: base0001
* **-n N** output only the weight for atom/basin N
* **-V** use Voronoi volumes instead of Bader (basins of attraction)
* **-v/-t** verbose or testing (extra verbose) modes

# References
Original algorithm described in the first reference; the latter two references include examples. *If you use the code, please cite the first reference* as well as [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18776.svg)](http://dx.doi.org/10.5281/zenodo.18776)

* M. Yu and D. R. Trinkle, "Accurate and efficient algorithm for Bader charge integration." *J. Chem. Phys.* **134**, 064111 (2011). [doi](http://dx.doi.org//10.1063/1.3553716)
* M. Yu, D. R. Trinkle, and R. M. Martin, "Energy density in density functional theory: Application to crystalline defects and surfaces." *Phys. Rev. B* **83**, 115113 (2011). [doi](http://dx.doi.org/10.1103/PhysRevB.83.115113)
* M. Yu and D. R. Trinkle, "Au/TiO2(110) interfacial reconstruction stability from ab initio." *J. Phys. Chem. C* **115**, 17799-17805 (2011). [doi](http://dx.doi.org/10.1021/jp2017133)

# Contributors
Min Yu and Dallas R. Trinkle, algorithm development and implementation

# Acknowledgments
The research was supported by NSF under grant number DMR-1006077 and through the Materials Computation Center at UIUC, NSF DMR-0325939, and with computational resources from NSF/TeraGrid provided by NCSA and TACC. We also thank G. Henkelman at U. Texas for helpful discussions, and  R. E. L. Deville at UIUC for helpful discussions.
