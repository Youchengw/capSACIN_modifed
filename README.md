#  
Thank you for your interest in our project! This repository contains all necessary files and scripts to build capsid surfaces and prepare them for molecular dynamics simulations using the CapSACIN workflow. Below, you'll find an overview of the repository structure and how to get started.

## Citation
If you use this repository in your research, please cite the following articles:

> Jonathan W. P. Zajac, Idris Tohidian, Praveen Muralikrishnan, Caryn L. Heldt, Sarah L. Perry, Sapna Sarupria. "Cracking the Capsid Code: A Computationally-Feasible Approach for Investigating Virus-Excipient Interactions in Biologics Design" (2026) J. Chem. Theory Comput. doi: 10.1021/acs.jctc.5c01810

## Repository Structure

### **Surface Construction and Restraint Generation** (`systemSetup/`)
This directory includes the key scripts required to generate a capsid surface model from an experimentally-resolved capsid structure. 

- **`sliceCapsid.py`** - Requires an input .pdb structure, desired symmetry, and a reference atom from the input .pdb file. Output is a capsid surface abstraction suitable for MD simulations.
- **`genRestraints.py`** - From the sliced configuration, generates position restraints for equilibration simulations. Atoms closer to z=0 are restrained fully, from z=1.5 to z=3.0 nm restraints are partial, and restraints are removed for atoms where z > 3.0.

## Quick Start
1. Clone the repository:
   ```bash
   git clone git@github.com:SAMPEL-Group/CapSACIN.git
   cd capSACIN
   ```
2. Create Conda environment:
   ```bash
   cd env
   source create-env.sh
   ```
3. Start slicing!
   ```bash
   cd ../systemSetup
   python sliceCapsid.py --pdb 1k3v --refindex 959 --symmetry 5 --weight 0.6
   ```
4. Create position restraints:
   ```bash
   python genRestraints.py --pdb 1k3v --weight 0.6
   ```

## License
All written and graphical materials here are made available under a CC-BY 4.0 license, and all source code/software is made available under an MIT license. Both of these allow broad reuse with attribution.

## Contact
For questions or feedback, please contact Jonathan Zajac at zajac028@umn.edu.
