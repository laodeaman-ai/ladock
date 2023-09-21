# LADOCK

LADOCK is an innovative and free tool designed for conducting simultaneous simulations in computer-aided drug discovery, encompassing molecular docking and molecular dynamics. In molecular docking, LADOCK excels in handling single or double ligands. It supports ligands from various online sources. In molecular dynamics, LADOCK efficiently manages protein-ligand interactions and even ligand-ligand interactions, accommodating scenarios with one or multiple proteins and ligands.

## Overview

LADOCK provides a flexible environment for conducting molecular docking and dynamics simulations. It supports several simulation types and offers a user-friendly interface for creating input files and running simulations.

## Installation
A. From GitHub

1. Clone the LADOCK repository to your local machine:
   ```
   git clone https://github.com/laodeaman.ai/LADOCK.git
   ```
2. From your job directory naivigate to the main.py file in the LADOCK directory, and run the following command to create the necessary input files for your desired simulation type (replace `simulation_type` with your choice):
   ```
      python /path/to/ladock/main.py --create <simulation_type>
   ```
5. Edit the input files in the newly created simulation directory to configure your simulation parameters.
6. Run the simulation using the following command (replace `simulation_type` with your choice):
   ```
      python /path/to/ladoc/main.py --run <simulation_type>
   ```
B. From PyPi
   1. Install the LADOCK package
   ```
      pip install ladock
   ```
   2. From your job directory, run the following command to create the necessary input files for your desired simulation type (replace `simulation_type` with your choice):
   ```
      ladock --create <simulation_type>
  ```
   3. Edit the input files in the newly created simulation directory to configure your simulation parameters.
   4. Run the simulation using the following command (replace `simulation_type` with your choice):
   ```
      ladock --run <simulation_type>
   ```

## Simulation Types
- `lavina`: Molecular docking using AutoDock Vina.
- `lavinagpu`: Molecular docking using AutoDock Vina GPU.
- `la2vina`: Molecular docking using AutoDock Vina with multiple ligands simultaneously.
- `ladock4`: Molecular docking using AutoDock4.
- `ladockgpu`: Molecular docking using AutoDock-GPU.
- `gmxprolig`: Molecular dynamics for single or multiple (chain) of protein with single or multiple ligands using GROMACS.
- `gmxliglig`: Molecular dynamics for multiple ligands using GROMACS.

## Developer Note

For detailed information on LADOCK, developer notes, and contact information, please refer to the `developerNote.txt` file in the source code.

## Dependencies

Make sure to install the following dependency packages on your system:
- autodock4
- autodock-gpu
- autodock vina
- vina-gpu
- mgltools
- gromacs
- acpype
- openbabel

## Usage

To see the available options and how to use LADOCK, run the following command:

```bash
./ladock.py --help
```

## License

This software is distributed under the MIT License. See the `LICENSE` file for details.

## Developer Contact

For questions or feedback, please contact the developer at laodeaman.ai@gmail.com, or laode_aman@ung.ac.id.

## Citation

If you use LADOCK in your research, please cite the relevant references provided in the `citation_list` in the source code.

