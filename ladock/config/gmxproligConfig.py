config_content = """# Silakan melakukan editing sesuai kebutuhan dan kenyataan Anda
# Nilai untuk -bt (box type):
box_type = triclinic

# Nilai untuk -d (distance)
distance = 1.2

# Nilai untuk -cs (coordinate file for solvent):
coordinate_file = spc216.gro

# This directory structure has been created using LADOCK. Follow these steps to prepare your input files and execute the simulations:

# 1. Place each receptor and ligand(s) (separate files) in the each "complex" directory.
# 2. Make necessary edits to the mdp files in the job directory  as needed.
# 3. Ensure that receptor file in .pdb format.
# 4. Ensure that ligand file(s) in .mol format.
# 5. Adjust the parameters in the 'config_gmxprolig.txt' file in the job directory as needed.
# 6. Verify that the required dependencies, such as acpype, openbabel, and gromacs, are installed and functioning properly.
# 7. Execute the following command: "ladock --run gmxprolig"
"""

