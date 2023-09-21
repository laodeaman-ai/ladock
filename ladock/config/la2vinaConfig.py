config_content = """# Silakan melakukan editing sesuai kebutuhan dan kenyataan Anda
size_x = 40
size_y = 40
size_z = 40
num_modes = 8
exhaustiveness = 8
cpu = 12
num_iterations = 1000
mgl_directory = /home/<user>/MGLTools-1.5.x

# This directory structure has been created using LADOCK. Follow these steps to prepare your input files and execute the simulations:
# 1. Place your test ligand files in the 'ligand_input' directory.
# 2. Place your main ligand files in the 'main_ligand' directory.
# 3. Put your receptor (model_*_rec.pdb) and reference ligand (model_*_lig.pdb) files in 'model_*' directories.
# 4. Ensure that your ligand files (test and main) are in one of the following formats: smi, smiles, mol, mol2, pdb, or sdf.
# 5. For ligands obtained from online sources, add the link to 'ligand_link.txt' and place it in the 'ligand_input' directory.
# 6. Adjust the parameters in the 'config_la2vina.txt' file in the job directory as needed. Make sure to correctly set the MGLTools installation directory and docking parameters in this file.
# 7. Execute the following command in the LADOCK_la2vina directory: `ladock --run la2vina`. 
# """

ligand_link_default = """# Add download links here. One link per line. Supported download link formats include: smi, smiles, mol, mol2, pdb, or sdf.
# Here's an example link (remove this link and add your own ligand file links):
# http://files.docking.org/2D/BA/BAAA.smi
# http://files.docking.org/2D/BA/BAAB.smi
"""
