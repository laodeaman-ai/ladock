config_content = """# Silakan melakukan editing sesuai kebutuhan dan kenyataan Anda
autodockgpu = /home/arga/AutoDock-GPU-develop/bin/autodock_gpu_128wi
mgl_directory = /home/<user>/MGLTools-1.5.x

# This directory structure has been created using LADOCK.
# Follow these steps to prepare your input files and execute the simulations:
# 1. Put your ligand files in the 'ligand_input' directory.
# 2. Put your receptor (*rec.pdb) dan reference ligand (*lig.pdb) files in 'model_*' directories.
# 3. Ensure your ligand input files are in formats: smi, smiles, mol, mol2, pdb, or sdf.
# 4. For ligands from online sources, add the link to 'ligand_link.txt' and put it in the 'ligand_input' directory.
# 5. Adjust the parameters in the 'ladock4gpu_config.txt' file in the job directory as needed.
# 6. Adjust the parameters in the 'config_ladockgpu.txt' file as needed. Set the MGLTools installation directory and docking parameters in this file correctly.
"""

ligand_link_default = """# Add download links here. One link per line. Supported download link formats include: smi, smiles, mol, mol2, pdb, or sdf.
# Here's an example link (remove this link and add your own ligand file links):
# http://files.docking.org/2D/BA/BAAA.smi
# http://files.docking.org/2D/BA/BAAB.smi
"""
