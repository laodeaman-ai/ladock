#!/usr/bin/env python3

import os
import sys
import argparse
import glob
import shutil
import subprocess

current_directory = os.getcwd()
source_directory = os.path.dirname(os.path.abspath(__file__))
config_directory = os.path.join(source_directory, "config")
share_directory = os.path.join(source_directory, "share")
mdp_directory = os.path.join(source_directory, "share", "mdp")
vina = os.path.join(share_directory, "vina_1.2.5_linux_x86_64")
vina_split = os.path.join(share_directory, "vina_split_1.2.5_linux_x86_64")

dn = os.path.join(source_directory, "developerNote.txt")
with open(dn, 'r') as file:
    content = file.read()

# Pisahkan isi file menjadi blok-blok yang sesuai dengan variabel-variabel
blocks = content.split('\n\n')
# Inisialisasi variabel
developer_note = ""
developer_contact = ""
citation_list = ""

# Loop melalui blok-blok dan mengisi variabel yang sesuai
for block in blocks:
    if block.startswith("developer_note ="):
        developer_note = block[block.find("=") + 1:].strip()
    elif block.startswith("developer_contact ="):
        developer_contact = block[block.find("=") + 1:].strip()
    elif block.startswith("citation_list ="):
        citation_list = block[block.find("=") + 1:].strip()

def create_lavina():
    # Membuat struktur direktori dan file konfigurasi di dalamnya
    if not os.path.exists('LADOCK_lavina'):
        os.mkdir('LADOCK_lavina')
    os.chdir('LADOCK_lavina')

    if not os.path.exists('ligand_input'):
        os.mkdir('ligand_input')
   
    for i in range(1, 4):  # Membuat direktori model_01, model_02, model_03, model_04
        model_dir = f'model_{i:02d}'
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

    # Mengambil isi 'config_content' dari file 'lavinaConfig.py'
    from config.lavinaConfig import config_content
    with open('config_lavina.txt', 'w') as config_file:
        config_file.write(config_content)

    # Mengambil isi 'ligand_link_default' dari file 'lavinaConfig.py'
    from config.lavinaConfig import ligand_link_default
    ligand_input_dir = os.path.join('ligand_input', 'ligand_link.txt')   
    with open(ligand_input_dir, 'w') as ligand_link_file:
        ligand_link_file.write(ligand_link_default) 
 
def create_lavinagpu():
    # Membuat struktur direktori dan file konfigurasi di dalamnya
    if not os.path.exists('LADOCK_lavinagpu'):
        os.mkdir('LADOCK_lavinagpu')
    os.chdir('LADOCK_lavinagpu')

    if not os.path.exists('ligand_input'):
        os.mkdir('ligand_input')
   
    for i in range(1, 4):  # Membuat direktori model_01, model_02, model_03, model_04
        model_dir = f'model_{i:02d}'
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

    # Mengambil isi 'config_content' dari file 'lavinagpuConfig.py'
    from config.lavinagpuConfig import config_content
    with open('config_lavinagpu.txt', 'w') as config_file:
        config_file.write(config_content)

    # Mengambil isi 'ligand_link_default' dari file 'lavinaConfig.py'
    from config.lavinagpuConfig import ligand_link_default
    ligand_input_dir = os.path.join('ligand_input', 'ligand_link.txt')   
    with open(ligand_input_dir, 'w') as ligand_link_file:
        ligand_link_file.write(ligand_link_default) 
        
def create_la2vina():
    # Membuat struktur direktori dan file konfigurasi di dalamnya
    if not os.path.exists('LADOCK_la2vina'):
        os.mkdir('LADOCK_la2vina')
    os.chdir('LADOCK_la2vina')

    if not os.path.exists('ligand_input'):
        os.mkdir('ligand_input')

    if not os.path.exists('main_ligand'):
        os.mkdir('main_ligand')
    
    for i in range(1, 4):  # Membuat direktori model_01, model_02, model_03, model_04
        model_dir = f'model_{i:02d}'
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

    # Mengambil isi 'config_content' dari file 'la2vinaConfig.py'
    from config.la2vinaConfig import config_content
    with open('config_la2vina.txt', 'w') as config_file:
        config_file.write(config_content)

    # Mengambil isi 'ligand_link_default' dari file 'lavinaConfig.py'
    from config.lavinaConfig import ligand_link_default
    
    ligand_input_dir = os.path.join('ligand_input', 'ligand_link.txt')   
    with open(ligand_input_dir, 'w') as ligand_link_file:
        ligand_link_file.write(ligand_link_default)
       
def create_ladock4():
    # Membuat struktur direktori dan file konfigurasi di dalamnya
    if not os.path.exists('LADOCK_ladock4'):
        os.mkdir('LADOCK_ladock4')
    os.chdir('LADOCK_ladock4')

    if not os.path.exists('ligand_input'):
        os.mkdir('ligand_input')
    
    for i in range(1, 4):  # Membuat direktori model_01, model_02, model_03, model_04
        model_dir = f'model_{i:02d}'
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

    # Mengambil isi 'config_content' dari file 'ladock4Config.py'
    from config.ladock4Config import config_content
    with open('config_ladock4.txt', 'w') as config_file:
        config_file.write(config_content)

    # Mengambil isi 'ligand_link_default' dari file 'ladock4Config.py'
    from config.ladock4Config import ligand_link_default
    ligand_input_dir = os.path.join('ligand_input', 'ligand_link.txt')   
    
    with open(ligand_input_dir, 'w') as ligand_link_file:
        ligand_link_file.write(ligand_link_default)
        
def create_ladockgpu():
    # Membuat struktur direktori dan file konfigurasi di dalamnya
    if not os.path.exists('LADOCK_ladockGPU'):
        os.mkdir('LADOCK_ladockGPU')
    os.chdir('LADOCK_ladockGPU')

    if not os.path.exists('ligand_input'):
        os.mkdir('ligand_input')
    
    for i in range(1, 4):  # Membuat direktori model_01, model_02, model_03, model_04
        model_dir = f'model_{i:02d}'
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

    # Mengambil isi 'config_content' dari file 'ladock4Config.py'
    from config.ladockgpuConfig import config_content
    with open('config_ladockgpu.txt', 'w') as config_file:
        config_file.write(config_content)

    # Mengambil isi 'ligand_link_default' dari file 'ladockGPUConfig.py'
    from config.ladockgpuConfig import ligand_link_default
    ligand_input_dir = os.path.join('ligand_input', 'ligand_link.txt')   
    
    with open(ligand_input_dir, 'w') as ligand_link_file:
        ligand_link_file.write(ligand_link_default)      

def create_gmxprolig():
    # Membuat struktur direktori dan file konfigurasi di dalamnya
    if not os.path.exists('LADOCK_gmxprolig'):
        os.mkdir('LADOCK_gmxprolig')
    os.chdir('LADOCK_gmxprolig')
   
    for i in range(1, 4):  # Membuat direktori model_01, model_02, model_03, model_04
        model_dir = f'complex_{i:02d}'
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

    # Menyalin file *.mdp dari direktori source ke direktori tujuan
    source_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "share", "mdp")
    destination_dir = os.getcwd()

    if not os.path.exists(destination_dir):
        os.mkdir(destination_dir)

    for filename in os.listdir(source_dir):
        if filename.endswith(".mdp"):
            source_file = os.path.join(source_dir, filename)
            destination_file = os.path.join(destination_dir, filename)
            shutil.copy(source_file, destination_file)

    # Mengambil isi 'config_content' dari file 'lagmxConfig.py'
    from config.gmxproligConfig import config_content
    with open('config_gmxprolig.txt', 'w') as config_file:
        config_file.write(config_content)
        
def create_gmxliglig():
    # Membuat struktur direktori dan file konfigurasi di dalamnya
    if not os.path.exists('LADOCK_gmxliglig'):
        os.mkdir('LADOCK_gmxliglig')
    os.chdir('LADOCK_gmxliglig')
    
    for i in range(1, 4):  # Membuat direktori model_01, model_02, model_03, model_04
        model_dir = f'complex_{i:02d}'
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)
            
    # Menyalin file *.mdp dari direktori source ke direktori tujuan
    source_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "share", "mdp")
    destination_dir = os.getcwd()

    if not os.path.exists(destination_dir):
        os.mkdir(destination_dir)

    for filename in os.listdir(source_dir):
        if filename.endswith(".mdp"):
            source_file = os.path.join(source_dir, filename)
            destination_file = os.path.join(destination_dir, filename)
            shutil.copy(source_file, destination_file)

    # Mengambil isi 'config_content' dari file 'ligligConfig.py'
    from config.ladyna3ligConfig import config_content
    with open('config_gmxliglig.txt', 'w') as config_file:
        config_file.write(config_content)

def main():
    parser = argparse.ArgumentParser(description='LADOCK is an innovative and free tool designed for conducting simultaneous simulations in computer-aided drug discovery, encompassing molecular docking and molecular dynamics. In molecular docking, LADOCK excels in handling single or double ligands. It supports ligands from various online sources. In molecular dynamics, LADOCK efficiently manages protein-ligand interactions and even ligand-ligand interactions, accommodating scenarios with one or multiple proteins and ligands.')
    parser.add_argument('--version', action='version', version='ladock 1.0')
    parser.add_argument('--create', choices=['lavina', 'lavinagpu', 'la2vina', 'ladock4', 'ladockgpu', 'gmxprolig', 'gmxliglig'], help='Create necessary input')
    parser.add_argument('--run', choices=['lavina', 'lavinagpu', 'la2vina', 'ladock4', 'ladockgpu', 'gmxprolig', 'gmxliglig'], help='Execute the specific simulation: lavina (AutoDock Vina), la2vina (MLSD with Autodock Vina), ladock4 (AutoDock4), ladockgpu (AutoDock-GPU), gmxprolig (gromacs for 1 protein with 1 or multiple ligands), gmxliglig (gromacs multiple ligands)' )
    args = parser.parse_args()

    if args.create:
        create_input(args.create)

    elif args.run:
        run_simulation(args.run)
    
    else:
        print("Invalid option. Please use one of the following:")
        print("'--create lavina' to create input for AutoDock Vina simulation.")
        print("'--create lavinagpu' to create input for AutoDock Vina GPU simulation.")
        print("'--create la2vina' to create input for AutoDock Vina simulation using multiple ligand simultaneously molecular-docking technique).")
        print("'--create ladock4' to create input for AutoDock4 simulation.")
        print("'--create ladockgpu' to create input for AutoDock4-GPU simulation.")
        print("'--create gmxprolig' to create input for single or multiple ligand simultaneously molecular dynamics simulation with GROMACS.")
        print("'--create gmxliglig' to create input for single/multiple ligand and single/multiple (chain) protein simultaneously molecular dynamics simulation with GROMACS")
        print("'--run lavina' for docking simulation with AutoDock Vina.")
        print("'--run la2vina' for docking simulation with AutoDock Vina using MLSD (multiple ligand simulation docking technique).")
        print("'--run ladock4' for docking simulation with AutoDock4.")
        print("'--run ladockgpu' for docking simulation with AutoDock-GPU.")
        print("'--run gmxprolig' for single or multiple ligand simultaneously molecular dynamics simulation with GROMACS.")
        print("'--run gmxliglig' for single/multiple ligand and single/multiple (chain) protein simultaneously molecular dynamics simulation with GROMACS")
        print("'--version' to show the LADOCK's version number.")
        print("'--help' or '-h' for help.")
        print("\nNote: Ensure that the following dependency packages are installed on your system:")
        print("* autodock4")
        print("* autodock-gpu")
        print("* autodock vina")
        print("* vina-gpu")
        print("* mgltools")
        print("* gromacs")
        print("* acpype")

def create_input(simulation_type):
    if simulation_type == 'lavina':
        create_lavina()
        print("LADOCK_lavina, subdirectories, 'config_lavina.txt', and 'ligand_link.txt' have been successfully created in the 'LADOCK_lavina' directory. Please edit them according to your needs.")
    
    if simulation_type == 'lavinagpu':   
        create_lavinagpu()
        print("LADOCK_lavinagpu, subdirectories, 'config_lavina.txt', and 'ligand_link.txt' have been successfully created in the 'LADOCK_lavinagpu' directory. Please edit them according to your needs.")

    elif simulation_type == 'la2vina':
        create_la2vina()
        print("LADOCK_la2vina, subdirectories, and 'config_la2vina.txt' have been successfully created in the 'LADOCK_la2vina' directory. Please edit them according to your needs.")

    elif simulation_type == 'ladock4':
        create_ladock4()
        print("LADOCK_ladock4, subdirectories, 'config_ladock4.txt', and 'ligand_link.txt' have been successfully created in the 'LADOCK_ladock4' directory. Please edit them according to your needs.")

    elif simulation_type == 'ladockgpu':
        create_ladockgpu()
        print("LADOCK_ladockgpu, subdirectories, 'config_ladockgpu.txt', and 'ligand_link.txt' have been successfully created in the 'LADOCK_ladock4' directory. Please edit them according to your needs.")

    elif simulation_type == 'gmxprolig':
        create_lagmx()
        print("LADOCK_gmxprolig and subdirectory created successfully.")
       
    elif simulation_type == 'gmxliglig':
        create_ladyna3lig()
        print("LADOCK_gmxliglig and subdirectory created successfully.")

def run_simulation(simulation_type):
    
    if simulation_type == 'lavina':
        import config.lavinaConfig
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, 'lavina.py')
        subprocess.run(['python3', path])
        
    if simulation_type == 'lavinagpu':
        import config.lavinagpuConfig
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, 'lavinagpu.py')
        subprocess.run(['python3', path])

    elif simulation_type == 'la2vina':
        import config.la2vinaConfig
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, 'la2vina.py')
        subprocess.run(['python3', path])

    elif simulation_type == 'ladock4':
        import config.ladock4Config
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, 'ladock4.py')
        subprocess.run(['python3', path])

    elif simulation_type == 'ladockgpu':
        import config.ladockgpuConfig
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, 'ladockgpu.py')
        subprocess.run(['python3', path])

    elif simulation_type == 'gmxprolig':
        import config.gmxproligConfig
        current_dir = os.path.dirname(os.path.abspath(__file__))   
        path = os.path.join(current_dir, 'gmxprolig.py')
        subprocess.run(['python3', path])        
       
    elif simulation_type == 'gmxliglig':
        import config.gmxligligConfig
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, 'gmxliglig.py')
        subprocess.run(['python3', path])
        
    elif simulation_type == 'test':
        import config.gmxligligConfig
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, 'test.py')
        subprocess.run(['python3', path])
        
if __name__ == "__main__":
    main()
