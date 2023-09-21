#!/usr/bin/env python

import os
import subprocess
import glob
import shutil
from datetime import datetime
from main import *
   
def process_mol2_files_in_directories(mol2_files) :                
    # Prepare the PDB file using acpype for gas phase    
    for mol2_file in mol2_files:
        lig_name = os.path.splitext(mol2_file)[0]
        subprocess.run(['obabel', mol2_file, '-opdb', '-h', '-O', f'{lig_name}.pdb'])
        subprocess.run(['acpype', '-i', f'{lig_name}.pdb', '-c', 'gas'])

def extract_ligand_id_from_mol2(mol2_file):
    lig_id = None
    with open(mol2_file, 'r') as file:
        for line in file:
            if line.startswith('@<TRIPOS>ATOM'):
                # Mencari baris yang berisi ID ligand
                lig_id_line = next(file)
                parts = lig_id_line.split()
                if len(parts) > 1:
                    lig_id = parts[7][:3]
                break
    print(f"Ligand ID extracted: {lig_id}")
    return lig_id

def create_complex_pdb(output_file, mol2_files):  
    complex_data = ""     
    for mol2_file in mol2_files:
        lig_name = os.path.splitext(mol2_file)[0]
        ligand_pdb = f"{lig_name}.acpype/{lig_name}_NEW.pdb"    
        with open(ligand_pdb, 'r') as ligand_file:
            ligand_data = ligand_file.read()
            complex_data += ligand_data
            
    with open(output_file, 'w') as complex_file:
        complex_file.write(complex_data)
        print(f"File {output_file} berhasil dibuat.")   

def create_complex_itp(mol2_files):
    atomtype_itp = ""    
    moleculetype_itp = ""
    complex_itp = ""
    
    for mol2_file in mol2_files:
        extract_ligand_id_from_mol2(mol2_file)            
        lig_name = os.path.splitext(mol2_file)[0]
        ligand_itp = f"{lig_name}.acpype/{lig_name}_GMX.itp"         
        with open(ligand_itp, 'r') as f:        
            flag_atomtypes = False            
            flag_moleculetype = False            
            
            for line in f:                
                if '[ atomtypes ]' in line:
                    flag_atomtypes = True
                if not line.strip():
                    flag_atomtypes = False                
                if flag_atomtypes:
                        if line not in complex_itp:
                            atomtype_itp += line

                if '[ moleculetype ]' in line:                    
                    flag_moleculetype = True                
                if flag_moleculetype:
                    moleculetype_itp += line               
            complex_itp = atomtype_itp + "\n"                                 
            complex_itp += moleculetype_itp + "\n"    
    
    # Ganti nama ligand dengan ID sesuai indeks
    for mol2_file in mol2_files:
        lig_name = os.path.splitext(mol2_file)[0]        
        lig_id = extract_ligand_id_from_mol2(mol2_file)
        print(f"{lig_id} : {lig_name}")
        complex_itp = complex_itp.replace(lig_name, lig_id)
    
    with open('complex.itp', 'w') as complex_itp_file:
        complex_itp_file.write(complex_itp)
    print("File complex.itp berhasil dibuat.")

def create_topol_top(mol2_files):
    current_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    combined_lig_id = "_".join(extract_ligand_id_from_mol2(file) for file in mol2_files)

    with open("topol.top", 'w') as topol_file:
        topol_file.write(f"; topol.top created by ladock (v: 2023) on {current_date}\n")
        topol_file.write("\n")
        topol_file.write("; Include forcefield parameters\n")
        topol_file.write("#include \"amber99sb-ildn.ff/forcefield.itp\"\n")
        topol_file.write("\n")
        topol_file.write("; Include complex topology\n")
        topol_file.write("#include \"complex.itp\"\n")
        topol_file.write("\n")
        topol_file.write("; Include water topology\n")
        topol_file.write("#include \"amber99sb-ildn.ff/tip3p.itp\"\n")
        topol_file.write("\n")
        topol_file.write("#ifdef POSRES_WATER\n")
        topol_file.write("; Position restraint for each water oxygen\n")
        topol_file.write("[ position_restraints ]\n")
        topol_file.write("   1    1       1000       1000       1000\n")
        topol_file.write("#endif\n")
        topol_file.write("\n")
        topol_file.write("; Include topology for ions\n")
        topol_file.write("#include \"amber99sb-ildn.ff/ions.itp\"\n")
        topol_file.write("\n")
        topol_file.write("[ system ]\n")
        topol_file.write(f"{combined_lig_id}\n")
        topol_file.write("\n")
        topol_file.write("[ molecules ]\n")
        for file in mol2_files:
            lig_id = extract_ligand_id_from_mol2(file)
            topol_file.write(f"{lig_id}                  1\n")

    print("File topol.top berhasil dibuat.")
    
def generate_box_pdb(input_pdb, box_pdb, box_type, distance):
    try:
        subprocess.run(["gmx", "editconf", "-f", input_pdb, "-o", box_pdb, "-bt", box_type, "-d", str(distance), "-c"])
        print(f"File {box_pdb} berhasil dibuat.")
    except subprocess.CalledProcessError:
        print("Gagal menjalankan perintah gmx editconf.") 

def solvate_system(box_pdb, spc216_gro, topol_top, solv_gro):
    try:
        subprocess.run(["gmx", "solvate", "-cp", box_pdb, "-cs", spc216_gro, "-p", topol_top, "-o", solv_gro])
        print(f"File {solv_gro} berhasil dibuat.")
    except subprocess.CalledProcessError:
        print("Gagal menjalankan perintah gmx solvate.")

def ionization(solv_gro, topol_top, ions_gro):
    try:
        # Generate TPR file
        subprocess.run(["gmx", "grompp", "-f", "ions.mdp", "-c", solv_gro, "-p", topol_top, "-o", "ions.tpr", "-maxwarn", "1"])

        # Run genion for ionization
        genion_cmd = f"gmx genion -s ions.tpr -o {ions_gro} -p {topol_top} -pname NA -nname CL -neutral"
        subprocess.run(["bash", "-c", genion_cmd])
        
        print(f"File {ions_gro} berhasil dibuat.")
    except subprocess.CalledProcessError:
        print("Gagal menjalankan perintah gmx genion.")
 
def minimization(ions_gro, topol_top):
    try:
        # Generate TPR file for minimization
        subprocess.run(["gmx", "grompp", "-f", "em.mdp", "-c", ions_gro, "-p", topol_top, "-o", "em.tpr", "-maxwarn", "1"])
        
        # Run energy minimization using mdrun
        subprocess.run(["gmx", "mdrun", "-v", "-deffnm", "em"])
        
        # Calculate potential energy
        subprocess.run(["echo", "-e", "10\n0", "|", "gmx", "energy", "-f", "em.edr", "-o", "potential.xvg"])
        
        print("Minimization completed.")
    except subprocess.CalledProcessError:
        print("Gagal menjalankan perintah GROMACS.")
 
def make_index_system(em_gro, index_ndx):
    try:
        # Create index file for the system
        subprocess.run(["gmx", "make_ndx", "-f", em_gro, "-o", index_ndx])
        
        print("Index file created.")
    except subprocess.CalledProcessError:
        print("Gagal menjalankan perintah GROMACS.")

def copy_mdp_files(current_directory, subdir):
    mdp_files = glob.glob(os.path.join(current_directory, "*.mdp"))    
    try:
        for mdp_file in mdp_files:
            mdp_filename = os.path.basename(mdp_file)
            destination_path = os.path.join(".", mdp_filename)
            shutil.copyfile(mdp_file, mdp_filename)
        
        print("MDP files copied successfully.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        # Handle the exception as needed
    finally:
        # Optionally, perform cleanup or additional actions here
        pass
            
def replace_text_in_files(mol2_files):
    mdp_files = glob.glob(os.path.join(".", "*.mdp"))    
    combined_lig_id = "_".join(extract_ligand_id_from_mol2(file) for file in mol2_files)    
    for mdp_file in mdp_files:
        with open(mdp_file, 'r') as file:
            file_data = file.read()
            file_data = file_data.replace("Protein_UNK", combined_lig_id)

        with open(mdp_file, 'w') as file:
            file.write(file_data)

def run_nvt_simulation():
    # Jalankan grompp untuk NVT
    grompp_cmd = [
        'gmx', 'grompp',
        '-f', 'nvt.mdp',
        '-c', 'em.gro',
        '-r', 'em.gro',
        '-p', 'topol.top',
        '-n', 'index.ndx',
        '-o', 'nvt.tpr',
        '-maxwarn', '5'
    ]
    subprocess.run(grompp_cmd, check=True)

    # Jalankan mdrun untuk NVT
    mdrun_cmd = [
        'gmx', 'mdrun',
        '-v',
        '-s', 'nvt.tpr',
        '-deffnm', 'nvt'
    ]
    subprocess.run(mdrun_cmd, check=True)

def run_npt_simulation():
    # Jalankan grompp untuk NPT
    grompp_cmd = [
        'gmx', 'grompp',
        '-f', 'npt.mdp',
        '-c', 'nvt.gro',
        '-t', 'nvt.cpt',
        '-r', 'nvt.gro',
        '-p', 'topol.top',
        '-n', 'index.ndx',
        '-o', 'npt.tpr',
        '-maxwarn', '5'
    ]
    subprocess.run(grompp_cmd, check=True)

    # Jalankan mdrun untuk NPT
    mdrun_cmd = [
        'gmx', 'mdrun',
        '-v',
        '-s', 'npt.tpr',
        '-deffnm', 'npt'
    ]
    subprocess.run(mdrun_cmd, check=True)

def run_production_simulation():
    # Jalankan grompp untuk produksi
    grompp_cmd = [
        'gmx', 'grompp',
        '-f', 'md.mdp',
        '-c', 'npt.gro',
        '-t', 'npt.cpt',
        '-p', 'topol.top',
        '-n', 'index.ndx',
        '-o', 'md.tpr',
        '-maxwarn', '5'
    ]
    subprocess.run(grompp_cmd, check=True)

    # Jalankan mdrun untuk produksi
    mdrun_cmd = [
        'gmx', 'mdrun',
        '-v',
        '-s', 'md.tpr',
        '-deffnm', 'md'
    ]
    subprocess.run(mdrun_cmd, check=True)

def print_dev(developer_note, developer_contact, citation_list):
    print("")
    print(developer_note)
    print("")
    print(developer_contact)
    print("")
    print(citation_list)
    print("")

def main():
    print_dev(developer_note, developer_contact, citation_list)
    # Membaca isi file 'config_tmp.txt'
    with open('config_gmxliglig.txt', 'r') as config_file:
        config_lines = config_file.readlines()
    
    # Membuat kamus (dictionary) untuk menyimpan variabel-variabel
    config_variables = {}

    # Memproses setiap baris dalam file 'config.txt'
    for line in config_lines:
            # Mengabaikan baris yang dimulai dengan '#' (komentar)
            if not line.strip().startswith('#'):
                # Memisahkan nama variabel dan nilai variabel
                try:
                    name, value = line.strip().split('=')
                    # Membersihkan spasi dan karakter lainnya
                    name = name.strip()
                    value = value.strip()
                    # Menambahkan variabel ke dalam kamus
                    config_variables[name] = value
                except ValueError:
                    pass

    box_type = config_variables['box_type']
    distance = config_variables['distance']
    spc216_gro = config_variables['coordinate_file']
    current_directory = os.getcwd()  
    input_pdb = "complex.pdb"  
    box_pdb = "box.pdb"    
    topol_top = "topol.top"    
    solv_gro = "solv.gro" 
    ions_gro = "ions.gro"
    em_gro = "em.gro"
    index_ndx = "index.ndx"
    
    for subdir in os.listdir(current_directory):
        if os.path.isdir(subdir):
            os.chdir(subdir)         
            mol2_files = glob.glob('*.mol2')
            copy_mdp_files(current_directory, subdir)       
            replace_text_in_files(mol2_files)
            process_mol2_files_in_directories(mol2_files)             
            create_complex_pdb("complex.pdb", mol2_files)     
            create_complex_itp(mol2_files)
            create_topol_top(mol2_files)
            generate_box_pdb(input_pdb, box_pdb, box_type, distance)
            solvate_system(box_pdb, spc216_gro, topol_top, solv_gro)
            ionization(solv_gro, topol_top, ions_gro)
            minimization(ions_gro, topol_top)
            make_index_system(em_gro, index_ndx)            
            os.chdir('..')
            
    for subdir in os.listdir(current_directory):
        if os.path.isdir(subdir):
            os.chdir(subdir)         
            run_nvt_simulation()
            run_npt_simulation()
            run_production_simulation()            
            os.chdir('..')
            
    print_dev(developer_note, developer_contact, citation_list)
    
if __name__ == "__main__":
    current_directory =  (os.path.join(os.getcwd(), "LADOCK_gmxliglig"))      
    
    if os.path.exists(current_directory):
        os.chdir(current_directory)
        main() 
    else:
        print("Your job directory (LADOCK_gmxliglig) is not ready. Please create it using:")
        print("ladock --create gmxliglig")
