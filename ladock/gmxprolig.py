#!/usr/bin/env python

import os
import subprocess
import glob
import shutil
from datetime import datetime
from main import *

def run_pdb2gmx(receptor_files):
    for receptor_file in receptor_files:
        rec_name = os.path.splitext(receptor_file)[0]        
        command = f'gmx pdb2gmx -f {receptor_file} -o {rec_name}_NEW.pdb -ignh'
        
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"Perintah GROMACS berhasil dijalankan untuk {receptor_file}.")
        except subprocess.CalledProcessError as e:
            print(f"Perintah GROMACS gagal dengan kode keluar {e.returncode} untuk {receptor_file}.")
        except Exception as e:
            print(f"Terjadi kesalahan: {str(e)}")
   
def process_mol2_files_in_directories(mol2_files):
    for mol2_file in mol2_files:
        lig_name = os.path.splitext(os.path.basename(mol2_file))[0]
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

def create_complex_pdb(complex_pdb, receptor_new_files, mol2_files):
    complex_data = ""
    
    try:
        for receptor_new_file in receptor_new_files:
            with open(receptor_new_file, 'r') as rec_new:
                rec_data = rec_new.read()
                complex_data += rec_data

        for mol2_file in mol2_files:
            lig_name = os.path.splitext(os.path.basename(mol2_file))[0]
            ligand_pdb = f"{lig_name}.acpype/{lig_name}_NEW.pdb"
            with open(ligand_pdb, 'r') as ligand_file:
                ligand_data = ligand_file.read()
                complex_data += ligand_data

        with open(complex_pdb, 'w') as complex_file:
            complex_file.write(complex_data)
        print(f"File {complex_pdb} berhasil dibuat.")
    except Exception as e:
        print(f"Terjadi kesalahan: {str(e)}")

def remove_lines_from_pdb_file(input_pdb_file, output_pdb_file):
    try:
        with open(input_pdb_file, 'r') as input_file:
            lines = input_file.readlines()

        # Filter out lines containing "MODEL", "ENDMDL", and "REMARK"
        filtered_lines = [line for line in lines if not any(keyword in line for keyword in ["MODEL", "ENDMDL", "TER", "REMARK"])]

        with open(output_pdb_file, 'w') as output_file:
            output_file.writelines(filtered_lines)
        
        print(f"File {output_pdb_file} berhasil dibuat tanpa baris yang mengandung 'MODEL', 'ENDMDL', dan 'REMARK'.")
    except Exception as e:
        print(f"Terjadi kesalahan: {str(e)}")

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
    
    with open('ligand.itp', 'w') as complex_itp_file:
        complex_itp_file.write(complex_itp)
    print("File ligand.itp berhasil dibuat.")

def create_topol_top(mol2_files, receptor_new_files):
    # Baris yang ingin ditambahkan
    new_lines1 = "; Include ligand parameters\n#include \"ligand.itp\"\n\n"
    new_lines2 = ""  # Inisialisasi string new_lines2
    new_lines3 = ""  # Inisialisasi string new_lines3

    for file in mol2_files:
        lig_id = extract_ligand_id_from_mol2(file)
        new_lines2 += f'#include "{lig_id}-posre.itp"\n'
        new_lines3 += f'{lig_id}\t\t\t1\n'

    with open("topol.top", 'r') as topol_file:
        lines = topol_file.readlines()

    index1 = None
    index2 = None

    for i, line in enumerate(lines):
        if line.strip() == "; Include forcefield parameters":
            index1 = i
        if line.strip() == "; Include Position restraint file":
            index2 = i

    if index1 is not None:
        # Menambahkan new_lines1 setelah baris yang sesuai
        lines.insert(index1 + 3, new_lines1)

    if index2 is not None:
        # Menambahkan new_lines2 setelah baris yang sesuai
        lines.insert(index2 + 4, new_lines2)

    # Menambahkan new_lines3 pada bagian terakhir
    lines.extend([new_lines3])

    # Menulis ulang isi file topol.top dengan baris tambahan
    with open("topol.top", 'w') as topol_file:
        topol_file.writelines(lines)
   
def generate_box_pdb(box_pdb, box_type, distance):
    try:
        subprocess.run(["gmx", "editconf", "-f", "complex_clean.pdb", "-o", box_pdb, "-bt", box_type, "-d", str(distance), "-c"])
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

def process_ligand_restraint(mol2_files):
    for mol2_file in mol2_files:
        lig_name = os.path.splitext(os.path.basename(mol2_file))[0]
        lig_id = extract_ligand_id_from_mol2(mol2_file)
        
        # Mencari semua file ligand_new_pdb dalam direktori yang sesuai
        ligand_new_pdbs = glob.glob(f"{lig_name}.acpype/{lig_name}_NEW.pdb")

        for ligand_new_pdb in ligand_new_pdbs:
            try:
                # Membuat indeks untuk ligand
                make_ndx_command = [
                    'gmx', 'make_ndx',
                    '-f', ligand_new_pdb,
                    '-o', f'{lig_id}-index.ndx'
                ]
                subprocess.run(make_ndx_command, check=True)

                # Menjalankan perintah gmx genrestr
                genrestr_command = [
                    'gmx', 'genrestr',
                    '-f', ligand_new_pdb,
                    '-n', f'{lig_id}-index.ndx',
                    '-o', f'{lig_id}-posre.itp',
                    '-fc', '1000', '1000', '1000'
                ]
                subprocess.run(genrestr_command, check=True)

                print(f"File {lig_id}-posre.itp berhasil dibuat untuk {ligand_new_pdb}")
            except subprocess.CalledProcessError as e:
                print(f"Terjadi kesalahan saat menjalankan perintah: {str(e)}")
            
def replace_text_in_files(mol2_files):
    mdp_files = glob.glob(os.path.join(".", "*.mdp"))    
    combined_lig_id = "_".join(extract_ligand_id_from_mol2(file) for file in mol2_files)    
    for mdp_file in mdp_files:
        with open(mdp_file, 'r') as file:
            file_data = file.read()            
            file_data = file_data.replace("Protein_UNK", f"Protein_{combined_lig_id}")

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
    with open('config_gmxprolig.txt', 'r') as config_file:
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
    complex_pdb = "complex.pdb"  
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
            receptor_files = glob.glob('*rec.pdb')            
            copy_mdp_files(current_directory, subdir)       
            replace_text_in_files(mol2_files)
            process_mol2_files_in_directories(mol2_files) 
            run_pdb2gmx(receptor_files)
            receptor_new_files = glob.glob('*_NEW.pdb')                     
            create_complex_pdb(complex_pdb, receptor_new_files, mol2_files)  
            remove_lines_from_pdb_file(complex_pdb, "complex_clean.pdb")            
            create_complex_itp(mol2_files)
            create_topol_top(mol2_files, receptor_new_files)
            generate_box_pdb(box_pdb, box_type, distance)
            solvate_system(box_pdb, spc216_gro, topol_top, solv_gro)
            ionization(solv_gro, topol_top, ions_gro)
            minimization(ions_gro, topol_top)
            process_ligand_restraint(mol2_files)
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
    current_directory =  (os.path.join(os.getcwd(), "LADOCK_gmxprolig"))      
    
    if os.path.exists(current_directory):
        os.chdir(current_directory)
        main() 
    else:
        print("Your job directory (LADOCK_gmxprolig) is not ready. Please create it using:")
        print("ladock --create gmxprolig")
