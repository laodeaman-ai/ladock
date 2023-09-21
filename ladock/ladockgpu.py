import os
import subprocess
import glob
import shutil
import urllib.request
from Bio.PDB import PDBParser
from tqdm import tqdm
import argparse
import gzip
import csv
from main import *
import concurrent.futures
from concurrent.futures import as_completed
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

def run_command(command):
    subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def read_energy_and_save_to_csv(counter, output_file, csv_path, lig_name):
    # Buka file untuk membaca
    with open(output_file, 'r') as file:
        content = file.read()

    # Cari baris yang mengandung "RMSD TABLE"
    in_rmsd_table = False
    for line in content.split('\n'):
        if "RMSD TABLE" in line:
            in_rmsd_table = True
            continue
        if in_rmsd_table:
            if line.strip():
                words = line.split()
                try:
                    energy = float(words[3])
                    break
                except (ValueError, IndexError):
                    pass
    else:
        print("Tabel RMSD tidak ditemukan atau tidak ada nilai energy pada baris pertama kolom keempat.")

    # Jika energy ditemukan, simpan ke dalam CSV
    if energy is not None:
        with open(csv_path, 'a') as csv_file:
            csv_file.write(f"{counter['value']}, {lig_name.replace('_minimized', '').upper()}, {energy:.3f}\n")

def process_docking_ligand_native(prepare_gpf, counter, max_workers, num_iterations, csv_path, prepare_receptor, receptor_name, prepare_ligand, lig_name, directory, ensamble, output_model_dir, nrun):

    # Create energy_summary.csv
    print("  Docking and save energy ...")
    if not os.path.exists(csv_path):
        with open(csv_path, 'w') as csv_file:
            csv_file.write("No., Ligand ID (name), Binding Affinity (kcal per mol)\n")
    else:
        pass

    counter ["value"] += 1
    
    # Prepare receptor & native, gpf, autogrid4, dpf, autodock4
    output_file = os.path.join(output_model_dir, f"{lig_name}.dlg")
    
    run_command(f'{prepare_ligand} -l {lig_name}.pdb')    
    run_command(f'{prepare_receptor} -r {receptor_name}.pdb -o {lig_name}_{receptor_name}.pdbqt -A checkhydrogens') 
    run_command(f'{prepare_gpf} -l {lig_name}.pdbqt -r {lig_name}_{receptor_name}.pdbqt -y -o {lig_name}_{receptor_name}.gpf')
    run_command(f'autogrid4 -p {lig_name}_{receptor_name}.gpf -l {lig_name}_{receptor_name}.glg')
    run_command(f'{ensamble} --ffile {lig_name}_{receptor_name}.maps.fld --lfile {lig_name}.pdbqt --nrun {nrun}')
    os.rename(f"{lig_name}.dlg", output_file)
    
    # Save to csv   
    read_energy_and_save_to_csv(counter, output_file, csv_path, lig_name)
         
def process_in_ligand_dir(ligand_file, ligand_dir, ligand_tmp_dir):    
    os.chdir(ligand_dir)
    ligand_name = os.path.basename(ligand_file)
    print(f"\n  Ligand: {ligand_name}")
    print(f"  Splitting {ligand_name} to single molecule")

    ligand_test = os.path.splitext(os.path.basename(ligand_file))[0]
    try:
        if ligand_file.endswith((".smiles", ".smi")):
            preparing_smi_file(ligand_file)            
            command = f"obabel {ligand_file} -osmi -O {os.path.join(ligand_tmp_dir, ligand_test)}.smi -m -h"
        elif ligand_file.endswith((".sdf", ".pdb", ".mol", ".mol2")):
            command = f"obabel {ligand_file} -osdf -O {os.path.join(ligand_tmp_dir, ligand_test)}.sdf -m -h"
        else:
            return  # Skip unsupported file formats

        subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while preparing ligand: {ligand_file}")

def process_and_optimize_smi(prepare_ligand, smi_file, ligand_tmp_dir):
    try:
        with open(smi_file, 'r') as file:
            parts = file.readline().split()[1]
            file_name = parts.split('/')[-1]
            new_name = file_name.split('.')[0].upper()

        # Convert .smi to .mol
        obabel_convert_command_mol = f'obabel {smi_file} -omol -h -O {new_name}.mol --gen2d'
        subprocess.run(obabel_convert_command_mol, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Perform geometry optimization (output: .pdb)
        obminimize_input = f"{new_name}.mol"
        obminimize_output = f"{new_name}_minimized.pdb"
        obminimize_command = f'obminimize -o "pdb" {obminimize_input} > {obminimize_output}'
        subprocess.run(obminimize_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        convert_to_pdbqt(obminimize_output, prepare_ligand)
        os.remove(smi_file)
        os.remove(obminimize_input)
        os.remove(obminimize_output)
    except subprocess.CalledProcessError as e:
        pass
        
def process_and_optimize_sdf(prepare_ligand, sdf_file, ligand_tmp_dir):
    try:
        with open(sdf_file, 'r') as file:
            new_name = file.readline().strip()
        # Perform geometry optimization (output: .pdb)
        obminimize_input = sdf_file
        obminimize_output = f"{new_name}_minimized.pdb"
        obminimize_command = f'obminimize -ff Ghemical -o "pdb" {obminimize_input} > {obminimize_output}'
        subprocess.run(obminimize_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        convert_to_pdbqt(obminimize_output, prepare_ligand)
        os.remove(obminimize_input)
        os.remove(obminimize_output)
    except subprocess.CalledProcessError as e:
        pass
        
def process_in_ligand_tmp_dir(prepare_ligand, ligand_tmp_file, ligand_tmp_dir, max_workers):
    os.chdir(ligand_tmp_dir)
    smi_futures = []
    sdf_futures = []
    smi_files = [f for f in os.listdir(ligand_tmp_dir) if f.endswith(".smi")]
    sdf_files = [f for f in os.listdir(ligand_tmp_dir) if f.endswith(".sdf")]
        
    if smi_files:  # Periksa apakah ada file .smi yang tersedia
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
                for smi_file in smi_files:
                    smi_futures.append(executor.submit(process_and_optimize_smi, prepare_ligand, smi_file, ligand_tmp_dir))
                
                for future in tqdm(concurrent.futures.as_completed(smi_futures), total=len(smi_futures), desc="  Geometry optimization"):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"An error occurred while processing : {e}")
    
    elif sdf_files:  # Periksa apakah ada file .sdf yang tersedia
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                for sdf_file in sdf_files:
                    sdf_futures.append(executor.submit(process_and_optimize_sdf, prepare_ligand, sdf_file, ligand_tmp_dir))
                
                for future in tqdm(concurrent.futures.as_completed(sdf_futures), total=len(sdf_futures), desc="  Geometry optimization"):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"An error occurred while processing {e}")

def process_ligand_link(prepare_dpf, prepare_gpf, ligand_url, ligand_dir, ligand_tmp_dir, lig_name, num_iterations, max_workers, csv_path, prepare_receptor, receptor_name, prepare_ligand, directory, ensamble, output_model_dir, counter, nrun):
    retry_count = 0
    max_retries = 10
    ligand_file = None  # Inisialisasi variabel ligand_file
    
    while retry_count < max_retries:
        try:
            # Download the ligand file
            ligand_file = os.path.join(ligand_dir, os.path.basename(ligand_url))
            urllib.request.urlretrieve(ligand_url, ligand_file)
            
            break
            
        except Exception as e:
            print(e)
            print(f"Retrying to download: {ligand_url}")
            retry_count += 1

    if retry_count == max_retries:
        print(f"Failed to download ligand after {max_retries} retries.")
    
    if ligand_file and ligand_file.endswith(".gz"):
        # Extract the .gz file
        extracted_file = os.path.splitext(ligand_file)[0]  # Remove .gz extension
        with gzip.open(ligand_file, 'rb') as gz_file, open(extracted_file, 'wb') as out_file:
            out_file.writelines(gz_file)
        
        # Remove the downloaded .gz file
        os.remove(ligand_file)
        
        # Update the ligand_file to the extracted file
        ligand_file = extracted_file
    
    # Read all lines into a list
    with open(ligand_file, 'r') as file:
        lines = file.readlines()
    
    # Write back the content to the file, excluding the first line
    with open(ligand_file, 'w') as file:
        file.writelines(lines[0:])
    
    if ligand_file.endswith(".smi"):
        preparing_smi_file(ligand_file)
    process_in_ligand_dir(ligand_file, ligand_dir, ligand_tmp_dir)            
    
    # Remove the downloaded ligand file
    os.remove(ligand_file)
    
    valid_extensions = (".sdf", ".smi")
    ligand_tmp_files = [f for f in os.listdir(ligand_tmp_dir) if any(f.endswith(ext) for ext in valid_extensions)]
    num_ligand_files = len(ligand_tmp_files)

    print(f"  Total molecules = {num_ligand_files}")                    
    print(f"  Using {lig_name}.pdb docking parameter as reference")
    
    for ligand_tmp_file in ligand_tmp_files:
        process_in_ligand_tmp_dir(prepare_ligand, ligand_tmp_file, ligand_tmp_dir, max_workers)        
        process_docking_ligand_test(counter, prepare_receptor, prepare_ligand, prepare_gpf, prepare_dpf, lig_name, output_model_dir, receptor_name, ligand_tmp_dir, csv_path,  ensamble, directory, max_workers, nrun)

def preparing_smi_file(smi_file_path):
    # Membaca file .smi ke dalam list
    with open(smi_file_path, 'r') as file:
        lines = file.readlines()

    # Membuat list untuk menyimpan baris-baris yang valid
    valid_lines = []

    for line in lines:
        # Memeriksa validitas SMILES dengan RDKit
        mol = Chem.MolFromSmiles(line.strip())
        if mol is not None:
            valid_lines.append(line)

    # Menyimpan kembali file .smi yang hanya berisi baris-baris yang valid
    with open(smi_file_path, 'w') as file:
        file.writelines(valid_lines)

def convert_to_pdbqt(pdb_file, prepare_ligand):
    filename = os.path.splitext(os.path.basename(pdb_file))[0]
    pdbqt_file = f"{filename}.pdbqt"
    run_command(f'{prepare_ligand} -l {pdb_file}')
    return pdbqt_file

def process_docking_ligand(directory, pdbqt_file, counter, prepare_receptor, receptor_name, prepare_gpf, lig_name, prepare_dpf,  ensamble, output_model_dir, csv_path, nrun):
    try:
        filename = os.path.splitext(os.path.basename(pdbqt_file))[0]
        output_file = os.path.join(output_model_dir, f'output_{filename}.dlg')
        dest_path = os.path.join(".", filename + ".pdbqt")     
        
        shutil.move(pdbqt_file, dest_path)
        
        # Prepare receptor & native, gpf, autogrid4, dpf, autodock4
        run_command(f'{prepare_receptor} -r {receptor_name}.pdb -o {filename}_{receptor_name}.pdbqt -A checkhydrogens') 
        run_command(f'{prepare_gpf} -l {filename}.pdbqt -r {filename}_{receptor_name}.pdbqt -i {lig_name}_{receptor_name}.gpf -o {filename}_{receptor_name}.gpf')
        run_command(f'autogrid4 -p {filename}_{receptor_name}.gpf -l {filename}_{receptor_name}.glg')        
        run_command(f'{ensamble} --ffile {filename}_{receptor_name}.maps.fld --lfile {filename}.pdbqt --nrun {nrun}')
        os.rename(f"{filename}.dlg", output_file)
        read_energy_and_save_to_csv(counter, output_file, csv_path, filename)
      
        # Remove the original and intermediate files
        file_to_remove = glob.glob(os.path.join(directory, f'{filename}*'))
        for file in file_to_remove:
            os.remove(file)

            
    except Exception as e:
        print(f"An error occurred while processing {filename}: {e}")

def process_docking_ligand_test(counter, prepare_receptor, prepare_ligand, prepare_gpf, prepare_dpf, lig_name, output_model_dir, receptor_name, ligand_tmp_dir, csv_path,  ensamble, directory, max_workers, nrun):
    os.chdir(directory)
    ligand_files = glob.glob(os.path.join(ligand_tmp_dir, '*.pdbqt'))
    num_ligand_files = len(ligand_files)
    
    counter ["value"] += 1
    if num_ligand_files == 0:        
        return

    futures = []
    
    for pdbqt_file in tqdm(ligand_files, desc="  Docking and save energy"):
        process_docking_ligand (directory, pdbqt_file, counter, prepare_receptor, receptor_name, prepare_gpf, lig_name, prepare_dpf,  ensamble, output_model_dir, csv_path, nrun)

def sort_and_rewrite_csv(csv_path):
    try:
        data = []
        with open(csv_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            header = next(csv_reader)  # Baca header (baris pertama)
            for row in csv_reader:
                data.append(row)            

        # Mengurutkan data berdasarkan kolom kedua (Binding Affinity)
        sorted_data = sorted(data, key=lambda x: float(x[2]))  # Menggunakan lambda untuk mengakses kolom kedua (indeks 2) sebagai kunci pengurutan

        # Menulis data yang sudah diurutkan kembali ke file CSV
        with open(csv_path, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_file.write("No., Ligand ID (name), Binding Affinity (kcal per mol)\n")  # Menulis header kembali
            for index, row in enumerate(sorted_data, start=1):
                csv_writer.writerow([index] + row[1:])  # Menulis nomor urutan, nama, dan binding affinity

        print(f"Output CSV file is in the {csv_path}")                
        return True
    except Exception as e:
        print(f"Error while sorting and rewriting CSV: {e}")
        return False

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
    
    os.chdir(current_directory)
    config_file = 'config_ladockgpu.txt'
    
    with open(config_file, 'r') as config_file:
        config_lines = config_file.readlines()
    
    config_variables = {}
    retry_count = 0
    
    for line in config_lines:
        if not line.strip().startswith('#'):
            try:
                name, value = line.strip().split('=')
                name = name.strip()
                value = value.strip()
                config_variables[name] = value
            except ValueError:
                pass

    ensamble = config_variables['autodockgpu']
    nrun = os.cpu_count()
    num_iterations = 1000 
    max_workers = os.cpu_count()         
       
    ligand_dir = os.path.join(current_directory, 'ligand_input')
    ligand_tmp_dir = os.path.join(current_directory, 'ligand_tmp')
    output_dir = os.path.join(current_directory, 'output')
    
    mgl_directory = config_variables['mgl_directory']
    prepare_ligand = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_ligand4.py")
    prepare_receptor = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_receptor4.py")
    prepare_gpf = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_gpf4.py")
    prepare_dpf = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_dpf4.py")
    prepare_lowest_energy = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "write_lowest_energy_ligand.py")
    prepare_summarize_result = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "summarize_results4.py")
    
    print("Removing unnecessary files and directories in job directory...")

    # Create ligand_tmp and output directory if it doesn't exist
    if os.path.exists(ligand_tmp_dir):
        shutil.rmtree(ligand_tmp_dir)
    os.makedirs(ligand_tmp_dir, exist_ok=True)    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True) 
    
    for directory in glob.glob(os.path.join(current_directory, 'model*')):
        if os.path.isdir(directory):
            counter = {"value": 0}
            dir_name = os.path.basename(directory)
            receptor_name = ""
            lig_name = ""
            print(f"\nProcessing {dir_name.upper()}")
            os.chdir(directory)
            receptor_files = glob.glob('*rec.pdb')
            if receptor_files:    
                receptor_name = os.path.splitext(os.path.basename(receptor_files[0]))[0]
            ligand_files = glob.glob('*lig.pdb')
            if ligand_files:    
                lig_name = os.path.splitext(os.path.basename(ligand_files[0]))[0]                
            output_model_dir = os.path.join(output_dir, dir_name)         
            os.makedirs(output_model_dir, exist_ok=True)
            csv_path = os.path.join(output_model_dir, 'energy_summary.csv')
                
            if os.path.exists(f"{receptor_name}.pdb") and os.path.exists(f"{lig_name}.pdb"):
                with open(csv_path, 'w') as csv_file:
                    csv_file.write("No., Ligand ID (name), Binding Affinity (kcal per mol)\n")
                os.chdir(directory)
                print("Docking reference ligand:")
                print(f"  Receptor: {receptor_name}.pdbqt")
                print(f"  Ligand: {lig_name}.pdbqt")
                process_docking_ligand_native(prepare_gpf, counter, max_workers, num_iterations, csv_path, prepare_receptor, receptor_name, prepare_ligand, lig_name, directory, ensamble, output_model_dir, nrun)
                print("\nDocking test ligands: ")
                print(f"  Receptor: {receptor_name}.pdbqt")                  
                os.chdir(ligand_dir)
                
                valid_extensions = (".txt")
                ligand_link_files = [f for f in os.listdir(ligand_dir) if any(f.endswith(ext) for ext in valid_extensions)]
                for ligand_link_path in ligand_link_files:
                    if os.path.exists(ligand_link_path):
                        with open(ligand_link_path, 'r') as f:
                            lines = f.readlines()                        
                        for line in lines:
                            if line.strip():
                                if not line.startswith("#"):
                                    ligand_url = line.strip()
                                    process_ligand_link(prepare_dpf, prepare_gpf, ligand_url, ligand_dir, ligand_tmp_dir, lig_name, num_iterations, max_workers, csv_path, prepare_receptor, receptor_name, prepare_ligand, directory, ensamble, output_model_dir, counter, nrun)
                    
                valid_extensions = (".sdf", ".smi", ".pdb", ".smiles", ".mol", ".mol2")
                ligand_files = [f for f in os.listdir(ligand_dir) if any(f.endswith(ext) for ext in valid_extensions)]
                for ligand_file in ligand_files: 
                    process_in_ligand_dir(ligand_file, ligand_dir, ligand_tmp_dir)
                    
                    valid_extensions = (".sdf", ".smi")
                    ligand_tmp_files = [f for f in os.listdir(ligand_tmp_dir) if any(f.endswith(ext) for ext in valid_extensions)]
                    num_ligand_files = len(ligand_tmp_files)
                    print(f"  Total molecules = {num_ligand_files}")                    
                    print(f"  Using {lig_name}.pdb docking parameter as reference")
                    for ligand_tmp_file in ligand_tmp_files:
                        process_in_ligand_tmp_dir(prepare_ligand, ligand_tmp_file, ligand_tmp_dir, max_workers)                        
                        process_docking_ligand_test(counter, prepare_receptor, prepare_ligand, prepare_gpf, prepare_dpf,  lig_name, output_model_dir, receptor_name, ligand_tmp_dir, csv_path,  ensamble, directory, max_workers, nrun)
                            
                sort_and_rewrite_csv(csv_path)               
                print(f"Success: {dir_name.upper()}")

            else:
                print(f"Skipping docking: receptor or reference ligand in {dir_name} not found.")
                os.chdir(current_directory)

    # Print developer's note, contact, and citation listdir
    print_dev(developer_note, developer_contact, citation_list)

if __name__ == "__main__":
    current_directory =  (os.path.join(os.getcwd(), "LADOCK_ladockgpu")) 
    if os.path.exists(current_directory):  
        main() 
    else:
        print("Your job directory (LADOCK_ladockgpu) is not ready. Please create it using:")
        print("ladock --create ladockgpu")
