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

def run_command(command):
    subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def calculate_gridbox_center(structure):
    print(f"  Calculating GridBox Center")    
    model = structure[0]  # Get the first model from the PDB structure
    atoms = model.get_atoms()  # Get the list of atoms in the model
    
    x_sum = y_sum = z_sum = 0.0
    num_atoms = 0
    
    for atom in atoms:
        x, y, z = atom.get_coord()
        x_sum += x
        y_sum += y
        z_sum += z
        num_atoms += 1
    
    x_center = round(x_sum / num_atoms, 3)
    y_center = round(y_sum / num_atoms, 3)
    z_center = round(z_sum / num_atoms, 3)
    
    return x_center, y_center, z_center

def read_energy_and_save_to_csv(counter, output_model_dir, lig_name, csv_path):    
    vina_output_file = os.path.join(output_model_dir, f"output_{lig_name}.pdbqt")
        
    with open(vina_output_file, 'r') as vina_output:
        lines = vina_output.readlines()

    energy = None
    for line in lines:
        if line.startswith("REMARK VINA RESULT:"):
            energy = float(line.split()[3])
            break  # Stop reading after the first matching line

    if energy is not None:
        with open(csv_path, 'a') as csv_file:
            csv_file.write(f"{counter['value']}, {lig_name.replace('output_', '').replace('_minimized', '').upper()}, {energy:.3f}\n")
    else:
        print(f"No REMARK VINA RESULT line found in the output of {os.path.basename(lig_pdbqt_path)}")

def process_docking_ligand_native(counter, size_x, size_y, size_z, thread, search_depth, max_workers, num_iterations, csv_path, prepare_receptor, receptor_name, prepare_ligand, lig_name, directory, ensamble, output_model_dir, cpu):

    # Get ligand's center coordinates
    parser = PDBParser()
    lig_pdb_path = os.path.join(current_directory, directory, f"{lig_name}.pdb")

    structure = parser.get_structure("ligand", lig_pdb_path)
    x_center, y_center, z_center = calculate_gridbox_center(structure)

    # Create config.txt with updated values
    with open('config.txt', 'w') as config_file:
        config_file.write(f"size_x = {size_x}\n")
        config_file.write(f"size_y = {size_y}\n")
        config_file.write(f"size_z = {size_z}\n")
        config_file.write(f"center_x = {x_center}\n")
        config_file.write(f"center_y = {y_center}\n")
        config_file.write(f"center_z = {z_center}\n")
        config_file.write(f"thread = {thread}\n")
        config_file.write(f"search_depth = {search_depth}\n")        
        config_file.write("# Script written by:\n")
        config_file.write("# La Ode Aman\n")
        config_file.write("# laodeaman.ai@gmail.com\n")
        config_file.write("# laode_aman@ung.ac.id\n")
        config_file.write("# Universitas Negeri Gorontalo, Indonesia\n")

    print("  Docking parameter:")
    print("\tsize_x =", size_x)
    print("\tsize_y =", size_y)
    print("\tsize_z =", size_z)
    print("\tcenter_x =", x_center)
    print("\tcenter_y =", y_center)
    print("\tcenter_z =", z_center)
    print("\tthread =", thread)
    print("\tsearch_depth =", search_depth)
    print("\tmax_workers =", max_workers)

    # Create energy_summary.csv
     
    if not os.path.exists(csv_path):
        with open(csv_path, 'w') as csv_file:
            csv_file.write("No., Ligand ID (name), Binding Affinity (kcal per mol)\n")
    else:
        pass

      
    # Prepare receptor and native ligand
    run_command(f'{prepare_receptor} -r {receptor_name}.pdb')
    run_command(f'{prepare_ligand} -l {lig_name}.pdb')
    
    # Run Vina
    print(f"  Docking and save energy")
    run_command(f'{ensamble} --receptor {receptor_name}.pdbqt --ligand {lig_name}.pdbqt --config config.txt --out output_{lig_name}.pdbqt')


    # Move files
    os.rename(f"output_{lig_name}.pdbqt", os.path.join(output_model_dir, f"output_{lig_name}.pdbqt"))

    
    # Save to csv   
    read_energy_and_save_to_csv(counter, output_model_dir, lig_name, csv_path)
         
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

def process_ligand_link(ligand_url, ligand_dir, ligand_tmp_dir, lig_name, num_iterations, max_workers, size_x, size_y, size_z, thread, search_depth, csv_path, prepare_receptor, receptor_name, prepare_ligand, directory, ensamble, output_model_dir, counter, cpu):
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
        process_docking_ligand_test(ligand_tmp_dir, counter, size_x, size_y, size_z, thread, search_depth, num_iterations, csv_path, prepare_receptor, receptor_name, prepare_ligand, lig_name, directory, ensamble, output_model_dir, cpu)

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

def process_docking_ligand(pdbqt_file, counter, prepare_ligand, receptor_name, ensamble, output_model_dir, csv_path):
    try:
        # preparing, docking and save energy        
        filename = os.path.splitext(os.path.basename(pdbqt_file))[0]      
        output_file = os.path.join(output_model_dir, f'output_{filename}.pdbqt')
        run_command(f'{ensamble} --receptor {receptor_name}.pdbqt --ligand {pdbqt_file} --config config.txt --out {output_file}')
        read_energy_and_save_to_csv(counter, output_model_dir, filename, csv_path)    
        os.remove(pdbqt_file)
    except Exception as e:
        print(f"An error occurred while processing a ligand: {e}")

def process_docking_ligand_test(ligand_tmp_dir, counter, size_x, size_y, size_z, thread, search_depth, num_iterations, csv_path, prepare_receptor, receptor_name, prepare_ligand, lig_name, directory, ensamble, output_model_dir, max_workers):
    os.chdir(directory)
    ligand_files = glob.glob(os.path.join(ligand_tmp_dir, '*.pdbqt'))
    num_ligand_files = len(ligand_files)
    
    counter["value"] += 1
    if num_ligand_files == 0:        
        return

    for pdbqt_file in tqdm(ligand_files, desc="  Docking and save energy"):
        process_docking_ligand(pdbqt_file, counter, prepare_ligand, receptor_name, ensamble, output_model_dir, csv_path)

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
    config_file = 'config_lavinagpu.txt'
    
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

    ensamble = config_variables['vina_gpu']    
    thread = config_variables['thread']
    search_depth = config_variables['search_depth']    
    num_iterations = 1000
    cpu = os.cpu_count()
    max_workers = os.cpu_count()    
       
    size_x = float(config_variables['size_x'])  
    size_y = float(config_variables['size_y'])  
    size_z = float(config_variables['size_z'])  
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
                process_docking_ligand_native(counter, size_x, size_y, size_z, thread, search_depth, max_workers, num_iterations, csv_path, prepare_receptor, receptor_name, prepare_ligand, lig_name, directory, ensamble, output_model_dir, cpu)
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
                                    process_ligand_link(ligand_url, ligand_dir, ligand_tmp_dir, lig_name, num_iterations, max_workers, size_x, size_y, size_z, thread, search_depth, csv_path, prepare_receptor, receptor_name, prepare_ligand, directory, ensamble, output_model_dir, counter, cpu)
                    
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
                        process_docking_ligand_test(ligand_tmp_dir, counter, size_x, size_y, size_z, thread, search_depth, num_iterations, csv_path, prepare_receptor, receptor_name, prepare_ligand, lig_name, directory, ensamble, output_model_dir, cpu)
                            
                sort_and_rewrite_csv(csv_path)               
                print(f"Success: {dir_name.upper()}")

            else:
                print(f"Skipping docking: receptor or reference ligand in {dir_name} not found.")
                os.chdir(current_directory)

    # Print developer's note, contact, and citation listdir
    print_dev(developer_note, developer_contact, citation_list)

if __name__ == "__main__":
    current_directory =  (os.path.join(os.getcwd(), "LADOCK_lavinagpu")) 
    if os.path.exists(current_directory):  
        main() 
    else:
        print("Your job directory (LADOCK_lavinagpu) is not ready. Please create it using:")
        print("ladock --create lavinagpu")
