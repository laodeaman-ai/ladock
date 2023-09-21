import os
import subprocess
import glob
import shutil
import urllib.request
from Bio.PDB import PDBParser
from tqdm import tqdm
import argparse
from main import vina, vina_split

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

def read_energy_and_save_to_csv(counter, output_model_dir, output_name, csv_path):    
    vina_output_file = os.path.join(output_model_dir, f"{output_name}.pdbqt")
        
    with open(vina_output_file, 'r') as vina_output:
        lines = vina_output.readlines()

    energy = None
    for line in lines:
        if line.startswith("REMARK VINA RESULT:"):
            energy = float(line.split()[3])
            break  # Stop reading after the first matching line

    if energy is not None:
        with open(csv_path, 'a') as csv_file:
            csv_file.write(f"{counter['value']}, {output_name.replace('output_', '').replace('_minimized', '').upper()}, {energy:.3f}\n")
    else:
        print(f"No REMARK VINA RESULT line found in the output of {os.path.basename(lig_pdbqt_path)}")

def process_docking_ligand_native(counter, lig_pdb_path, size_x, size_y, size_z, num_modes, exhaustiveness, cpu, num_iterations, csv_path, prepare_receptor_script, receptor_pdb_path, receptor_pdbqt_path, prepare_ligand_script, lig_pdbqt_path, directory, vina, output_model_dir):
    # Get ligand's center coordinates
    parser = PDBParser()
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
        config_file.write(f"num_modes = {num_modes}\n")
        config_file.write(f"exhaustiveness = {exhaustiveness}\n")
        config_file.write(f"cpu = {cpu}\n")
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
    print("\tnum_modes =", num_modes)
    print("\texhaustiveness =", exhaustiveness)
    print("\tcpu =", cpu)

    # Create energy_summary.csv
    print("  Create energy summary as csv file")
    
    if not os.path.exists(csv_path):
        with open(csv_path, 'w') as csv_file:
            csv_file.write("No., Ligand ID (name), Binding Affinity (kcal per mol)\n")
    else:
        print("CSV file for summary energies is available")
      
    # Prepare receptor and native ligand
    run_command(f'{prepare_receptor_script} -r {receptor_pdb_path} -o {receptor_pdbqt_path}')
    run_command(f'{prepare_ligand_script} -l {lig_pdb_path} -o {lig_pdbqt_path}')
   
    counter ["value"] += 1
    
    # Run Vina
    print(f"  Docking and save energy")
    output_name = f"output_{directory}_lig"
    run_command(f'{vina} --receptor {receptor_pdbqt_path} --ligand {lig_pdbqt_path} --config config.txt --out {output_name}.pdbqt')

    # Move files
    os.rename(f"{output_name}.pdbqt", os.path.join(output_model_dir, f"{output_name}.pdbqt"))

    # Save to csv   
    read_energy_and_save_to_csv(counter, output_model_dir, output_name, csv_path)
            
def process_in_ligand_dir(ligand_dir, ligand_tmp_dir):
    os.chdir(ligand_dir)
    ligand_list = os.listdir(ligand_dir)
    ligand_files = [filename for filename in ligand_list if filename.endswith((".smi",".sdf", ".pdb", ".smiles", ".mol", ".mol2"))]
     
    for filename in ligand_files:
        ligand_test = os.path.splitext(os.path.basename(filename))[0]
        try:
            command = f"obabel {filename} -osmi -d -O {os.path.join(ligand_tmp_dir, ligand_test)}.smi -m"
            subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while preparing ligand: {filename}")

def process_in_main_ligand(main_ligand_dir, main_ligand_tmp):
    os.chdir(main_ligand_dir)
    ligand_list = os.listdir(main_ligand_dir)
    ligand_files = [filename for filename in ligand_list if filename.endswith((".smi",".sdf", ".pdb", ".smiles", ".mol", ".mol2"))]
     
    for filename in ligand_files:
        ligand_test = os.path.splitext(os.path.basename(filename))[0]
        try:
            command = f"obabel {filename} -osmi -d -O {os.path.join(main_ligand_tmp, ligand_test)}.smi -m"
            subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while preparing ligand: {filename}")

def process_in_main_ligand_tmp(main_ligand_tmp, num_iterations):
    os.chdir(main_ligand_tmp)    
    smi_files = [f for f in os.listdir(main_ligand_tmp) if f.endswith(".smi")]
    
    for smi_file in tqdm(smi_files, desc="Generating 3D conformation of main ligands"):
        with open(smi_file, 'r') as file:
            parts = file.readline().split()[1]
            file_name = parts.split('/')[-1]
            new_name = file_name.split('.')[0].upper()

        try:
            obabel_convert_command_mol = f'obabel {smi_file} -omol -h -O {new_name}_2d.mol --gen2d'
            subprocess.run(obabel_convert_command_mol, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            obabel_convert_command_mol2 = f'obabel {new_name}_2d.mol -omol -h -O {new_name}.mol --gen3d'
            subprocess.run(obabel_convert_command_mol2, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
           
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while generating conformation of {new_name}")
        finally:
            os.remove(smi_file)
            os.remove(f"{new_name}_2d.mol")
     
    pdb_files = [f for f in os.listdir(main_ligand_tmp) if f.endswith(".mol")]
    
    for pdb_file in tqdm(pdb_files, desc="Geometry optimization of main ligands"):
        obminimize_input = pdb_file
        obminimize_output = f"{os.path.splitext(pdb_file)[0]}_minimized.pdb"
         
        try:
            obminimize_command = f'obminimize -ff Ghemical -o "pdb" {obminimize_input} > {obminimize_output}'
            subprocess.run(obminimize_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)             
            
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while energy minimization of {obminimize_input}")
        finally:
            os.remove(obminimize_input)
  
def process_in_ligand_tmp_dir(ligand_tmp_dir, num_iterations):
    os.chdir(ligand_tmp_dir)    
    smi_files = [f for f in os.listdir(ligand_tmp_dir) if f.endswith(".smi")]
    
    for smi_file in tqdm(smi_files, desc="Generating 3D conformation of test ligands"):
        with open(smi_file, 'r') as file:
            parts = file.readline().split()[1]
            file_name = parts.split('/')[-1]
            new_name = file_name.split('.')[0].upper()

        try:
            obabel_convert_command_mol = f'obabel {smi_file} -omol -h -O {new_name}_2d.mol --gen2d'
            subprocess.run(obabel_convert_command_mol, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            obabel_convert_command_mol2 = f'obabel {new_name}_2d.mol -omol -h -O {new_name}.mol --gen3d'
            subprocess.run(obabel_convert_command_mol2, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
           
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while generating conformation of {new_name}")
        finally:
            os.remove(smi_file)
            os.remove(f"{new_name}_2d.mol")
     
    pdb_files = [f for f in os.listdir(ligand_tmp_dir) if f.endswith(".mol")]
    
    for pdb_file in tqdm(pdb_files, desc="Geometry optimization of test ligands"):
        obminimize_input = pdb_file
        obminimize_output = f"{os.path.splitext(pdb_file)[0]}_minimized.pdb"
         
        try:
            obminimize_command = f'obminimize -ff Ghemical -o "pdb" {obminimize_input} > {obminimize_output}'
            subprocess.run(obminimize_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)             
            
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while energy minimization of {obminimize_input}")
        finally:
            os.remove(obminimize_input)
     
def process_ligand_link(ligand_link_path, ligand_dir):
    for link_path in ligand_link_path:
        if os.path.exists(link_path):
            with open(link_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("http"):
                        ligand_url = line
                        try:
                            # Download the ligand file
                            ligand_file = os.path.join(ligand_dir, os.path.basename(ligand_url))
                            urllib.request.urlretrieve(ligand_url, ligand_file)
                            
                            # Membaca semua baris ke dalam list
                            with open(ligand_file, 'r') as file:
                                lines = file.readlines()

                            # Menulis kembali isi list ke file, mengabaikan baris pertama
                            with open(ligand_file, 'w') as file:
                                file.writelines(lines[1:])
                        except Exception as e:
                            print(f"An error occurred while downloading the url ligand: {ligand_url}")
                            print(f"Error details: {e}")
                        
def process_docking_ligand_test(main_ligand_tmp, ligand_tmp_dir, output_model_dir, counter, size_x, size_y, size_z, num_modes, exhaustiveness, cpu, num_iterations, csv_path, receptor_pdbqt_path, prepare_ligand_script, directory, vina):
    ligand_files = glob.glob(os.path.join(ligand_tmp_dir, '*.pdb'))
    num_ligand_files = len(ligand_files)
    
    main_ligand_files = glob.glob(os.path.join(main_ligand_tmp, '*.pdb'))
    num_main_files = len(main_ligand_files)
    print(f"  Total main ligand = {num_main_files}")
    print(f"  Total test ligand = {num_ligand_files}")
    print("  Using docking parameter from reference ligand")
    
    for main_ligand in main_ligand_files:
        main_name = os.path.splitext(os.path.basename(main_ligand))[0]
        main_dest_path = os.path.join(".", main_name + ".pdb")
        print(f"  Main ligand: {main_name}")
        for pdb_file in tqdm(ligand_files, desc="  Docking & save energy"):
            counter["value"] += 1
            filename = os.path.splitext(os.path.basename(pdb_file))[0]

            # Construct the destination path
            dest_path = os.path.join(".", filename + ".pdb")

            try:
                if not os.path.exists(dest_path):
                    # Copy pdb_file to the current directory
                    shutil.copy(main_ligand, main_dest_path)
                    shutil.copy(pdb_file, dest_path)

                # Convert *_lig.pdb to pdbqt
                main_pdbqt = f"{main_name}.pdbqt"
                run_command(f'{prepare_ligand_script} -l {main_ligand} -o {main_pdbqt}')
                
                pdbqt_file = f"{filename}.pdbqt"
                run_command(f'{prepare_ligand_script} -l {pdb_file} -o {pdbqt_file}')

                # Run Vina            
                output_name = f"output_{main_name}_{filename}"                            
                run_command(f'{vina} --receptor {receptor_pdbqt_path} --ligand {main_pdbqt} {pdbqt_file} --config config.txt --out {output_name}.pdbqt')
                
                # Move files
                os.rename(f"{output_name}.pdbqt", os.path.join(output_model_dir, f"{output_name}.pdbqt"))

                #Remove the original pdb_file and pdbqt_file
                # os.remove(main_dest_path)
                # os.remove(dest_path)
                # os.remove(main_pdbqt)
                # os.remove(pdbqt_file)

                # Save to CSV
                read_energy_and_save_to_csv(counter, output_model_dir, output_name, csv_path)            
                
            except Exception as e:
                print(f"An error occurred writing energy of {filename.replace('_minimized', '')}")

            except Exception as e:
                if os.path.exists(main_dest_path):
                    os.remove(main_dest_path)
                if os.path.exists(os.path.join(directory, f'{main_name}.pdbqt')):
                    os.remove(os.path.join(directory, f'{main_name}.pdbqt'))
                if os.path.exists(dest_path):
                    os.remove(dest_path)
                if os.path.exists(os.path.join(directory, f'{filename}.pdbqt')):
                    os.remove(os.path.join(directory, f'{filename}.pdbqt'))

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
    current_directory = os.getcwd()
    os.chdir(current_directory)
          
    # Membaca isi file 'config_tmp.txt'
    with open('config_la2vina.txt', 'r') as config_file:
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

    # Sekarang, Anda memiliki variabel-variabel dalam kamus config_variables yang dapat digunakan dalam skrip 'lavina.py'
    size_x = config_variables['size_x']
    size_y = config_variables['size_y']
    size_z = config_variables['size_z']
    num_modes = config_variables['num_modes']
    exhaustiveness = config_variables['exhaustiveness']
    cpu = config_variables['cpu']
    num_iterations = config_variables['num_iterations']    
    ligand_dir = os.path.join(current_directory, "ligand_input")
    ligand_tmp_dir = os.path.join(current_directory, "ligand_tmp")
    main_ligand_dir = os.path.join(current_directory, "main_ligand")
    main_ligand_tmp = os.path.join(current_directory, "main_ligand_tmp")    
    ligand_link_files = glob.glob(os.path.join(ligand_dir, "*ligand_link.txt"))
    
    mgl_directory = config_variables['mgl_directory']
    prepare_ligand = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_ligand4.py")
    prepare_receptor = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_receptor4.py")
    prepare_gpf = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_gpf4.py")
    prepare_dpf = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_dpf4.py")
    prepare_lowest_energy = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "write_lowest_energy_ligand.py")
    prepare_summarize_result = os.path.join(mgl_directory, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "summarize_results4.py")
    
    if ligand_link_files:
        ligand_link_path = ligand_link_files
        
    else:
        ligand_link_path = []

    if os.path.exists(config_variables['vina']):
        vina = config_variables['vina']
    else:
        vina = "vina"
   
    prepare_ligand_script = os.path.join(mgltools_dir, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py")
    prepare_receptor_script = os.path.join(mgltools_dir, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py")
    ligand_tmp_dir = os.path.join(current_directory, 'ligand_tmp')
    main_ligand_tmp = os.path.join(current_directory, 'main_ligand_tmp')
    output_dir = os.path.join(current_directory, 'output')

    print("Removing unnecessary files and directories in the output target")

    # Create ligand_tmp directory if it doesn't exist
    if os.path.exists(ligand_tmp_dir):
        for item in os.listdir(ligand_tmp_dir):
            item_path = os.path.join(ligand_tmp_dir, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
    else:
        os.makedirs(ligand_tmp_dir, exist_ok=True)
        
    # Create main_ligand_tmp directory if it doesn't exist
    if os.path.exists(main_ligand_tmp):
        for item in os.listdir(main_ligand_tmp):
            item_path = os.path.join(main_ligand_tmp, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
    else:
        os.makedirs(main_ligand_tmp, exist_ok=True)

    # Create output directory if it doesn't exist
    if os.path.exists(output_dir):
        for item in os.listdir(output_dir):
            item_path = os.path.join(output_dir, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
    else:
        os.makedirs(output_dir, exist_ok=True)

    print("Getting ligands from online sources")
    process_ligand_link(ligand_link_path, ligand_dir)

    print("Preparing the test ligands")
    process_in_ligand_dir(ligand_dir, ligand_tmp_dir)
    process_in_ligand_tmp_dir(ligand_tmp_dir, num_iterations)
    process_in_main_ligand(main_ligand_dir, main_ligand_tmp)
    process_in_main_ligand_tmp(main_ligand_tmp, num_iterations)

    
    os.chdir(current_directory)
    
    for directory in glob.glob('model*'):
        if os.path.isdir(directory):
            counter = {"value": 0}
            print(f"\nProcessing {directory.upper()}")
            os.chdir(directory)
            
            receptor_pdb_path = f"{os.path.basename(directory)}_rec.pdb"
            receptor_pdbqt_path = f"{os.path.basename(directory)}_rec.pdbqt"
            lig_pdb_path = f"{os.path.basename(directory)}_lig.pdb"
            lig_pdbqt_path = f"{os.path.basename(directory)}_lig.pdbqt"
            lig_name = f"{os.path.basename(directory)}_lig"
            model_dir = f"output_{os.path.basename(directory)}"
            output_model_dir = os.path.join(output_dir, model_dir)
            os.makedirs(output_model_dir, exist_ok=True)
            csv_path = os.path.join(output_model_dir, 'energy_summary.csv')
                    
            
            if os.path.exists(receptor_pdb_path) or os.path.exists(lig_pdb_path) or os.path.exists(receptor_pdbqt_path) or os.path.exists(lig_pdbqt_path):
                print("Docking reference ligand:")
                print(f"  Receptor: {receptor_pdbqt_path}")
                print(f"  Ligand: {lig_pdbqt_path}")
                process_docking_ligand_native(counter, lig_pdb_path, size_x, size_y, size_z, num_modes, exhaustiveness, cpu, num_iterations, csv_path, prepare_receptor_script, receptor_pdb_path, receptor_pdbqt_path, prepare_ligand_script, lig_pdbqt_path, directory, vina, output_model_dir)
                print("Docking test ligands: ")
                print(f"  Receptor: {receptor_pdbqt_path}")        
                process_docking_ligand_test(main_ligand_tmp, ligand_tmp_dir, output_model_dir, counter, size_x, size_y, size_z, num_modes, exhaustiveness, cpu, num_iterations, csv_path, receptor_pdbqt_path, prepare_ligand_script, directory, vina)
            else:
                print(f"Skipping docking: receptor or reference ligand in {directory} not found.")
            
            os.chdir(current_directory)
            print(f"Success: {directory}")

    # Print developer's note, contact, and citation list
    print_dev(developer_note, developer_contact, citation_list)

if __name__ == "__main__":
  
    current_directory = (os.path.join(os.getcwd(), "LADOCK_la2vina"))
    
    if os.path.exists(current_directory):
        os.chdir(current_directory)
        main() 
    else:
        print("Your job directory (LADOCK_la2vina) is not ready. Please create it using:")
        print("ladock --create la2vina")
