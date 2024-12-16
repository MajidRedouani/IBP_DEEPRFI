import os
import requests
import time
import subprocess

# Set up directories
os.makedirs('pdb_files', exist_ok=True)
missing_pdb_path = 'Proteins_without_PDB.txt'
results_path = 'results_DeepFRI.tsv'

# Read UniProt IDs from the fasta file
fasta_path = 'proteins.fasta'
with open(fasta_path, 'r') as fasta_file:
    protein_ids = [line.strip()[1:] for line in fasta_file if line.startswith(">")]

# Function to download and modify AlphaFold PDB file
def download_and_modify_pdb(uniprot_id):
    pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(pdb_url)
    
    if response.status_code == 200:
        file_path = f"pdb_files/{uniprot_id}.pdb"
        with open(file_path, 'wb') as file:
            file.write(response.content)
        
        # Modify the HEADER line to add a dummy PDB ID (e.g., "0000")
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        with open(file_path, 'w') as file:
            for line in lines:
                if line.startswith("HEADER"):
                    file.write(f"HEADER    DUMMY PDB ID 0000\n")
                else:
                    file.write(line)
        
        print(f"Downloaded and modified {uniprot_id}.pdb from AlphaFold")
        return True
    else:
        print(f"Failed to download PDB for UniProt ID {uniprot_id} from AlphaFold")
        return False

# Function to run DeepFRI on a modified PDB file and save results
def run_deepfri(pdb_file, uniprot_id):
    
    try:
        result = subprocess.run(
            ["python", "predict.py", "-pdb", pdb_file, "-ont", "mf", "-v"],
            capture_output=True, text=True, check=True
        )
        
        # Parse and save results
        with open(results_path, 'a') as results_file:
            for line in result.stdout.splitlines():
                if line.startswith("query_prot"):
                    results_file.write(f"{uniprot_id}\t{line}\n")
        
        print(f"Prediction completed for {uniprot_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error running DeepFRI on {pdb_file}: {e}")
        
    try:
        result = subprocess.run(
            ["python", "predict.py", "-pdb", pdb_file, "-ont", "cc", "-v"],
            capture_output=True, text=True, check=True
        )
        
        # Parse and save results
        with open(results_path, 'a') as results_file:
            for line in result.stdout.splitlines():
                if line.startswith("query_prot"):
                    results_file.write(f"{uniprot_id}\t{line}\n")
        
        print(f"Prediction completed for {uniprot_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error running DeepFRI on {pdb_file}: {e}")
    try:
        result = subprocess.run(
            ["python", "predict.py", "-pdb", pdb_file, "-ont", "bp", "-v"],
            capture_output=True, text=True, check=True
        )
        
        # Parse and save results
        with open(results_path, 'a') as results_file:
            for line in result.stdout.splitlines():
                if line.startswith("query_prot"):
                    results_file.write(f"{uniprot_id}\t{line}\n")
        
        print(f"Prediction completed for {uniprot_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error running DeepFRI on {pdb_file}: {e}")

# Main process to download and predict
with open(missing_pdb_path, 'w') as missing_file:
    for protein_id in protein_ids:
        if download_and_modify_pdb(protein_id):
            pdb_file_path = f"pdb_files/{protein_id}.pdb"
            run_deepfri(pdb_file_path, protein_id)
        else:
            missing_file.write(f"{protein_id}\n")
        time.sleep(1)  # Delay to avoid overloading the server
