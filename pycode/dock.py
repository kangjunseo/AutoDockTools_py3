import os
import subprocess
import requests
import time
from urllib.parse import quote

from rdkit import Chem
from rdkit.Chem import SaltRemover


current_dir = os.getcwd()
# 사용자 에이전트 헤더 추가
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}

# PubChem 요청 시 서버 바쁨 오류 처리 및 재시도 로직 추가
def get_smiles(INCI, retries=5, delay=5):

    # Helper function to make requests with headers
    def make_request(identifier, namespace):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{namespace}/{quote(identifier)}/property/CanonicalSMILES/TXT"
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.text.strip()


    # INCI Name으로 검색
    for attempt in range(retries):
        try:
            smiles = make_request(INCI, 'name')
            if smiles: return smiles
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                print(f"Not found: {INCI} (404 error)")
                break
            elif e.response.status_code == 400:
                print(f"Bad request: {INCI} (400 error)")
                break
            elif 'PUGREST.ServerBusy' in str(e):
                print(f"Server busy: {e}. Retrying ({attempt + 1}/{retries}) in {delay} seconds...")
                time.sleep(delay)
            else:
                raise e
        except (requests.ConnectionError, requests.Timeout) as e:
            print(f"Connection error: {e}. Retrying ({attempt + 1}/{retries}) in {delay} seconds...")
            time.sleep(delay)
    return None
    #INCI_SMILES_df['SMILES'] = INCI_SMILES_df.apply(lambda row: get_smiles(row["INCI"]), axis=1)


def salt_remove(smiles):
    remover = SaltRemover.SaltRemover()
    
    # Mol 객체 생성
    mol = Chem.MolFromSmiles(smiles)
    # 염 제거
    clean_mol = remover.StripMol(mol, dontRemoveEverything=True)
    # SMILES로 변환
    clean_smiles = Chem.MolToSmiles(clean_mol)
    
    return clean_smiles


def is_organic(smiles):
    metals = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            has_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
            is_inorganic = any(atom.GetSymbol() in metals for atom in mol.GetAtoms())
            return has_carbon and not is_inorganic
    except:
        return False
    return False

def create_pdb(name, smiles, output_dir):
    
    def sanitize_filename(filename): return filename.replace(' ', '_')

    output_file = os.path.join(output_dir, f"{name}.pdb")
    output_file = sanitize_filename(output_file)

    command = f'obabel -:"{smiles}" -O {output_file} --gen3d'

    try:
        subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        print(f"Converted {name} to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error converting {name}: {e}")


def name2pdb(name, output_dir):
    res = get_smiles(name)
    res = salt_remove(res)
    if is_organic(res): 
        create_pdb(name, res, output_dir)
    else:
        print("Not organic compound")

## prepare pdbqt
#"python3 /home/kjs/Downloads/AutoDockTools_py3/AutoDockTools/Utilities24/prepare_ligand4.py -h"

def prepare_pdbqt(pdb_file, output_file):
    output_file = os.path.join(output_file, os.path.basename(pdb_file).replace('.pdb', '.pdbqt'))
    try:
        prepare_ligand_path = "../AutoDockTools/Utilities24/prepare_ligand4.py"
        command = f'python3 {prepare_ligand_path} -l {pdb_file} -o {output_file} -A hydrogens'
        result = subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout_output = result.stdout.decode('utf-8')
        stderr_output = result.stderr.decode('utf-8')

        print(f"Processed {pdb_file} to {output_file}")
        print("Standard Output:")
        print(stdout_output)

        if stderr_output:
            print("Standard Error:")
            print(stderr_output)

    except subprocess.CalledProcessError as e:
        print(f"Error converting {pdb_file} to PDBQT: {e.stderr.decode('utf-8')}")
        
#name = 'RTX'
#prepare_pdbqt(os.path.join(current_dir, f'{name}.pdb'))

## Rigid Docking
def dock_ligand(receptor, ligand, output):
    vina_executable = "../../autodock_vina_1_1_2_linux_x86/bin/vina"  # AutoDock Vina 실행 파일 경로 설정
    try:
        command = [
            vina_executable,
            '--receptor', receptor,
            '--ligand', ligand,
            '--out', output,
            '--log', output.replace('.pdbqt', '.log'),
            '--center_x', '-23.075', '--center_y', '-2.183', '--center_z', '-7.899',  # 이 값을 조정하여 도킹 박스의 중심을 설정하십시오.
            '--size_x', '30', '--size_y', '30', '--size_z', '30'  # 이 값을 조정하여 도킹 박스의 크기를 설정하십시오.
        ]
        result = subprocess.run(command, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        print(f"Docked {ligand} to {receptor}, output saved to {output}")
    except subprocess.CalledProcessError as e:
        print(f"Error docking {ligand} to {receptor}: {e.stderr.decode('utf-8')}")


def main(receptor, name):
    output_dir = os.path.join(current_dir, f'dock_result_{name}')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        os.makedirs(output_path, exist_ok=True)
    
    name2pdb(name, output_dir) # convert name to SMILES, and SMILES to 3D pdb
    
    prepare_pdbqt(os.path.join(output_dir, f'{name}.pdb'), output_dir) # convert pdb to pdbqt
    
    dock_ligand(receptor, os.path.join(output_dir, f'{name}.pdbqt'), os.path.join(output_dir, f'r_docked_{name}.pdbqt'))

if __name__ == '__main__':
    import argparse
    import sys
    start_time = time.time()
    parser = argparse.ArgumentParser(description="(ex. python3 ./dock.py -r receptor.pdbqt -l ligand_name)")
    parser.add_argument("-r", "--receptor", type=str, required=True, help="receptor pdbqt file")
    parser.add_argument("-l", "--ligand_name", type=str, required=True, help="name of ligand(will be searched at pubchem)")
    #parser.add_argument("-o", "--output", type=str, required=True, help="output path, must be .pdbqt")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    main(args.receptor, args.ligand_name)

    end_time = time.time()
    elapsed_time = end_time - start_time
    log_path = os.path.join(current_dir, f'dock_result_{args.ligand_name}', f'r_docked_{args.ligand_name}.log')
    try:
        subprocess.call(['cat', log_path])
        print(f"Execution time: {elapsed_time:.3f} seconds")
    except FileNotFoundError:
        print(f"Log not Found: {log_path}, Docking is not completely done.")
    