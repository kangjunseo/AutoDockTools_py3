import os
import subprocess
import requests
import time
#from timeout_decorator import timeout, TimeoutError
from urllib.parse import quote

from rdkit import Chem
from rdkit.Chem import SaltRemover
import pandas as pd

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

def sanitize_filename(filename): return filename.replace(' ', '_')

def create_pdb(name, smiles, output_dir):
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
    if res == None: return 1
    res = salt_remove(res)
    if is_organic(res): 
        create_pdb(name, res, output_dir)
    else:
        print("Not organic compound")

def smiles2pdb(name, smiles, output_dir):
    res = salt_remove(smiles)
    if is_organic(res): 
        create_pdb(name, res, output_dir)
    else:
        print("Not organic compound")


def parse_log_file(log_file):
    """Parse the log file to extract the compound name and affinity value."""
    with open(log_file, 'r') as file:
        lines = file.readlines()
    
    compound_name = os.path.basename(log_file).replace('r_docked_', '').replace('.log', '')
    affinity = None
    
    for line in lines:
        if line.strip().startswith('1'):
            parts = line.split()
            if len(parts) > 1:
                try:
                    affinity = float(parts[1])
                except ValueError:
                    continue
            break
    
    return compound_name, affinity


def process_log_files(log_dir):
    """Process all log files in the specified directory and return a DataFrame."""
    data = []

    for root, dirs, files in os.walk(log_dir):
        for log_file in files:
            if log_file.endswith('.log'):
                log_path = os.path.join(root, log_file)
                compound_name, affinity = parse_log_file(log_path)
                if affinity is not None:
                    data.append((compound_name, affinity))
    
    df = pd.DataFrame(data, columns=['INCI', 'Affinity'])
    return df