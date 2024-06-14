import os
import subprocess
import time

import utils

current_dir = os.getcwd()

## prepare pdbqt
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
        os.makedirs(output_dir, exist_ok=True)
    
    utils.name2pdb(name, output_dir) # convert name to SMILES, and SMILES to 3D pdb
    
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
    