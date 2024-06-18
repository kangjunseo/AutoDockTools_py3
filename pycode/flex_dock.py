import os
import subprocess
import time

import utils, dock

current_dir = os.getcwd()

##main 만들고 pdb pdbqt flexdocking 해야함
def main(r_receptor, f_receptor, ligand_pdb, output):
    name = os.path.basename(ligand_pdb)[:-4]
    output_dir = os.path.join(current_dir, output, f'f_dock_result_{name}')
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        dock.prepare_pdbqt(ligand_pdb, output_dir)  # convert pdb to pdbqt
        dock.dock_ligand_flex(r_receptor, f_receptor, os.path.join(output_dir, f'{name}.pdbqt'), os.path.join(output_dir, f'f_docked_{name}.pdbqt'))

    except Exception as e:
        print(f"Error: {e}")
        exit(1)
        


if __name__ == '__main__':
    import argparse
    import sys
    start_time = time.time()
    parser = argparse.ArgumentParser(description="(ex. python3 ./flex_dock.py -r rigid_receptor.pdbqt -f flex_receptor.pdbqt -l ligand_pdb)")
    parser.add_argument("-r", "--r_receptor", type=str, required=True, help="rigid receptor pdbqt file")
    parser.add_argument("-f", "--f_receptor", type=str, required=True, help="flex receptor pdbqt file")
    parser.add_argument("-l", "--ligand_pdb", type=str, required=True, help="pdb file of ligand")
    parser.add_argument("-o", "--output", type=str, required=True, help="output directory(please use relative path)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    args.ligand_pdb = utils.sanitize_filename(args.ligand_pdb)
    main(args.r_receptor, args.f_receptor, args.ligand_pdb, args.output)

    end_time = time.time()
    elapsed_time = end_time - start_time
    op = os.path.join(args.output if args.output else current_dir, f'f_dock_result_{os.path.basename(args.ligand_pdb)[:-4]}')
    log_path = os.path.join(op, f'f_docked_{os.path.basename(args.ligand_pdb)[:-4]}.log')
    try:
        subprocess.call(['cat', log_path])
        print(f"Execution time: {elapsed_time:.3f} seconds")
    except FileNotFoundError:
        print(f"Log not Found: {log_path}, Docking is not completely done.")
    