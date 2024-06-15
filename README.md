# AutoDockTools_py3
Python3 translation of AutoDockTools
+one click docking

## Prerequisites(IMPORTANT)

**you must have excutable vina file**
download at https://vina.scripps.edu/downloads/  
or https://github.com/ccsb-scripps/AutoDock-Vina

major dependency of this repository is:  
- python=3.12.3=h996f2a0_1
- pandas=2.2.1=py312h526ad5a_0
- babel=2.11.0=py312h06a4308_0
- pip:
  - rdkit==2023.9.6
  - matplotlib==3.9.0
 
OR  
See 'environment.yml' for the required packeges.

## Installation

```sh
git clone https://github.com/kangjunseo/AutoDockTools_py3.git
cd AutoDockTools_py3
conda env create -f environment.yml
conda activate MD
```


## Usage

- One receptor + ligand_name  

```shell

python3 pycode/dock.py -r {path/receptor.pdbqt} -l {ligand_name}

```



## Known Errors (not critical)

![image](https://github.com/kangjunseo/AutoDockTools_py3/assets/88201041/cf43adcc-3432-4aec-b71e-488ff17fe514)

This error is minor error because of incomplete converting ADT py files from py2 to py3.
