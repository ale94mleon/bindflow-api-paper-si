from bindflow.preparation.solvent import index_for_soluble_system
from pathlib import Path
from glob import glob

index_files = glob('mmpbsa/*/*/input/complex/index.ndx')

for index_file in index_files:
    index_file = Path(index_file)
    parent_dir = index_file.parent
    print(index_file)
    index_for_soluble_system(
        configuration_file=parent_dir/'complex.gro',
        ndxout=index_file,
        ligand_name='LIG',
        host_name='Protein or resname COF',
        load_dependencies=['export GMX_MAXBACKUP=-1']
    )