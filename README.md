# pocket_sucos_code
Code for calculating sucos-pocket similarity from an arbitrary protein-ligand complex. Currently, similarity metrics are
calculated with respect to PLINDER systems.

Setup the conda environment:
```
conda env create -f environment.yml
```

Install the PLINDER database, and set the following environment variables:
```
export PLINDER_OFFLINE=true
export PLINDER_MOUNT=/path/to/your/plinder
export PLINDER_ITERATION=v2
export PLINDER_RELEASE=2024-06
```

To run:
```
bash 01_Sh_run_foldseek.sh $RECEPTOR_PDB $LIGAND_SDF
```

Known Issues:

-Foldseek alignment fails for chains with nonstandard residues

-Code only tested on single chain proteins
