# pocket_sucos_code
Code for calculating sucos-pocket similarity from an arbitrary protein-ligand complex. Currently, similarity metrics are
calculated with respect to PLINDER systems.

### Setup the conda environment:
```
conda env create -f environment.yml
```

### Install the PLINDER database, and set the following environment variables:
```
export PLINDER_OFFLINE=true
export PLINDER_MOUNT=/path/to/your/plinder
export PLINDER_ITERATION=v2
export PLINDER_RELEASE=2024-06
```

### Create a custom foldseek database from PLINDER receptor structures.

Copy receptor.cif structures installed in the `$PLINDER_MOUNT/plinder/2024-06/v2/systems/`, to a new directory, then run
the `foldseek createdb` command as specified in the [Foldseek documantation](https://github.com/steineggerlab/foldseek)

### To run:

Edit the `FOLDSEEK_DB` variable in `01_Sh_run_foldseek.sh` to point to the path of tour custom foldseek database
```
bash 01_Sh_run_foldseek.sh $RECEPTOR_PDB $LIGAND_SDF $OUTDIR
```

### Known Issues:

-Foldseek alignment fails for chains with nonstandard residues

-Code only tested on single chain proteins
