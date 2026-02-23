# pocket_sucos_code
Code for calculating sucos-pocket similarity from an arbitrary protein-ligand complex

To run:
```
bash 01_Sh_run_foldseek.sh $RECEPTOR_PDB $LIGAND_SDF
```

Known Issues:
-Foldseek alignment fails for chains with nonstandard residues
-Code only tested on single chain proteins
