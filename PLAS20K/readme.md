# Input extended_PLAS20K.csv:
#### usePDBnames.go to extract all the pdb_id of protein-ligand complex used in PLAS20K (stored in PLAS20K_pdb_ids.txt)
#### use batch_download.sh to grab all the pdb files from RSCB
#### use splitPDB.go to seperate proteins and ligands (output two pdb files for proteins and ligands)
#### use convert_pdb_to_mol2.sh (calls Open Babel) to convert all pdb files to mol2 files.


#### mol2 files were send using Google Drive.
