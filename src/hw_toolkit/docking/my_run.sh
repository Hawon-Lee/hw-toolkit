#!/bin/bash

##########################################################################################################

# python my_receptor_preparation.py \
# --receptor /home/tech/Hawon/docking_workspace/raw/receptors/6qcb.pdb \
# --out ../processed/6qcb.pdbqt \
# --autobox_path ../processed/6qcb_box.txt \
# --autobox_size 25

##########################################################################################################

# python my_ligand_preparation.py \
# --ligand /home/tech/Hawon/docking_workspace/processed/deg2mol_mw210_random.csv \
# --out_dir ../processed/ \
# --create_batch

##########################################################################################################

# python my_docking.py \
# -r /home/tech/Hawon/docking_workspace/processed/6qdb.pdbqt \
# -l /home/tech/Hawon/docking_workspace/processed/04_erlo_renamed/batch.txt \
# --autobox /home/tech/Hawon/docking_workspace/processed/6qdb.txt \
# --out_dir /home/tech/Hawon/docking_workspace/result \
# --exhaustiveness 12 \
# --project_id "erlo_to_6qdb"

##########################################################################################################

python my_docking_multi.py \
-r ../processed/6qcb.pdbqt \
-l /home/tech/Hawon/docking_workspace/processed/04_erlo_renamed/batch.txt \
--autobox ../processed/6qcb_box.txt \
--out_dir /home/tech/Hawon/docking_workspace/result \
--exhaustiveness 12 \
--num_processes 3 \
--project_id "erlo_to_6qcb"