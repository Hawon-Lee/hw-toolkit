#!/bin/bash

# echo "CC(C)[C@H]1C(=O)NCCCCCCCCCCCC(=O)N[C@H]([C@H](CC(=O)N1)O)Cc2ccccc2" > ligand_of_interest.smi

# obabel -ismi ligand_of_interest.smi \
# -osdf -O ligand_of_interest.sdf \
# --gen3d --best -p 7.4

# mk_prepare_ligand.py -i ligand_of_interest.sdf \
# -o ligand_of_interest.pdbqt

/home/tech/Hawon/DockingTool/autodock_vina_1_1_2_linux_x86/bin/vina \
--receptor /home/tech/Hawon/docking_workspace/processed/6qdb.pdbqt \
--ligand ligand_of_interest.pdbqt \
--config /home/tech/Hawon/docking_workspace/processed/6qdb.txt \
--exhaustiveness 12 \
--out ./singlerun.sdf

rm ligand_of_interest.sdf ligand_of_interest.pdbqt ligand_of_interest.smi