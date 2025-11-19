import os
import pandas as pd

def p2rank_res_to_pdb(pred_csv="./1a30_protein.pdb_predictions.csv", src_pdb="1a30_protein.pdb", out_dir="./p2rank_pocket"):
    '''
    Extract atom records from p2rank prediction and export .pdb format
    '''

    os.makedirs(out_dir, exist_ok=True)

    # read p2rank results
    pred = pd.read_csv(pred_csv)
    pred.columns = [i.strip() for i in pred.columns]
    residue_ids = pred.residue_ids[0].split()

    # match surf_atom_ids from p2rank to the source .pdb file
    poc_res_ls = []
    with open(src_pdb, "r") as fp:
        pdb_lines = fp.read().splitlines()

        for l in pdb_lines:
            if not l.startswith("ATOM" or "HETATM"):
                continue
            
            for res_id in residue_ids:
                res, r_id = res_id.split("_")
                
                if l[21] == res and l[22:27].strip() == r_id.strip():
                    poc_res_ls.append(l)   

    # export
    f_name = os.path.basename(src_pdb).split("_")[0] + "_pocket.pdb"
    with open(f"{out_dir}/{f_name}", "w") as fp:
        fp.write("\n".join(poc_res_ls))   
    
    return 0

def p2rank_atom_to_pdb(pred_csv="./1a30_protein.pdb_predictions.csv", src_pdb="1a30_protein.pdb", out_dir="./p2rank_pocket"):
    '''
    Extract atom records from p2rank prediction and export .pdb format
    '''

    os.makedirs(out_dir, exist_ok=True)

    # read p2rank results
    pred = pd.read_csv(pred_csv)
    pred.columns = [i.strip() for i in pred.columns]
    surf_atom_ids = {int(i) for i in pred.surf_atom_ids[0].split()}

    # match surf_atom_ids from p2rank to the source .pdb file
    poc_atom_ls = []
    with open(src_pdb, "r") as fp:
        pdb_lines = fp.read().splitlines()

        for l in pdb_lines:
            if not l.startswith("ATOM" or "HETATM"):
                continue
            if int(l[6:11]) in surf_atom_ids:
                poc_atom_ls.append(l)

    # export
    f_name = os.path.basename(src_pdb).split("_")[0] + "_pocket_atoms.pdb"
    with open(f"{out_dir}/{f_name}", "w") as fp:
        fp.write("\n".join(poc_atom_ls))   
    
    return 0

if __name__ == "__main__":
    p2rank_res_to_pdb("1bcu_protein.pdb_predictions.csv", "1bcu_protein.pdb", "./p2rank_pocket")