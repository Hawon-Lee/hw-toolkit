import os
import numpy as np
import pandas as pd
import argparse
import inquirer
import subprocess
import random
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem 

def smiles_into_pdbqt(smiles: str, mol_id: str, output_dir: str) -> str:
    command = ["scrub.py", smiles, "-o", f"{output_dir}/{mol_id}.sdf"]
    subprocess.run(command)
    
    command = [
        "obabel",
        "-isdf", f"{output_dir}/{mol_id}.sdf",
        "-osdf", "-O", f"{output_dir}/{mol_id}_p.sdf",
        "-p", "7.4",
        "--minimize"]
    subprocess.run(command)

    command = ["mk_prepare_ligand.py", 
        "-i", f"{output_dir}/{mol_id}_p.sdf",
        "-o", f"{output_dir}/{mol_id}.pdbqt"]
    subprocess.run(command)

    os.remove(f"{output_dir}/{mol_id}.sdf")
    os.remove(f"{output_dir}/{mol_id}_p.sdf")

    if os.path.exists(f"{output_dir}/{mol_id}.pdbqt"):
        with open(f"{output_dir}/{mol_id}.pdbqt", "r") as fp:
            pdbqt = fp.read()
    else:
        print(f"Failed to generate pdbqt from mol id {mol_id}")
        return None
    
    return f"{output_dir}/{mol_id}.pdbqt"


def validate_mol2_coordinates(mol2_file):
    """MOL2 파일의 좌표가 유효한지 확인"""
    try:
        with open(mol2_file, "r") as f:
            lines = f.readlines()

        in_atom_section = False
        atom_coords = []

        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
                continue

            if in_atom_section and line.strip():
                parts = line.strip().split()
                if len(parts) >= 5:
                    try:
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        atom_coords.append([x, y, z])
                    except ValueError:
                        continue

        # 모든 좌표가 0인지 확인
        if len(atom_coords) == 0:
            return False

        non_zero_coords = [coord for coord in atom_coords if coord != [0.0, 0.0, 0.0]]
        return (
            len(non_zero_coords) > len(atom_coords) * 0.5
        )  # 50% 이상이 0이 아닌 좌표여야 함

    except Exception:
        return False


def parse_arguments():

    parser = argparse.ArgumentParser(description="ligand preparation for docking.")

    parser.add_argument(
        "--ligand_csv",
        "-l",
        type=str,
        required=True,
        help="mol_id, smiles 로 구성된 csv 파일 경로",
    )
    parser.add_argument(
        "--out_dir",
        "-o",
        type=str,
        required=True,
        help="전처리 결과 sdf 파일들을 저장할 폴더",
    )
    parser.add_argument(
        "--create_batch",
        "-b",
        action="store_true",
        help="(Optional) 전처리된 smiles의 경로들을 txt 파일로 기록합니다.",
    )
    parser.add_argument("--tool", "-t", type=str, help="mgl or mk")

    args = parser.parse_args()
    return args


def ligand_preparation(out_dir, ligand_csv, mol_id_col, smiles_col, create_batch=False, prep_type="mgl"):
    if isinstance(ligand_csv, str):
        project_dir = os.path.join(out_dir, os.path.basename(ligand_csv).split(".")[0])
    else:
        from datetime import datetime
        project_dir = f"{out_dir}/{datetime.now()}"

    os.makedirs(project_dir, exist_ok=True)

    if isinstance(ligand_csv, str):
        df = pd.read_csv(ligand_csv)
    else:
        df = ligand_csv

    # mol_ids, smiles = df.iloc[:, 0].values, df.iloc[:, 1].values
    mol_ids, smiles = df.loc[:, mol_id_col].values, df.loc[:, smiles_col].values

    pdbqt_files = []
    blacklist = []

    curr_dir = os.getcwd()
    os.chdir(project_dir)
    for mol_id, smile in zip(mol_ids, smiles):
        # temp_smi_file = f"{project_dir}/{mol_id}.smi"
        # temp_mol2_file = f"{project_dir}/{mol_id}.mol2"
        # pdbqt_file = f"{project_dir}/{mol_id}.pdbqt"
        temp_smi_file = f"{mol_id}.smi"
        temp_mol2_file = f"{mol_id}.mol2"
        temp_sdf_file = f"{mol_id}.sdf"
        pdbqt_file = f"{mol_id}.pdbqt"

        # smiles를 smi 파일로 저장
        with open(temp_smi_file, "w") as fp:
            fp.write(f"{smile}")

        # smi gen3d
        command = (
            f"obabel -ismi {temp_smi_file} -omol2 -O {temp_mol2_file} --gen3d -p 7.4"
        )
        try:
            result = subprocess.run(
                command.split(),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

            # 성공적으로 실행된 경우

            if not validate_mol2_coordinates(temp_mol2_file):
                os.remove(temp_mol2_file)
                m = Chem.MolFromSmiles(smile)
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m, useBasicKnowledge=False)
                AllChem.MMFFOptimizeMolecule(m, maxIters=150)
                Chem.MolToMolFile(m, temp_sdf_file)
                command = (
                    f"obabel -isdf {temp_sdf_file} -omol2 -O {temp_mol2_file} -p 7.4"
                )
                result = subprocess.run(
                    command.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                if not validate_mol2_coordinates(temp_mol2_file):
                    print(
                        f"Warning: Obabel -> Rdkit fallback으로도 smi -> mol2 변환 실패: {mol_id}"
                    )
                    blacklist.append(mol_id)
            print(f"SMI -> MOL2 변환이 성공적으로 완료되었습니다: {mol_id}")
            # 여기에 성공 시 수행할 작업 추가
        except Exception as e:
            print(f"OpenBabel 실행 중 오류: {e}")

            # 통합 오류 처리 로직

        # sdf -> pdbqt
        if prep_type.lower() == "mk":
            command = f"mk_prepare_ligand.py -i {temp_mol2_file} \
                -o {pdbqt_file}"
            subprocess.run(command.split())
        elif prep_type.lower() == "mgl":
            command = f"python2 \
                /home/tech/anaconda3/envs/autodock/bin/prepare_ligand4.py \
                -l {temp_mol2_file} \
                -o {pdbqt_file}"
            subprocess.run(command.split())

        # temp 파일 삭제
        os.remove(temp_smi_file)
        os.remove(temp_mol2_file)

        if not os.path.exists(pdbqt_file):
            print(
                f"warning: {mol_id} was not converted successfully into .pdbqt - added to blacklist."
            )
            # blacklist.append(mol_id)
            with open("./blacklist.txt", "a") as f:
                f.write(f"{mol_id}\n")
            continue
        with open(pdbqt_file, "r") as fp:
            if fp.read() is not None:
                pdbqt_files.append(os.path.abspath(pdbqt_file))

    if create_batch:
        with open(f"{project_dir}/batch.txt", "w") as fp:
            fp.write("\n".join(pdbqt_files))
            print(f"batch file was created at {project_dir}/batch.txt")


def main():
    args = parse_arguments()
    ligand_preparation(args.out_dir, args.ligand_csv, args.create_batch)


if __name__ == "__main__":
    main()
