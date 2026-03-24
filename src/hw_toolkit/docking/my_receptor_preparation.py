import os
import tempfile
import numpy as np
import argparse
import inquirer
import subprocess
from rdkit import Chem

from .box_calculator import DockingBoxCalculator

def parse_arguments():
    parser = argparse.ArgumentParser(description='protein preparation for docking.')
    
    parser.add_argument(
        '--receptor', '-r',
        type=str,
        required=True,
        help='단백질 PDB 파일 경로'
        )
    parser.add_argument(
        '--out_dir', '-o',
        type=str,
        required=True,
        help='결과물 파일을 저장할 디렉토리 경로'
        )
    parser.add_argument(
        '--autobox_path', '-a',
        type=str,
        default=None,
        help='리간드 정보로부터 자동 계산한 box 정보를 저장할 경로 (option)'
    )
    parser.add_argument(
        '--autobox_size', '-s',
        type=float,
        default=25,
        help="리간드 정보로부터 자동 계산할 box의 크기 ('fixed' 방식에만 사용)"
    )
    parser.add_argument(
        '--autobox_method',
        type=str,
        default='fixed',
        choices=['fixed', 'optimal'],
        help="자동 계산할 box 크기 결정 방식: 'fixed' (지정된 크기), 'optimal' (리간드 구조 기반 최적화)"
    )
    parser.add_argument(
        '--ligand_out', '-l',
        type=str,
        default=None,
        help='추출한 리간드를 저장할 SDF 파일 경로 (option)'
    )
    
    args = parser.parse_args()
    return args

def select_with_inquirer(tuple_options):
    # 표시용 문자열 생성 (resname, chainID, resID 포함)
    display_options = [f"resname '{tup[0]}' - chainID '{tup[1]}' - resID {tup[2]}" for tup in tuple_options]
    
    # 표시할 옵션과 원래 튜플을 매핑
    option_map = {display_options[i]: tuple_options[i] for i in range(len(tuple_options))}
    
    questions = [
        inquirer.List('chosen_option',
                      message="Box를 계산할 레퍼런스 리간드를 선택하세요",
                      choices=display_options,
                     )
    ]
    answers = inquirer.prompt(questions)
    
    # 사용자가 선택한 문자열을 원래 튜플로 변환
    selected_display = answers['chosen_option']
    selected_tuple = option_map[selected_display]
    
    return selected_tuple

def get_ligand_info(pdb_path, lig_name, chain2parse, res_id):
    with open(pdb_path, 'r') as f:
        pdb_lines = f.readlines()

    ligand_hetatm_lines = []
    ligand_atom_serials = set()

    # Extract HETATM lines for the specific ligand (resname, chain, resID 모두 일치해야 함)
    for l in pdb_lines:
        if l.startswith('HETATM') and l[17:20].strip() == lig_name.strip():
            if l[21] == chain2parse:
                # residue sequence number 확인 (columns 23-26, 0-indexed: 22:26)
                line_res_id = int(l[22:26].strip())
                if line_res_id == res_id:
                    ligand_hetatm_lines.append(l)
                    ligand_atom_serials.add(int(l[6:11]))

    if not ligand_hetatm_lines:
        return None, None

    # Extract CONECT records relevant to our ligand
    ligand_conect_lines = []
    for l in pdb_lines:
        if l.startswith('CONECT'):
            try:
                tokens = l.split()
                # Check if the first atom in the CONECT record is part of our ligand
                if int(tokens[1]) in ligand_atom_serials:
                    ligand_conect_lines.append(l)
            except (ValueError, IndexError):
                continue

    # Create a PDB block for the ligand
    ligand_pdb_block = "".join(ligand_hetatm_lines) + "".join(ligand_conect_lines)

    # Load molecule from PDB block. RDKit will use CONECT records for bonding.
    mol = Chem.MolFromPDBBlock(ligand_pdb_block, removeHs=False, sanitize=True)

    if mol is None:
        print(f"Warning: RDKit could not process the ligand {lig_name} (chain {chain2parse}, resID {res_id}).")
        # Fallback to calculate center from coordinates
        coords = []
        for l in ligand_hetatm_lines:
                x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
                coords.append([x,y,z])
        center = np.array(coords).mean(axis=0)
        return center, None
    
    # Calculate ligand center
    try:
        conformer = mol.GetConformer()
        positions = conformer.GetPositions()
        ligand_center = positions.mean(axis=0)
    except Exception as e:
        print(f"Could not calculate ligand center for {lig_name} (chain {chain2parse}, resID {res_id}): {e}")
        return None, None


    # Generate SDF block
    mol.SetProp("_Name", lig_name.strip())
    sdf_text = Chem.MolToMolBlock(mol)

    return ligand_center, sdf_text

def protein_preparation(receptor, ligand_out=None, autobox_path=None, autobox_method=None, autobox_size=None, out_dir=None):
    basename = os.path.basename(receptor).replace(".pdb", "")
    BUFFERS = ['GOL', 'PEG', 'EDO', 'SO4', 'BIS', 'HOH']
    METALS = ['ZN', 'MG', 'MN', 'CA', 'FE', 'CO', 'NI', 'CU']
    
    with open(receptor, 'r') as fp:
        pdb_raw = fp.read().splitlines()
        
    # get hetatm name for binding box setting
    hetatms = []
    prot_lines = []
    for line in pdb_raw:
        if line.startswith("HET   "):
            # HET 레코드 포맷: columns 8-10 (resname), 13 (chainID), 14-17 (seqNum)
            # 0-indexed: 7:10 (resname), 12 (chainID), 13:17 (seqNum)
            lig_name = line[7:10]
            lig_chain_id = line[12]
            lig_res_id = int(line[13:17].strip())
            if lig_name not in BUFFERS:
                hetatms.append((lig_name, lig_chain_id, lig_res_id))
        
        elif line.startswith("ATOM  ") or line.startswith("TER   "):
            prot_lines.append(line)
        
        elif line.startswith("HETATM"):
            res_name = line[17:20].strip()
            if res_name in METALS:
                prot_lines.append(line)
        
    prot_lines.append("END")
    prot_text = "\n".join(prot_lines)

    box_resname, box_chain, box_res_id = select_with_inquirer(hetatms)
    ligand_center, ligand_sdf = get_ligand_info(receptor, box_resname, box_chain, box_res_id)

    if ligand_center is None:
        print(f"Error: Ligand '{box_resname}' with chain '{box_chain}' and resID {box_res_id} not found.")
        return

    if ligand_out and ligand_sdf:
        with open(ligand_out, 'w') as f:
            f.write(ligand_sdf)
        print(f"Ligand conformation saved to {ligand_out}")

    if autobox_path:
        os.makedirs(autobox_path, exist_ok=True)
        config_path = os.path.join(autobox_path, f"{basename}_config.txt")

        # 'fixed' 방식 (기존 방식)
        if autobox_method == 'fixed':
            if ligand_center is None:
                print(f"Error: Ligand center could not be determined. Cannot generate config file.")
            else:
                with open(config_path, 'w') as fp:
                    fp.write(f'''center_x = {ligand_center[0]}
center_y = {ligand_center[1]}
center_z = {ligand_center[2]}

size_x = {autobox_size}
size_y = {autobox_size}
size_z = {autobox_size}''')
                print(f"Box config file was successfully generated at {config_path}")

        # 'optimal' 방식 (box_calculator 사용)
        elif autobox_method == 'optimal':
            if not ligand_sdf:
                print("Error: Could not generate SDF for ligand, cannot calculate optimal box size.")
            else:
                calculator = DockingBoxCalculator()
                with tempfile.NamedTemporaryFile(mode='w', suffix='.sdf', delete=False) as tmp_lig:
                    tmp_lig.write(ligand_sdf)
                    tmp_ligand_path = tmp_lig.name
                try:
                    if calculator.load_ligand_file(tmp_ligand_path):
                        try:
                            calculator.generate_autodock_config(
                                output_file=config_path,
                                use_optimal_size=True
                            )
                            print(f"Optimal box config file was successfully generated at {config_path}")
                        except ValueError as e:
                            print(f"Error generating optimal config: {e}")
                    else:
                        print(f"Error: Could not process ligand with BoxCalculator from {tmp_ligand_path}")
                finally:
                    os.remove(tmp_ligand_path)

    # prepare_receptor4.py의 경로를 찾기
    import shutil
    prepare_receptor4_path = shutil.which('prepare_receptor4.py')

    if prepare_receptor4_path is None:
        raise FileNotFoundError("prepare_receptor4.py를 찾을 수 없습니다. PATH에 AutoDockTools가 설치되어 있는지 확인하세요.")

    # 출력 디렉토리 생성
    os.makedirs(out_dir, exist_ok=True)
    output_pdbqt_path = os.path.join(out_dir, f"{basename}.pdbqt")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_rec:
        tmp_rec.write(prot_text)
        tmp_rec_path = tmp_rec.name
    try:
        command = f"python3 {prepare_receptor4_path} -r {tmp_rec_path} \
                    -A bonds_hydrogens \
                    -U nphs_lps_waters \
                    -o {output_pdbqt_path}"
        subprocess.run(command.split())
    finally:
        os.remove(tmp_rec_path)

    print(f"Successfully converted into pdbqt. Path: {output_pdbqt_path}")
            

if __name__ == '__main__':
    args = parse_arguments()
    protein_preparation(
        args.receptor,
        args.out_dir
        )