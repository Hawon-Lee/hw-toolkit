import os
import numpy as np
import pandas as pd
import argparse
import inquirer
import subprocess
import random
import multiprocessing
from datetime import datetime
from tqdm import tqdm
from functools import partial

def parse_arguments():
    
    parser = argparse.ArgumentParser(description='ligand preparation for docking.')
    
    parser.add_argument(
        '--receptor_pdbqt', '-r',
        type=str,
        required=True,
        help='전처리가 완료된 pdbqt 포맷 단백질 파일'
        )
    parser.add_argument(
        '--ligand_batch', '-l',
        type=str,
        required=True,
        help='ligand batch 파일을 읽어서 도킹을 수행합니다'
    )
    
    parser.add_argument(
        '--out_dir', '-o',
        type=str,
        required=True,
        help='도킹 리간드를 저장할 파일 경로'
    )
    parser.add_argument(
        '--autobox', '-a',
        type=str,
        required=True,
        help='단백질에 해당하는 box config 파일'
        )
    parser.add_argument(
        '--exhaustiveness', '-e',
        type=int,
        default=12,
        help='exhaustiveness'  
    )
    
    parser.add_argument(
        '--project_id', '-p',
        type=str,
        default=None,
        help="실험을 구분할 id 명입니다. 이 이름으로 도킹 결과가 저장됩니다."
    )
    
    parser.add_argument(
    '--num_processes', '-n',
    type=int,
    default=1,
    help='멀티프로세싱에 사용할 프로세스 수 (기본값: 1)'
    )
    
    args = parser.parse_args()
    return args

def process_ligand(ligand, receptor_pdbqt, autobox, exhaustiveness, project_dir):
    ligand_id = os.path.basename(ligand).split('.')[0]
    
    command = f"/home/tech/Hawon/DockingTool/autodock_vina_1_1_2_linux_x86/bin/vina \
        --receptor {receptor_pdbqt} \
        --ligand {ligand} \
        --config {autobox} \
        --exhaustiveness {exhaustiveness} \
        --out {project_dir}/{ligand_id}.pdbqt"
    subprocess.run(command.split())

    try:
        command = f"obabel -ipdbqt {project_dir}/{ligand_id}.pdbqt \
            -osdf -O {project_dir}/d{ligand_id}.sdf"
        subprocess.run(command.split())
    
        os.remove(f"{project_dir}/{ligand_id}.pdbqt")
        
        return ligand_id
    
    except Exception as e:
        print(f"Error processing {ligand}: {str(e)}")
        return ligand_id    
    

def vina_docking(args):
    with open(args.ligand_batch, 'r') as fp:
        ligand_batch = fp.read().splitlines()
    
    if args.project_id:
        project_dir = os.path.join(args.out_dir, args.project_id)
    else:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        project_dir = os.path.join(args.out_dir, timestamp)
    
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
    
    # 프로세스 수 설정
    num_processes = min(args.num_processes, multiprocessing.cpu_count())
    
    # process_ligand 함수의 일부 인자를 고정
    process_func = partial(
        process_ligand,
        receptor_pdbqt=args.receptor_pdbqt,
        autobox=args.autobox,
        exhaustiveness=args.exhaustiveness,
        project_dir=project_dir
    )
    
    # 멀티프로세싱 풀 생성 및 작업 실행
    with multiprocessing.Pool(processes=num_processes) as pool:
        # imap_unordered를 사용하여 작업이 완료되는 대로 결과 반환
        for _ in tqdm(pool.imap_unordered(process_func, ligand_batch), total=len(ligand_batch)):
            pass
        
if __name__ == '__main__':
    vina_docking(parse_arguments())