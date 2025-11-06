import subprocess
import os
from Bio import SeqIO

def run_blastp(query_file, db_path, output_file):
    """
    BLASTP를 실행하고 결과를 txt 형식으로 저장합니다.
    """
    cmd = f"blastp -query {query_file} -db {db_path} -out {output_file} -evalue 1e-5 -num_threads 4"
    subprocess.run(cmd, shell=True, check=True)

def parse_blast_result(blast_output):
    """
    BLAST XML 결과를 파싱하고 최상위 히트의 UniProt ID를 반환합니다.
    """
    with open(blast_output, 'r') as result_handle:
        blast_records = result_handle.readlines()
        for record in blast_records:
            if record.startswith('sp|'):
                # 첫 번째 정렬에서 UniProt ID 추출
                uniprot_id = record.split('|')[1]
                return uniprot_id
    return None

def main(query_sequence, db_path):
    # 임시 파일 생성
    temp_query_file = "temp_query.fasta"
    temp_output_file = "temp_output.txt"

    try:
        # 쿼리 시퀀스를 임시 FASTA 파일로 저장
        with open(temp_query_file, "w") as f:
            f.write(">query\n" + query_sequence)

        # BLASTP 실행
        run_blastp(temp_query_file, db_path, temp_output_file)

        # 결과 파싱
        uniprot_id = parse_blast_result(temp_output_file)

    finally:
        # 임시 파일 삭제
        for file in (temp_query_file, temp_output_file):
            if os.path.exists(file):
                os.remove(file)
                
    return uniprot_id

if __name__ == "__main__":
    # 예시 사용
    query_sequence = 'MAESAGASSFFPLVVLLLAGSGGSGPRGVQALLCACTSCLQANYTCETDGACMVSIFNLDGMEHHVRTCIPKVELVPAGKPFYCLSSEDLRNTHCCYTDYCNRIDLRVPSGHLKEPEHPSMWGPVELVGIIAGPVFLLFLIIIIVFLVINYHQRVYHNRQRLDMEDPSCEMCLSKDKTLQDLVYDLSTSGSGSGLPLFVQRTVARTIVLQEIIGKGRFGEVWRGRWRGGDVAVKIFSSREERSWFREAEIYQTVMLRHENILGFIAADNKDNGTWTQLWLVSDYHEHGSLFDYLNRYTVTIEGMIKLALSAASGLAHLHMEIVGTQGKPGIAHRDLKSKNILVKKNGMCAIADLGLAVRHDAVTDTIDIAPNQRVGTKRYMAPEVLDETINMKHFDSFKCADIYALGLVYWEIARRCNSGGVHEEYQLPYYDLVPSDPSIEEMRKVVCDQKLRPNIPNWWQSYEALRVMGKMMRECWYANGAARLTALRIKKTLSQLSVQEDVKI'
    db_path = "/nas/binding_affinity/uniprot_sprot_db"  # UniProt 데이터베이스 경로를 지정하세요
    uniprot_id = main(query_sequence, db_path)
    if uniprot_id:
        print(f"찾은 UniProt ID: {uniprot_id}")
    else:
        print("매치되는 UniProt ID를 찾지 못했습니다.")