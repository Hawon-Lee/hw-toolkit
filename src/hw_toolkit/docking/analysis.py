from pathlib import Path
import sys
from typing import List, Union


def parse_pdbqt(pdbqt: Union[str, Path]) -> List[float]:
    """
    Parse PDBQT file and extract VINA docking scores.
    
    Args:
        pdbqt: Path to the PDBQT file (can be string or Path object)
        
    Returns:
        List of docking scores as floats
    """
    pdbqt_path = Path(pdbqt)
    if not pdbqt_path.exists():
        print(f"PDBQT 파일이 존재하지 않습니다: {pdbqt}")
        return []

    text = pdbqt_path
    
    content = text.read_text()
    lines = content.splitlines()
    
    scores = []
    for i, line in enumerate(lines):
        if line.startswith("REMARK VINA RESULT"):
            score = float(line.split()[3])
            scores.append(score)

    return scores
    

def main():
    pdbqt = sys.argv[1]
    return parse_pdbqt(pdbqt)
    

if __name__ == "__main__":
    result = main()
    print(result)