import numpy as np
import re

class DockingBoxCalculator:
    """
    RDKit 없이 순수 파이썬으로 ligand 파일을 파싱하여 docking box 중심 좌표를 계산하는 클래스
    """
    
    def __init__(self):
        self.ligand_coords = None
        self.atom_symbols = None
        self.center = None
        self.box_size = None
        self.atomic_masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
            'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904,
            'B': 10.811, 'Si': 28.086, 'Se': 78.96, 'Mg': 24.305, 'Ca': 40.078,
            'Fe': 55.845, 'Zn': 65.38, 'Cu': 63.546, 'Mn': 54.938, 'Co': 58.933,
            'Ni': 58.693, 'Al': 26.982, 'K': 39.098, 'Na': 22.990
        }
    
    def parse_mol2_file(self, mol2_file):
        """
        MOL2 파일을 직접 파싱하여 원자 좌표를 추출
        
        Args:
            mol2_file (str): MOL2 파일 경로
            
        Returns:
            bool: 성공 여부
        """
        try:
            with open(mol2_file, 'r') as file:
                lines = file.readlines()
            
            coords = []
            atom_symbols = []
            in_atom_section = False
            
            for line in lines:
                line = line.strip()
                
                # @<TRIPOS>ATOM 섹션 시작
                if line == "@<TRIPOS>ATOM":
                    in_atom_section = True
                    continue
                
                # 다른 섹션 시작시 원자 섹션 종료
                if line.startswith("@<TRIPOS>") and line != "@<TRIPOS>ATOM":
                    in_atom_section = False
                    continue
                
                # 원자 정보 파싱
                if in_atom_section and line:
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            # MOL2 형식: atom_id atom_name x y z atom_type ...
                            x = float(parts[2])
                            y = float(parts[3])
                            z = float(parts[4])
                            atom_type = parts[5].split('.')[0]  # 원소 기호만 추출
                            
                            coords.append([x, y, z])
                            atom_symbols.append(atom_type)
                        except (ValueError, IndexError):
                            continue
            
            if not coords:
                raise ValueError("MOL2 파일에서 원자 좌표를 찾을 수 없습니다.")
            
            self.ligand_coords = np.array(coords)
            self.atom_symbols = atom_symbols
            print(f"MOL2에서 {len(coords)}개 원자의 좌표를 로드했습니다.")
            return True
            
        except Exception as e:
            print(f"MOL2 파일 파싱 오류: {e}")
            return False
    
    def parse_sdf_file(self, sdf_file):
        """
        SDF 파일을 직접 파싱하여 원자 좌표를 추출
        
        Args:
            sdf_file (str): SDF 파일 경로
            
        Returns:
            bool: 성공 여부
        """
        try:
            with open(sdf_file, 'r') as file:
                lines = file.readlines()
            
            # SDF 헤더 정보 찾기 (4번째 줄에 원자/결합 수)
            if len(lines) < 4:
                raise ValueError("SDF 파일 형식이 올바르지 않습니다.")
            
            counts_line = lines[3].strip()
            if len(counts_line) >= 6:
                try:
                    num_atoms = int(counts_line[:3])
                    num_bonds = int(counts_line[3:6])
                except ValueError:
                    raise ValueError("SDF 파일의 원자/결합 수를 읽을 수 없습니다.")
            else:
                raise ValueError("SDF 파일 형식이 올바르지 않습니다.")
            
            coords = []
            atom_symbols = []
            
            # 원자 정보 파싱 (5번째 줄부터)
            for i in range(4, 4 + num_atoms):
                if i >= len(lines):
                    break
                
                line = lines[i].strip()
                parts = line.split()
                
                if len(parts) >= 4:
                    try:
                        # SDF 형식: x y z atom_symbol ...
                        x = float(parts[0])
                        y = float(parts[1])
                        z = float(parts[2])
                        atom_symbol = parts[3]
                        
                        coords.append([x, y, z])
                        atom_symbols.append(atom_symbol)
                    except (ValueError, IndexError):
                        continue
            
            if not coords:
                raise ValueError("SDF 파일에서 원자 좌표를 찾을 수 없습니다.")
            
            self.ligand_coords = np.array(coords)
            self.atom_symbols = atom_symbols
            print(f"SDF에서 {len(coords)}개 원자의 좌표를 로드했습니다.")
            return True
            
        except Exception as e:
            print(f"SDF 파일 파싱 오류: {e}")
            return False
    
    def parse_pdb_file(self, pdb_file):
        """
        PDB 파일을 직접 파싱하여 원자 좌표를 추출 (HETATM 레코드)
        
        Args:
            pdb_file (str): PDB 파일 경로
            
        Returns:
            bool: 성공 여부
        """
        try:
            with open(pdb_file, 'r') as file:
                lines = file.readlines()
            
            coords = []
            atom_symbols = []
            
            for line in lines:
                # HETATM 또는 ATOM 레코드 찾기
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    try:
                        # PDB 형식: 고정 폭 필드
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        element = line[76:78].strip()
                        
                        # element가 없으면 atom name에서 추출
                        if not element:
                            atom_name = line[12:16].strip()
                            element = re.sub(r'[^A-Za-z]', '', atom_name)[:2]
                        
                        coords.append([x, y, z])
                        atom_symbols.append(element)
                    except (ValueError, IndexError):
                        continue
            
            if not coords:
                raise ValueError("PDB 파일에서 원자 좌표를 찾을 수 없습니다.")
            
            self.ligand_coords = np.array(coords)
            self.atom_symbols = atom_symbols
            print(f"PDB에서 {len(coords)}개 원자의 좌표를 로드했습니다.")
            return True
            
        except Exception as e:
            print(f"PDB 파일 파싱 오류: {e}")
            return False
    
    def load_ligand_file(self, file_path):
        """
        파일 확장자에 따라 적절한 파서를 선택하여 ligand 로드
        
        Args:
            file_path (str): ligand 파일 경로
            
        Returns:
            bool: 성공 여부
        """
        file_ext = file_path.lower().split('.')[-1]
        
        if file_ext == 'mol2':
            return self.parse_mol2_file(file_path)
        elif file_ext == 'sdf':
            return self.parse_sdf_file(file_path)
        elif file_ext == 'pdb':
            return self.parse_pdb_file(file_path)
        else:
            print(f"지원하지 않는 파일 형식입니다: {file_ext}")
            print("지원 형식: mol2, sdf, pdb")
            return False
    
    def calculate_geometric_center(self):
        """
        Ligand의 기하학적 중심 계산
        
        Returns:
            numpy.ndarray: [x, y, z] 중심 좌표
        """
        if self.ligand_coords is None:
            raise ValueError("먼저 ligand를 로드해주세요.")
        
        self.center = np.mean(self.ligand_coords, axis=0)
        return self.center
    
    def calculate_mass_center(self):
        """
        Ligand의 질량 중심 계산 (원자량 고려)
        
        Returns:
            numpy.ndarray: [x, y, z] 질량 중심 좌표
        """
        if self.ligand_coords is None or self.atom_symbols is None:
            raise ValueError("먼저 ligand를 로드해주세요.")
        
        total_mass = 0
        weighted_coords = np.zeros(3)
        
        for i, (coord, symbol) in enumerate(zip(self.ligand_coords, self.atom_symbols)):
            mass = self.atomic_masses.get(symbol, 12.011)  # 기본값: 탄소
            total_mass += mass
            weighted_coords += coord * mass
        
        self.center = weighted_coords / total_mass
        return self.center
    
    def calculate_bounding_box(self, padding=5.0):
        """
        Ligand 주변의 bounding box 크기 계산
        
        Args:
            padding (float): 각 방향으로 추가할 여백 (Angstrom)
            
        Returns:
            tuple: (box_size_x, box_size_y, box_size_z)
        """
        if self.ligand_coords is None:
            raise ValueError("먼저 ligand를 로드해주세요.")
        
        min_coords = np.min(self.ligand_coords, axis=0)
        max_coords = np.max(self.ligand_coords, axis=0)
        
        # 각 차원의 크기 + 패딩
        self.box_size = (max_coords - min_coords) + (2 * padding)
        
        return tuple(self.box_size)
    
    def calculate_radius_of_gyration(self):
        """
        Ligand의 회전반경 계산 (2.9배 규칙 적용용)
        
        Returns:
            float: 회전반경 값
        """
        if self.ligand_coords is None:
            raise ValueError("먼저 ligand를 로드해주세요.")
        
        if self.center is None:
            self.calculate_geometric_center()
        
        # 각 원자의 중심으로부터 거리의 제곱평균의 제곱근
        distances_squared = np.sum((self.ligand_coords - self.center)**2, axis=1)
        rog = np.sqrt(np.mean(distances_squared))
        
        return rog
    
    def calculate_optimal_box_size(self):
        """
        최적화된 box 크기 계산 (회전반경의 2.9배)
        
        Returns:
            float: 최적 box 크기
        """
        rog = self.calculate_radius_of_gyration()
        return rog * 2.9
    
    def generate_autodock_config(self, output_file="docking_config.txt",
                               use_optimal_size=False):
        """
        AutoDock Vina 설정 파일 생성
        
        Args:
            output_file (str): 출력 파일명
            use_optimal_size (bool): 최적화된 box 크기 사용 여부
        """
        if self.center is None:
            self.calculate_mass_center()  # 질량 중심 사용 (더 정확)
        
        if self.center is None:
            raise ValueError("Center could not be calculated. Cannot generate config.")
        
        if use_optimal_size:
            optimal_size = self.calculate_optimal_box_size()
            box_x = box_y = box_z = optimal_size
        else:
            if self.box_size is None:
                self.calculate_bounding_box()
            
            if self.box_size is None:
                raise ValueError("Box size could not be calculated. Cannot generate config.")
            box_x, box_y, box_z = self.box_size
        
        config_content = f"""# AutoDock Vina 설정 파일

center_x = {self.center[0]:.3f}
center_y = {self.center[1]:.3f}
center_z = {self.center[2]:.3f}

size_x = {box_x:.1f}
size_y = {box_y:.1f}
size_z = {box_z:.1f}
"""
        
        with open(output_file, 'w') as f:
            f.write(config_content)
        
        print(f"AutoDock Vina 설정 파일이 생성되었습니다: {output_file}")
    
    def generate_glide_input(self, output_file="glide_grid.inp"):
        """
        Glide용 그리드 입력 파일 생성
        
        Args:
            output_file (str): 출력 파일명
        """
        if self.center is None:
            self.calculate_mass_center()
        
        if self.box_size is None:
            self.calculate_bounding_box()
        
        if self.center is None or self.box_size is None:
            raise ValueError("Center or Box Size could not be calculated. Cannot generate Glide input.")
        
        glide_content = f"""# Glide 그리드 생성 입력 파일
GRIDFILE   grid.zip
RECEP_FILE receptor.mae

# 그리드 중심 좌표
GRID_CENTER {self.center[0]:.3f}, {self.center[1]:.3f}, {self.center[2]:.3f}

# 내부/외부 박스 크기 (Å)
INNERBOX {self.box_size[0]:.1f}, {self.box_size[1]:.1f}, {self.box_size[2]:.1f}
OUTERBOX {self.box_size[0]+10:.1f}, {self.box_size[1]+10:.1f}, {self.box_size[2]+10:.1f}
"""
        
        with open(output_file, 'w') as f:
            f.write(glide_content)
        
        print(f"Glide 그리드 입력 파일이 생성되었습니다: {output_file}")
    
    def print_summary(self):
        """
        계산 결과 요약 출력
        """
        if self.ligand_coords is None:
            print("Ligand가 로드되지 않았습니다.")
            return
        
        print("\n" + "="*50)
        print("          DOCKING BOX 설정 요약")
        print("="*50)
        
        if self.center is not None:
            print(f"중심 좌표: ({self.center[0]:.3f}, {self.center[1]:.3f}, {self.center[2]:.3f})")
        
        if self.box_size is not None:
            print(f"Box 크기:  {self.box_size[0]:.1f} × {self.box_size[1]:.1f} × {self.box_size[2]:.1f} Å")
        
        rog = self.calculate_radius_of_gyration()
        optimal_size = rog * 2.9
        print(f"회전반경:  {rog:.3f} Å")
        print(f"최적 Box: {optimal_size:.1f} Å (2.9 × RoG)")
        
        print(f"원자 수:   {len(self.ligand_coords)}")
        print("="*50 + "\n")


def main():
    """
    사용 예시 및 데모
    """
    print("=== RDKit 없는 Docking Box Calculator ===\n")
    
    # 실제 사용법 예시
    print("실제 사용법:")
    print("""
# 1. 인스턴스 생성
calculator = DockingBoxCalculator()

# 2. Ligand 파일 로드 (자동 형식 인식)
success = calculator.load_ligand_file("ligand.sdf")
# 또는
success = calculator.load_ligand_file("ligand.mol2")
# 또는
success = calculator.load_ligand_file("ligand.pdb")

# 3. 중심 좌표 계산
center = calculator.calculate_mass_center()  # 권장 (질량 중심)
# 또는
center = calculator.calculate_geometric_center()  # 기하학적 중심

# 4. Box 크기 계산
box_size = calculator.calculate_bounding_box(padding=5.0)

# 5. 설정 파일 생성
calculator.generate_autodock_config("vina_config.txt")

# 6. 결과 확인
calculator.print_summary()
    """)
    
    # 데모용 가상 데이터
    print("\n=== 데모 실행 (가상 벤젠 분자) ===")
    calculator = DockingBoxCalculator()
    
    # 벤젠 고리의 좌표와 원소 정보
    demo_coords = np.array([
        [0.000, 1.400, 0.000],   # C
        [1.212, 0.700, 0.000],   # C
        [1.212, -0.700, 0.000],  # C
        [0.000, -1.400, 0.000],  # C
        [-1.212, -0.700, 0.000], # C
        [-1.212, 0.700, 0.000],  # C
        [0.000, 2.480, 0.000],   # H
        [2.148, 1.240, 0.000],   # H
        [2.148, -1.240, 0.000],  # H
        [0.000, -2.480, 0.000],  # H
        [-2.148, -1.240, 0.000], # H
        [-2.148, 1.240, 0.000]   # H
    ])
    
    demo_symbols = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    
    calculator.ligand_coords = demo_coords
    calculator.atom_symbols = demo_symbols
    
    # 중심 계산
    geo_center = calculator.calculate_geometric_center()
    mass_center = calculator.calculate_mass_center()
    
    print(f"기하학적 중심: ({geo_center[0]:.3f}, {geo_center[1]:.3f}, {geo_center[2]:.3f})")
    print(f"질량 중심:     ({mass_center[0]:.3f}, {mass_center[1]:.3f}, {mass_center[2]:.3f})")
    
    # Box 크기 계산
    box_size = calculator.calculate_bounding_box(padding=5.0)
    
    # 설정 파일 생성
    calculator.generate_autodock_config("demo_vina.txt")
    calculator.generate_autodock_config("demo_vs.txt", use_optimal_size=True)
    calculator.generate_glide_input("demo_glide.inp")
    
    # 요약 출력
    calculator.print_summary()


if __name__ == "__main__":
    main()