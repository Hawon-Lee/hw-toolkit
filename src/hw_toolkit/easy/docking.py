# requirements : openbabel
import py3Dmol
import subprocess
import tempfile
import os
import datetime
from rdkit import Chem

def neutralize_smiles(smiles: str) -> str:
    '''
    Neutralize all charges from SMILES before execute protonation.
    '''
    command = f'obabel -:{smiles} -osmi --neutralize'
    command_result = subprocess.run(command.split(), capture_output=True, text=True).stdout
    smiles_n = command_result.rstrip()
    
    if not smiles_n:
        raise ValueError(f'smiles cannot be neutralized: {smiles}')
    
    return smiles_n

def protonate_smiles(smiles_n: str, pH=7.4) -> str:
    '''
    Protonates SMILES before generate conformer.
    '''
    
    command = f'obabel -:{smiles_n} -osmi -p {pH}'
    command_result = subprocess.run(command.split(), capture_output=True, text=True).stdout
    smiles_h = command_result.rstrip()
    
    if not smiles_h:
        raise ValueError(f'smiles cannot be protonated: {smiles_n}')
    
    return smiles_h

def generate_conformer_from_smiles(smiles_h: str, extension: str = 'sdf', timeout: int = 15) -> str:
    '''
    generate 3D conformer from protonated smiles, and return format is string whose type is choosed extension.
    
    Args:
        smiles_h: Any smiles to be 3D conformer. appropriately protonated smiles are recommended.
        extension: File extension to convert smiles which in sdf, mol2, pdb, pdbqt
        timeout: Converting time limit (seconds)
        
    Returns:
        Strings as appropriate formats
    '''
    if extension not in ['sdf', 'mol2', 'pdb', 'pdbqt']:
        raise ValueError(f'Input extension must be sdf, mol2, pdb or pdbqt.')
    command = f'obabel -:{smiles_h} --gen3D -o{extension}'
    command_result = subprocess.run(command.split(), capture_output=True, timeout=timeout, text=True).stdout
    conformer = command_result.rstrip()

    if not conformer:
        raise ValueError(f'This smiles cannot be converted to 3d conformer: {smiles_h}')    
    
    return conformer

def protonate_proteins(protein_pdb: str, pH: float = 7.4, as_pdbqt: bool = False) -> str:
    """
    Protonate a protein structure at a given pH using OpenBabel.

    Args:
    protein_pdb (str): The protein structure in PDB format as a string.
    pH (float): The pH at which to protonate the protein. Default is 7.4.
    as_pdbqt (bool): If True, output as PDBQT format. If False, output as PDB. Default is False.

    Returns:
    str: The protonated protein structure as a string in PDB or PDBQT format.

    Raises:
    subprocess.CalledProcessError: If the OpenBabel command fails.
    """
    # Determine output format
    extension = 'pdbqt' if as_pdbqt else 'pdb'

    # Construct OpenBabel command
    command = f'obabel -ipdb -o{extension} -p {pH}'

    try:
        # Run OpenBabel command
        result = subprocess.run(
            command.split(),
            input=protein_pdb,
            capture_output=True,
            text=True,
            check=True  # This will raise CalledProcessError if the command fails
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        # Handle OpenBabel execution errors
        print(f"Error running OpenBabel: {e}")
        print(f"Error output: {e.stderr}")
        raise  # Re-raise the exception after logging
    
def obabel_conversion(file_text: str, in_format: str, out_format: str):
    command = f'obabel -i{in_format} -o{out_format}'
    try:
        result = subprocess.run(command.split(),
                    input = file_text.encode(),
                    capture_output=True,
                    text=False)
        return result.stdout.rstrip().decode()
    except subprocess.CalledProcessError as e:
        print(f"Error running Openbabel: {e}")
        print(f"Error output: {e.stderr}")
        raise

def run_docking_smina(
        protein: str,
        ligand: str,
        autobox_ligand: str,
        autobox_add: int = 4,
        protein_extension: str = 'pdbqt',
        ligand_extension: str = 'sdf',
        exhaustiveness: int = 12,
        num_modes: int = 9,
        logfile_path: str = './smina_log.log',
        job_id: int = 0) -> str:
    '''
    Perform molecular docking using Smina.

    This function takes protein and ligand structures as input, performs molecular docking
    using the Smina docking software, and returns the docked ligand structure. It also logs
    the docking process and any errors that occur.

    Arguments:
    -----------
    protein : str
        The protein structure in PDBQT format.
    ligand : str
        The ligand structure in SDF format.
    autobox_ligand : str
        The ligand structure used for autoboxing in SDF format.
    autobox_add : int, optional
        The amount to add to the autobox size (default is 4).
    protein_extension : str, optional
        The file extension for the protein file (default is 'pdbqt').
    ligand_extension : str, optional
        The file extension for the ligand file (default is 'sdf').
    exhaustiveness : int, optional
        The exhaustiveness of the search (default is 12).
    num_modes : int, optional
        The number of binding modes to generate (default is 9).
    logfile_path : str, optional
        The path to the log file (default is './smina_log.log').
    job_id : int, optional
        An identifier for the docking job and be reported at the logfile.
        Will not be reported if job_id is zero. (default is 0).

    Returns:
    --------
    str or None
        The docked ligand structure in SDF format if docking is successful, None otherwise.

    Notes:
    ------
    This function requires the Smina docking software to be installed and accessible
    in the system path. It uses temporary files for the docking process and logs the
    results and any errors to the specified log file.
    '''
    
    # Generate timestamp and separator for logging
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_separator = f"\n{'='*50}\n"
    job_id_str = f'job_id: {job_id}' if job_id else ''

    # Create and use a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Set paths for temporary files
        temp_protein = os.path.join(temp_dir, f"protein.{protein_extension}")
        temp_ligand = os.path.join(temp_dir, f"ligand.{ligand_extension}")
        temp_autobox_ligand = os.path.join(temp_dir, f"autobox_ligand.{ligand_extension}")
        temp_output = os.path.join(temp_dir, f"docked.{ligand_extension}")
        
        # Write input structures to temporary files
        with open(temp_protein, 'w') as fp:
            fp.write(protein)
        with open(temp_ligand, 'w') as fp:
            fp.write(ligand)
        with open(temp_autobox_ligand, 'w') as fp:
            fp.write(autobox_ligand)
        
        # Construct Smina docking command
        command = f'smina -r {temp_protein} \
                -l {temp_ligand} \
                -o {temp_output} \
                --autobox_ligand {temp_autobox_ligand} \
                --autobox_add {autobox_add} \
                --exhaustiveness {exhaustiveness} \
                --num_modes {num_modes} \
                --quiet'
                
        try:
            # Execute Smina docking
            result = subprocess.run(command.split(), check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Log successful docking results            
            with open(logfile_path, 'a') as fp:
                fp.write(f"{log_separator}DOCKING START: {timestamp} {job_id_str}\n")
                fp.write(result.stdout.rstrip()+'\n')
                fp.write(f"DOCKING END")
                
        except subprocess.CalledProcessError as e:
            # Handle and log docking errors
            print(f'Error : {e}, {job_id_str}')
            
            with open(logfile_path, 'a') as fp:
                fp.write(f"{log_separator}DOCKING START: {timestamp} {job_id_str}\n")
                fp.write(f"ERROR: DOCKING FAILED, {job_id_str}\n")
                fp.write(e.stderr)
            return None
        
        # Read and return docked ligand structure
        with open(temp_output, 'r') as fp:
            docked_sdf = fp.read()
            # print(docked_sdf)
            
        return docked_sdf
    
def visualize(protein_path, original_ligand_path, docked_ligand_path):
    '''
    May not works well. This needs to be update
    '''
    
    # Open empty pallete
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline',
                    'color':'black',
                    'width':0.1})

    # Add the protein structure
    view.addModel(open(protein_path, 'r').read(), format='pdb')
    Prot = view.getModel()
    Prot.setStyle({'cartoon':{'arrows':True,
                            'tubes':False,
                            'style':'oval',
                            'color':'white'}})
    # view.addSurface(py3Dmol.VDW, {'opacity':0.4, 'color':'white'})

    # Add the original ligand conformation
    view.addModel(open(original_ligand_path, 'r').read(), format='sdf')
    ref_m = view.getModel()
    ref_m.setStyle({}, {'stick':{'colorscheme':'redCarbon','radius':0.2}})

    # Add docked conformation
    results = Chem.SDMolSupplier(docked_ligand_path)
    p = Chem.MolToMolBlock(results[0], False)
    print('Reference : Red | Vina Pose : Cyan')
    print('Score: {}'.format(results[0].GetProp('minimizedAffinity')))

    view.addModel(p, 'mol')
    docked = view.getModel()
    docked.setStyle({}, {'stick':{'colorscheme':'cyanCarbon','radius':0.2}})

    view.zoomTo()
    view.show()
