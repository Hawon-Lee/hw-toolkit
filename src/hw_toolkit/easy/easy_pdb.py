# python version >= 3.8 

from rdkit import Chem
import numpy as np
import os
import json
import requests
import pickle
import time, random
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from collections import Counter

def open_file(path, is_pickle=False):
    if is_pickle:
        with open(path, 'rb') as fp:
            return pickle.load(fp)
    with open(path, 'r') as fp:
        return fp.read()
    
def write_file(file, path):
    with open(path, 'w') as fp:
        fp.write(file)

def query_pdb_normal(query):
    # HTTP request headers
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    req = json.dumps(query)
    response = requests.post(url, req, headers=headers)
    
    # check response
    if response.status_code == 200:
        return response.json()
    else:
        return f"Response Error: {response.status_code}, {response.text}"

import requests
import json

def search_api_pdb_graphql(query, variables=None):
    """
    Function to send queries to the RCSB PDB API and return the results.
    Supports both GraphQL queries and PDB JSON query format.
    
    Parameters:
    -----------
    query : str or dict
        GraphQL query string or PDB JSON query object
    variables : dict, optional
        Dictionary containing variables to be used in the GraphQL query (default: None)
        
    Returns:
    --------
    dict
        Dictionary containing the API response
    """
    # Check if the input is a JSON string or dict (PDB search query format)
    is_json_query = False
    
    # If it's a dictionary, it's a JSON query
    if isinstance(query, dict):
        is_json_query = True
        query_data = query
    # If it's a string, try to detect format
    elif isinstance(query, str):
        query = query.strip()
        # Try to parse it as JSON
        try:
            query_data = json.loads(query)
            # If parse successful and has typical PDB query structure, treat as JSON query
            if isinstance(query_data, dict) and "query" in query_data:
                is_json_query = True
        except json.JSONDecodeError:
            # Not valid JSON, probably a GraphQL query
            is_json_query = False
            
    if is_json_query:
        # PDB JSON query format - use search API endpoint
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        # Prepare request data - already parsed to dict in the detection step
        data = query_data

    else:
        # GraphQL query format - use GraphQL API endpoint
        url = "https://data.rcsb.org/graphql"
        
        # Prepare request data
        data = {
            "query": query
        }
        
        # Add variables if provided
        if variables:
            data["variables"] = variables
    
    # HTTP request headers
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    try:
        # Send API request
        response = requests.post(url, headers=headers, data=json.dumps(data))
        response.raise_for_status()  # Raise exception for HTTP errors
        
        # Parse JSON response
        result = response.json()
        
        # Check for and print errors in GraphQL response
        if not is_json_query and "errors" in result:
            print("GraphQL query errors:", result["errors"])
        
        return result
    
    except requests.exceptions.RequestException as e:
        print(f"API request error: {e}")
        return {"error": str(e)}
    except json.JSONDecodeError:
        print("JSON parsing error: Invalid response format")
        return {"error": "JSON parsing error"}
    except Exception as e:
        print(f"Error occurred: {e}")
        return {"error": str(e)}

def data_api_pdb_graphql(query, pdb_id):
    # GraphQL endpoint URL
    url = "https://data.rcsb.org/graphql"

    # Set variables
    variables = {
        "id": pdb_id
    }

    # HTTP request headers
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    # send POST request
    response = requests.post(url, json={"query": query, "variables": variables}, headers=headers)

    # check response
    if response.status_code == 200:
        return response.json()
    else:
        return f"Error: {response.status_code}, {response.text}"
    
    

    
def pdb_search_from_uniprot_id(uniprot_id, mutation_ok=False, max_resolution=2.7):

    if mutation_ok:
        mut_operator = 'greater'
    else:
        mut_operator = 'equals'
        
    query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "in",
                        "negation": False,
                        "value": [
                            uniprot_id
                        ]
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                        "operator": "exact_match",
                        "value": "UniProt",
                        "negation": False
                        }
                    }
                    ],
                    "label": "nested-attribute"
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "attribute": "entity_poly.rcsb_mutation_count",
                    "operator": mut_operator,
                    "negation": False,
                    "value": 0
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less",
                    "negation": False,
                    "value": max_resolution
                    }
                }
                ],
                "label": "text"
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                "start": None,
                "rows": None
                },
                "results_content_type": [
                "experimental"
                ],
                "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
                ],
                "scoring_strategy": "combined"
            }
            }

    result_set = []
    
    # for pagenation
    start, rows = 0, 25
    
    # query pagenation option update
    query['request_options']['paginate']['start'], query['request_options']['paginate']['rows'] = start, rows
    
    # get response with treating pagination
    response = query_pdb_normal(query)
    
    if type(response) == str:
        return {}
    
    total_count = response['total_count']
    result_set += response['result_set']
    
    while start + rows < total_count:
        start += rows
        # query pagenation update again to get full data
        query['request_options']['paginate']['start'], query['request_options']['paginate']['rows'] = start, rows
        response = query_pdb_normal(query)
        result_set += response['result_set']
        
    response['result_set'] = result_set
    return response

def download_pdb_files(pdb_ids: list[str], output_dir: str = "pdb_files", max_workers:int = 5):
    """
    Download PDB files for given PDB IDs with rate limiting.

    Args:
        pdb_ids (List[str]): List of PDB IDs to download.
        output_dir (str): Directory to save downloaded files. Defaults to "pdb_files".
        max_workers (int): Maximum number of concurrent downloads. Defaults to 5.

    Returns:
        Dict[str, Optional[str]]: Dictionary mapping PDB IDs to their file paths, or None if download failed.
    """
    
    def download_pdb(pdb_id, output_dir):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
        
        # Add rate limiting
        time.sleep(random.uniform(0.5, 1.5))
        
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            return pdb_id, output_path
        else:
            return pdb_id, None
    
    os.makedirs(output_dir, exist_ok=True)
    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(download_pdb, pdb_id, output_dir) for pdb_id in pdb_ids]
        
        for future in tqdm(as_completed(futures), total=len(pdb_ids), desc="Downloading PDB files"):
            pdb_id, file_path = future.result()
            results[pdb_id] = file_path
            
    print('Downloading PDB files has done.')
    return results

def get_parsed_chains(pdb_file_path: str) -> dict[str, str]:
    """
    Parse a PDB file to extract chain and ligand information.

    param:
        pdb_file_path: Path to the PDB file
    return:
        Dictionary containing parsed information
    """
    result_parsed_dict = {}
    refined_text = []
    ter_indice = [-1]
    groups = []

    # Read PDB file and extract relevant lines
    with open(pdb_file_path, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM', 'TER')):
                refined_text.append(line)

    # Find indices of TER lines
    for i, r_line in enumerate(refined_text):
        if r_line.startswith('TER'):
            ter_indice.append(i)
    ter_indice.append(len(refined_text))

    # Create groups between TER lines
    for i in range(len(ter_indice) - 1):
        start = ter_indice[i] + 1  # Exclude TER line
        end = ter_indice[i+1]
        group = refined_text[start:end]
        groups.append(group)

    # Process each group
    for group in groups:
        key_identifier = ""
        headers = []
        residue_names = []
        chains = []

        # Collect header, residue_name, and chain for each line
        for line in group:
            header, residue_name, chain = line[0:6], line[17:20], line[21]
            headers.append(header)
            residue_names.append(residue_name)
            chains.append(chain)

        # Determine if it's a protein or non-protein
        header_count = Counter(headers)
        if len(header_count) == 1:
            add_1 = 'protein' if headers[0] == 'ATOM  ' else 'non_protein'
        elif len(header_count) == 2:
            add_1 = 'protein' if list(header_count.keys())[0] == 'ATOM  ' else 'non_protein'
        else:
            raise ValueError('HEADER column has another string other than ATOM or HETATM')
        key_identifier += add_1

        # Determine chain identifier
        add_2 = '_' + chains[0] if len(set(chains)) == 1 else '_X'
        key_identifier += add_2

        # Process non-protein groups
        if 'non_protein' in key_identifier:
            if len(set(residue_names)) == 1:
                add_3 = '_' + residue_names[0]
                key_identifier += add_3
            elif 'HOH' in set(residue_names):
                for ligand_name in set(residue_names):
                    ligand_index = [i for i, j in enumerate(residue_names) if j == ligand_name]
                    start = ligand_index[0]
                    end = ligand_index[-1] + 1
                    ligand_group = group[start:end]
                    result_parsed_dict[key_identifier + '_' + ligand_name] = "".join(ligand_group)
                continue
            else:
                add_3 = '_UNK'
                key_identifier += add_3

        result_parsed_dict[key_identifier] = "".join(group)

    return result_parsed_dict

# def get_filtered_chains(parsed_chains: dict) -> dict:
#     """
#     Filter the parsed PDB chains to retain only one standard chain and one HETATM chain which represents ligand.

#     This function takes the output of func 'parse_pdb_chains' and returns a new dictionary
#     with at most two entries: one standard chain (the first one encountered) and
#     one HETATM chain (the longest one excepts HOH).

#     Args:
#         parsed_chains: A dictionary as returned by parse_pdb_chains.

#     Returns:
#         A new dictionary with at most two entries: one standard chain and one HETATM chain.
#     """
#     filtered_chains = {}
#     standard_chain = None
#     hetatm_chain = None
#     max_hetatm_length = 0

#     for key, value in parsed_chains.items():
#         if '_' not in key and standard_chain is None:
#             # This is the first standard chain encountered
#             standard_chain = (key, value)
#         elif ('_' in key) and ('HOH' not in key):
#             # This is a HETATM chain
#             if len(value) > max_hetatm_length:
#                 hetatm_chain = (key, value)
#                 max_hetatm_length = len(value)

#     if standard_chain:
#         filtered_chains[standard_chain[0]] = standard_chain[1]
#     if hetatm_chain:
#         filtered_chains[hetatm_chain[0]] = hetatm_chain[1]

#     return filtered_chains

def extract_pocket(protein_pdb:str, ligand_pdb:str, distance_cutoff=5) -> str:
    '''
    Extract pocket atoms from protein_pdb.
    
    The pocket is composed of protein atom pdb lines whose distance from the ligand atom is within a cutoff.
    '''
    p_total_coords = []
    l_total_coords = []

    for p_line in protein_pdb.splitlines()[:-1]:
        p_coords = p_line[30:38], p_line[38:46], p_line[46:54]
        p_coords = [float(v) for v in p_coords]
        p_total_coords.append(p_coords)
        
    for l_line in ligand_pdb.splitlines()[:-1]:
        l_coords = l_line[30:38], l_line[38:46], l_line[46:54]
        l_coords = [float(v) for v in l_coords]
        l_total_coords.append(l_coords)

    p_total_coords, l_total_coords = np.array(p_total_coords), np.array(l_total_coords)
    distance_mat = np.sqrt(np.square(p_total_coords[:, np.newaxis, :] - l_total_coords[np.newaxis, :, :]).sum(axis=2))
    p_index = np.nonzero(distance_mat < distance_cutoff)[0]
    p_index = sorted(list(set(p_index)))

    pocket_lines = [protein_pdb.splitlines()[i] for i in p_index] + ['TER']
    
    return '\n'.join(pocket_lines)

# # use example
# if __name__ == '__main__':
#     pdb_file_path = './8V2F.pdb'
#     parsed_chains = get_parsed_chains(pdb_file_path=pdb_file_path)
#     filtered_chains = get_filter_chains(parsed_chains=parsed_chains)
#     p, l = list(filtered_chains.items())[0], list(filtered_chains.items())[1]
#     extract_pocket(p, l, distance_cutoff=5)