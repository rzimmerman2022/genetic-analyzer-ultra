import json
import hashlib
from datetime import datetime

ANALYSIS_VERSION = '3.2.0' # Incremented minor version for QA hardening

DB_VERSIONS = {
    'ClinVar': '2025-03-01', # Placeholder, update with actual dates used
    'GWAS_Catalog': 'e110_r2025-02-15', # Placeholder
    'dbSNP': 'Build 156', # Placeholder
    'gnomAD': 'v4.0', # Placeholder
    'PharmGKB': '2025-03-10', # Placeholder
    'RefSeq': 'Release 220', # Placeholder
    # Add versions for any specific PRS models or other databases used
    'SSGAC_EA_PRS_Model': 'Okbay_et_al_2022_NatureGenetics',
    'CardiogramC4D_CAD_PRS_Model': 'Inouye_et_al_2018_JACC'
}

# To be called at the start of an analysis run
def get_initial_provenance():
    return {
        'analysis_script_version': ANALYSIS_VERSION,
        'database_versions_used': DB_VERSIONS,
        'analysis_start_time_utc': datetime.utcnow().isoformat()
    }

# To be called at the end of an analysis run
def finalize_provenance(provenance_dict: dict, results: dict) -> dict:
    """
    Adds a reproducibility hash and end time to the provenance dictionary.
    """
    # Ensure results are consistently ordered for hashing
    # Convert complex objects like numpy arrays or pandas Series/DataFrames to basic types if they exist in results
    # For simplicity, this example assumes results are already JSON serializable
    # In a real scenario, you might need a more robust way to serialize results
    try:
        serialized_results = json.dumps(results, sort_keys=True, default=str) # default=str for non-serializable
    except TypeError:
        # Fallback if json.dumps fails even with default=str (e.g. custom objects without __str__)
        # This is a very basic fallback and might not be perfectly reproducible if object order changes
        serialized_results = str(sorted(results.items()))

    provenance_dict['reproducibility_hash_sha256'] = hashlib.sha256(serialized_results.encode('utf-8')).hexdigest()
    provenance_dict['analysis_end_time_utc'] = datetime.utcnow().isoformat()
    return provenance_dict

def provenance_hash(results: dict) -> str:
    """
    Provided for direct use if only the hash is needed, matching the checklist.
    This is essentially what finalize_provenance does for the hash part.
    """
    try:
        serialized_results = json.dumps(results, sort_keys=True, default=str)
    except TypeError:
        serialized_results = str(sorted(results.items()))
    return hashlib.sha256(serialized_results.encode('utf-8')).hexdigest()
