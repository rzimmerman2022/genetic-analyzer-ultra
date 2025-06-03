import numpy as np
from sklearn.neighbors import KNeighborsClassifier

# These files would need to be present in the same directory or a specified path
# For this example, we assume they are in a 'data' subdirectory or accessible path
# and that the main script handles loading them appropriately.
# If running this standalone, these paths would need to be adjusted.
try:
    REF = np.load("data/1000g_pca.npy")          # (N_ref, 10) PCA coordinates from 1000 Genomes
    LABELS = np.load("data/1000g_superpop.npy")  # Super-population labels (EU/AFR/EAS/SAS/AMR)
    KNN = KNeighborsClassifier(n_neighbors=7).fit(REF, LABELS)
except FileNotFoundError:
    print("WARNING: Ancestry reference data (1000g_pca.npy, 1000g_superpop.npy) not found. Ancestry inference will be disabled.")
    KNN = None
except Exception as e:
    print(f"WARNING: Error loading ancestry reference data: {e}. Ancestry inference will be disabled.")
    KNN = None

def infer_superpop(sample_pcs: np.ndarray) -> str:
    """
    Infers the super-population label (EU, AFR, EAS, SAS, AMR) for a sample
    given its principal components (PCs).

    Args:
        sample_pcs: A NumPy array of shape (n_features,) or (1, n_features)
                    representing the sample's PCs (e.g., first 10 PCs).

    Returns:
        A string representing the inferred super-population, or "UNKNOWN" if
        inference cannot be performed (e.g., KNN model not loaded).
    """
    if KNN is None:
        return "UNKNOWN"
    
    if sample_pcs.ndim == 1:
        sample_pcs_reshaped = sample_pcs.reshape(1, -1)
    elif sample_pcs.ndim == 2 and sample_pcs.shape[0] == 1:
        sample_pcs_reshaped = sample_pcs
    else:
        raise ValueError("sample_pcs must be a 1D array or a 2D array with one row.")

    try:
        prediction = KNN.predict(sample_pcs_reshaped)
        return str(prediction[0])
    except Exception as e:
        print(f"Error during ancestry prediction: {e}")
        return "UNKNOWN"
