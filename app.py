import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import os
import json

from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive
from google.oauth2.service_account import Credentials

# --- Page Configuration ---
st.set_page_config(
    page_title="UMAP Plotter from Google Drive",
    page_icon="🧬",
    layout="wide"
)

# --- Authentication (Service Account) ---
@st.cache_resource
def get_gdrive_service():
    """Authenticates with Google Drive using Service Account credentials."""
    # Load credentials from Streamlit secrets
    creds_json = st.secrets.google_credentials.service_account_json
    
    # The Pydrive2 library expects a file, so we write the JSON string to a temporary file
    with open("service_creds.json", "w") as f:
        f.write(creds_json)
    
    gauth = GoogleAuth()
    scope = ["https://www.googleapis.com/auth/drive.readonly"]
    gauth.credentials = Credentials.from_service_account_file("service_creds.json", scopes=scope)
    drive = GoogleDrive(gauth)
    
    # Clean up the temporary credentials file
    os.remove("service_creds.json")
    return drive

# --- Data Fetching Functions ---
@st.cache_data
def list_h5ad_files_in_folder(_drive, folder_id):
    """Lists .h5ad files in a specific Google Drive folder."""
    query = f"'{folder_id}' in parents and trashed=false"
    file_list = _drive.ListFile({'q': query}).GetList()
    return {file['title']: file['id'] for file in file_list if file['title'].endswith('.h5ad')}

def download_h5ad_from_drive(_drive, file_id, file_name):
    """Downloads an h5ad file from drive to a temporary local path."""
    h5ad_file = _drive.CreateFile({'id': file_id})
    h5ad_file.GetContentFile(file_name)
    return file_name

# --- App UI ---
st.title("🔬 Interactive UMAP from Curated Single-Cell Data")
st.write("Select a dataset from the dropdown below to generate and explore its UMAP visualization.")

try:
    drive = get_gdrive_service()
    folder_id = st.secrets.folder_id
    
    # List available files
    h5ad_files = list_h5ad_files_in_folder(drive, folder_id)

    if not h5ad_files:
        st.warning("No `.h5ad` files found. Please check your Folder ID and sharing permissions with the service account.")
        st.stop()
    
    # --- Main Panel for Plotting ---
    st.header("📊 UMAP Visualization")

    selected_file_title = st.selectbox("Select a dataset:", options=list(h5ad_files.keys()))

    if selected_file_title:
        # Load data if it's not already in session state or if a new file is selected
        if 'loaded_file' not in st.session_state or st.session_state.loaded_file != selected_file_title:
            with st.spinner(f"Loading '{selected_file_title}'..."):
                file_id = h5ad_files[selected_file_title]
                temp_file_path = download_h5ad_from_drive(drive, file_id, selected_file_title)
                
                try:
                    st.session_state['adata'] = sc.read_h5ad(temp_file_path)
                    st.session_state['loaded_file'] = selected_file_title
                    st.success(f"Successfully loaded {st.session_state.adata.n_obs} cells.")
                finally:
                    # Clean up the downloaded file
                    if os.path.exists(temp_file_path):
                        os.remove(temp_file_path)
    
    # --- Plotting Section ---
    if 'adata' in st.session_state:
        adata = st.session_state.adata
        
        categorical_obs = adata.obs.select_dtypes(include=['category', 'object']).columns.tolist()

        if not categorical_obs:
            st.warning("No categorical data found in `adata.obs` to color the plot.")
        else:
            color_by = st.selectbox(
                "Color UMAP plot by:",
                options=categorical_obs,
                help="Select a column from `adata.obs` to color the cells."
            )
            
            st.write(f"### UMAP for `{st.session_state.loaded_file}` colored by `{color_by}`")
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.umap(adata, color=color_by, ax=ax, show=False, legend_loc='on data')
            st.pyplot(fig)

except Exception as e:
    st.error(f"An error occurred: {e}")
    st.error("Please ensure your secrets are configured correctly and the Google Drive folder is shared with the service account's email.")
