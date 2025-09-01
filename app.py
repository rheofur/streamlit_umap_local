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
    page_title="UMAP & Gene Expression Plotter",
    page_icon="🧬",
    layout="wide"
)

# --- Authentication (Service Account) ---
@st.cache_resource
def get_gdrive_services():
    """Authenticates with Google Drive using Service Account credentials."""
    creds_json = st.secrets.google_credentials.service_account_json
    with open("service_creds.json", "w") as f:
        f.write(creds_json)
    
    gauth = GoogleAuth()
    scope = ["https://www.googleapis.com/auth/drive.readonly"]
    gauth.credentials = Credentials.from_service_account_file("service_creds.json", scopes=scope)
    drive = GoogleDrive(gauth)
    
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
st.title("🔬 Interactive UMAP & Gene Expression Plotter")
st.write("Select a dataset to visualize cell metadata or individual gene expression.")

try:
    drive = get_gdrive_services()
    folder_id = st.secrets.folder_id
    h5ad_files = list_h5ad_files_in_folder(drive, folder_id)

    if not h5ad_files:
        st.warning("No `.h5ad` files found. Check your Folder ID and sharing permissions.")
        st.stop()

    selected_file_title = st.selectbox("Select a dataset:", options=list(h5ad_files.keys()))

    if selected_file_title:
        if 'loaded_file' not in st.session_state or st.session_state.loaded_file != selected_file_title:
            with st.spinner(f"Loading '{selected_file_title}'..."):
                file_id = h5ad_files[selected_file_title]
                temp_file_path = download_h5ad_from_drive(drive, file_id, selected_file_title)
                try:
                    st.session_state['adata'] = sc.read_h5ad(temp_file_path)
                    st.session_state['loaded_file'] = selected_file_title
                    st.success(f"Successfully loaded {st.session_state.adata.n_obs} cells from {st.session_state.adata.n_vars} genes.")
                finally:
                    if os.path.exists(temp_file_path):
                        os.remove(temp_file_path)

    # --- Plotting Section ---
    if 'adata' in st.session_state:
        adata = st.session_state.adata
        st.header(f"Visualizing: `{st.session_state.loaded_file}`")

        col1, col2 = st.columns(2)

        # ----- Column 1: Plot by Metadata -----
        with col1:
            st.subheader("Plot by Metadata")
            categorical_obs = adata.obs.select_dtypes(include=['category', 'object']).columns.tolist()
            if not categorical_obs:
                st.warning("No categorical metadata found in `adata.obs`.")
            else:
                color_by_meta = st.selectbox(
                    "Select metadata to color by:",
                    options=categorical_obs
                )
                fig1, ax1 = plt.subplots()
                sc.pl.umap(adata, color=color_by_meta, ax=ax1, show=False, legend_loc='on data')
                st.pyplot(fig1)

        # ----- Column 2: Plot by Gene Expression -----
        with col2:
            st.subheader("Plot by Gene Expression")
            gene_name = st.text_input("Enter a gene name (e.g., 'CD14', 'NKG7')", "").strip()
            
            if gene_name:
                if gene_name in adata.var_names:
                    with st.spinner(f"Plotting '{gene_name}'..."):
                        fig2, ax2 = plt.subplots()
                        sc.pl.umap(adata, color=gene_name, ax=ax2, show=False, 
                                   cmap='viridis',
                                   )
                        st.pyplot(fig2)
                else:
                    st.error(f"Gene '{gene_name}' not found in dataset. Please check the name (it's case-sensitive).")

except Exception as e:
    st.error(f"An error occurred: {e}")
