import streamlit as st
import pandas as pd
import io
import scanpy as sc
import json
import matplotlib.pyplot as plt

# --- Google Drive Authentication & File Handling Functions ---
# (These functions remain the same as the previous version)

from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

@st.cache_resource(ttl=600)
def authenticate_gdrive():
    st.write("Authenticating with Google Drive...")
    try:
        json_credentials = st.secrets["google_credentials"]["service_account_json"]
        dict_credentials = json.loads(json_credentials)
        creds = Credentials.from_service_account_info(
            dict_credentials,
            scopes=['https://www.googleapis.com/auth/drive.readonly']
        )
        service = build('drive', 'v3', credentials=creds)
        st.write("Authentication successful.")
        return service
    except Exception as e:
        st.error(f"Google Drive authentication failed: {e}")
        return None

@st.cache_data(ttl=600)
def list_files_in_folder(_service, folder_id):
    st.write(f"Fetching file list from folder...")
    try:
        files_list = []
        page_token = None
        while True:
            response = _service.files().list(
                q=f"'{folder_id}' in parents and mimeType != 'application/vnd.google-apps.folder'",
                spaces='drive',
                fields='nextPageToken, files(id, name)',
                pageToken=page_token
            ).execute()
            
            files_list.extend(response.get('files', []))
            page_token = response.get('nextPageToken', None)
            if page_token is None:
                break
        st.write(f"Found {len(files_list)} file(s).")
        return files_list
    except Exception as e:
        st.error(f"An error occurred while listing files: {e}")
        return []

def download_file_from_drive(service, file_id):
    st.write(f"Downloading selected file...")
    try:
        request = service.files().get_media(fileId=file_id)
        file_io = io.BytesIO()
        downloader = MediaIoBaseDownload(file_io, request)
        
        done = False
        while not done:
            status, done = downloader.next_chunk()
        
        file_io.seek(0)
        st.write("File download complete.")
        return file_io
    except Exception as e:
        st.error(f"An error occurred while downloading the file: {e}")
        return None

# --- Main Streamlit App ---

st.set_page_config(layout="wide")
st.title("🔬 Interactive UMAP Visualizer for Pre-processed Data")

if 'gene' not in st.session_state:
    st.session_state.gene = ''

drive_service = authenticate_gdrive()

if drive_service:
    folder_id = st.secrets.get("folder_id")
    if not folder_id:
        st.error("`folder_id` not found in your secrets.toml file.")
    else:
        files = list_files_in_folder(drive_service, folder_id)
        
        if not files:
            st.warning(f"No files found in the specified Google Drive folder (ID: {folder_id}).")
        else:
            file_dict = {file['name']: file['id'] for file in files}
            
            selected_filename = st.selectbox(
                "Choose a file from your Google Drive folder:",
                options=list(file_dict.keys())
            )
            
            if st.button("Load and Visualize File"):
                if selected_filename:
                    selected_file_id = file_dict[selected_filename]
                    
                    with st.spinner(f"Loading `{selected_filename}`..."):
                        file_buffer = download_file_from_drive(drive_service, selected_file_id)
                        if file_buffer:
                            adata = sc.read_h5ad(filename=file_buffer)
                            st.session_state.adata = adata
                            st.success(f"✅ Successfully loaded `{selected_filename}`!")

if 'adata' in st.session_state:
    adata = st.session_state.adata
    
    st.header("🎨 UMAP Visualization")
    
    if 'X_umap' not in adata.obsm:
        st.error("Error: The loaded file does not contain UMAP coordinates in `.obsm['X_umap']`.")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            # --- CHANGE: REMOVED CURATED LIST ---
            # The dropdown is now populated with all columns from the .obs dataframe.
            available_obs = list(adata.obs.columns)
            
            color_by_meta = st.selectbox(
                "Color by cell metadata:",
                options=available_obs
            )
            
            dot_size = st.slider("Dot size:", min_value=1, max_value=100, value=15)
            
            st.write(f"Generating UMAP plot colored by **{color_by_meta}**...")
            with st.spinner('Generating plot...'):
                fig_meta = sc.pl.umap(adata, color=[color_by_meta], 
                                      frameon=False, size=dot_size, 
                                      return_fig=True, show=False)
                st.pyplot(fig_meta)

        with col2:
            st.session_state.gene = st.text_input(
                "Search for a gene to color by:",
                value=st.session_state.gene
            )
            
            if st.session_state.gene:
                gene_to_plot = st.session_state.gene.strip()
                if gene_to_plot in adata.var_names:
                    st.write(f"Generating UMAP plot for gene **{gene_to_plot}**...")
                    with st.spinner('Generating plot...'):
                        fig_gene = sc.pl.umap(adata, color=[gene_to_plot], 
                                              frameon=False, size=dot_size,
                                              return_fig=True, show=False,
                                              cmap='viridis')
                        st.pyplot(fig_gene)
                else:
                    st.warning(f"Gene '{gene_to_plot}' not found in the dataset.")

        with st.expander("Show Data Preview"):
            st.subheader("Cell Metadata (`.obs`)")
            st.dataframe(adata.obs.head())
            st.subheader("Gene Metadata (`.var`)")
            st.dataframe(adata.var.head())
