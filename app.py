import streamlit as st
import pandas as pd
import io
import scanpy as sc
import json
import matplotlib.pyplot as plt
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

# --- Page Configuration ---
st.set_page_config(layout="wide")

# --- Password Protection ---
def check_password():
    """
    Returns True if the user has entered the correct password.
    Otherwise, returns False.
    """
    def password_entered():
        """Checks whether the password entered by the user is correct."""
        if st.session_state["password"] == st.secrets["password_code"]:
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # Don't store password in session state.
        else:
            st.session_state["password_correct"] = False

    # Initialize session state if it's the first run
    if "password_correct" not in st.session_state:
        st.session_state["password_correct"] = False

    # Show the password input field if the password is not correct
    if not st.session_state["password_correct"]:
        st.header("🔑 Secure Access")
        st.write("Please enter the password to use this application.")
        st.text_input(
            "Password", type="password", on_change=password_entered, key="password"
        )
        if "password" in st.session_state and not st.session_state["password_correct"]:
            st.error("😕 Incorrect password. Please try again.")
        return False
    else:
        return True

# --- Main Application Logic ---
def main_app():
    """
    Contains the core logic of the Streamlit app.
    This function only runs after the user has been authenticated.
    """
    st.title("🔬 Interactive UMAP Visualizer for Pre-processed Data")

    # --- Google Drive Authentication & File Handling Functions ---
    @st.cache_resource(ttl=600)
    def authenticate_gdrive():
        try:
            json_credentials = st.secrets["google_credentials"]["service_account_json"]
            dict_credentials = json.loads(json_credentials)
            creds = Credentials.from_service_account_info(
                dict_credentials,
                scopes=['https://www.googleapis.com/auth/drive.readonly']
            )
            service = build('drive', 'v3', credentials=creds)
            return service
        except Exception as e:
            st.error(f"Google Drive authentication failed: {e}")
            return None

    @st.cache_data(ttl=600)
    def list_files_in_folder(_service, folder_id):
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
            return files_list
        except Exception as e:
            st.error(f"An error occurred while listing files: {e}")
            return []

    def download_file_from_drive(service, file_id):
        try:
            request = service.files().get_media(fileId=file_id)
            file_io = io.BytesIO()
            downloader = MediaIoBaseDownload(file_io, request)
            done = False
            while not done:
                status, done = downloader.next_chunk()
            file_io.seek(0)
            return file_io
        except Exception as e:
            st.error(f"An error occurred while downloading the file: {e}")
            return None

    # --- Sidebar for Controls ---
    st.sidebar.header("Controls")
    drive_service = authenticate_gdrive()

    if drive_service:
        folder_id = st.secrets.get("folder_id")
        if not folder_id:
            st.sidebar.error("`folder_id` not found in your secrets.toml file.")
        else:
            files = list_files_in_folder(drive_service, folder_id)
            if not files:
                st.sidebar.warning("No files found in the folder.")
            else:
                file_dict = {file['name']: file['id'] for file in files}
                selected_filename = st.sidebar.selectbox(
                    "1. Choose a file:",
                    options=list(file_dict.keys())
                )
                
                if st.sidebar.button("Load Data"):
                    with st.spinner(f"Loading `{selected_filename}`..."):
                        selected_file_id = file_dict[selected_filename]
                        file_buffer = download_file_from_drive(drive_service, selected_file_id)
                        if file_buffer:
                            # When new data is loaded, it replaces the old data in session_state.
                            # Python's garbage collector then frees the memory of the old object.
                            # This fulfills the requirement to "free" the previous file's data.
                            adata = sc.read_h5ad(filename=file_buffer)
                            st.session_state.original_adata = adata
                            st.session_state.adata = adata # This is the one we will filter and plot
                            st.success(f"✅ Successfully loaded `{selected_filename}`!")

    if 'original_adata' in st.session_state:
        st.sidebar.markdown("---")
        st.sidebar.header("Filtering Options")
        
        adata_original = st.session_state.original_adata
        
        filter_column = st.sidebar.selectbox(
            "2. Filter by metadata column:",
            options=["None"] + list(adata_original.obs.columns)
        )

        if filter_column != "None":
            unique_values = list(adata_original.obs[filter_column].unique())
            selected_values = st.sidebar.multiselect(
                f"3. Select values from '{filter_column}' to keep:",
                options=unique_values,
                default=unique_values
            )
            
            adata_filtered = adata_original[adata_original.obs[filter_column].isin(selected_values)]
            st.session_state.adata = adata_filtered
        else:
            st.session_state.adata = st.session_state.original_adata
            
        if st.sidebar.button("Reset All Filters"):
            st.session_state.adata = st.session_state.original_adata
            st.rerun()

    # --- Main Visualization Area (Plotting logic is unchanged) ---
    if 'adata' in st.session_state:
        adata = st.session_state.adata
        
        st.info(f"Displaying **{len(adata)} cells** out of {len(st.session_state.original_adata)} total cells.")
        st.header("🎨 UMAP Visualization")
        
        if 'X_umap' not in adata.obsm:
            st.error("Error: Loaded file has no UMAP coordinates in `.obsm['X_umap']`.")
        else:
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Color by Cell Metadata")
                available_obs = list(adata.obs.columns)
                color_by_meta = st.selectbox("Select metadata column:", available_obs)
                dot_size_meta = st.slider("Dot size:", 1, 100, 15, key="meta_dot_size")
                
                with st.spinner('Generating plot...'):
                    fig_meta = sc.pl.umap(adata, color=[color_by_meta], frameon=False, size=dot_size_meta, return_fig=True, show=False)
                    st.pyplot(fig_meta)
                    
                    buf = io.BytesIO()
                    fig_meta.savefig(buf, format="png", bbox_inches='tight', dpi=300)
                    st.download_button("Download Plot", buf.getvalue(), f"umap_{color_by_meta}.png", "image/png")

            with col2:
                st.subheader("Color by Gene Expression")
                genes_to_plot = st.multiselect("Search and select gene(s):", list(adata.var_names))
                
                if genes_to_plot:
                    dot_size_gene = st.slider("Dot size:", 1, 100, 15, key="gene_dot_size")
                    with st.spinner('Generating plot(s)...'):
                        fig_gene = sc.pl.umap(adata, color=genes_to_plot, frameon=False, size=dot_size_gene, return_fig=True, show=False, cmap='viridis')
                        st.pyplot(fig_gene)
                        
                        buf = io.BytesIO()
                        fig_gene.savefig(buf, format="png", bbox_inches='tight', dpi=300)
                        gene_filename = "_".join(genes_to_plot)
                        st.download_button("Download Plot", buf.getvalue(), f"umap_gene_{gene_filename}.png", "image/png")

            with st.expander("Show Data Preview"):
                st.subheader("Filtered Cell Metadata (`.obs`)")
                st.dataframe(adata.obs)


# --- App Entry Point ---
# First, run the password check. If it passes, run the main application.
if check_password():
    main_app()

