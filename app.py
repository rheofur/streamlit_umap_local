import streamlit as st
import pandas as pd
import io
import altair as alt
import umap.umap_ as umap
import scanpy as sc
from sklearn.preprocessing import StandardScaler
import json

# --- Google Drive Authentication & File Handling Functions ---

from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

@st.cache_resource(ttl=600)
def authenticate_gdrive():
    """
    Uses Streamlit secrets to create credentials and build the 
    Google Drive service object.
    """
    # This function remains unchanged
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
    """
    Lists all non-folder files in a specific Google Drive folder.
    """
    # This function remains unchanged
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

# NOTE: The download function is now specific to the selected file,
# so we remove the global caching to ensure it always downloads the correct file.
def download_file_from_drive(service, file_id):
    """
    Downloads a file from Google Drive into an in-memory bytes buffer.
    """
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
st.title("🧬 UMAP Clustering from Google Drive Data")

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
            
            # --- CHANGE 1: ADD A BUTTON TO TRIGGER ANALYSIS ---
            # The download and analysis code is now inside this 'if' block.
            # It will only run when the user clicks the button.
            if st.button("Analyze Selected File"):
                if selected_filename:
                    selected_file_id = file_dict[selected_filename]
                    
                    file_buffer = download_file_from_drive(drive_service, selected_file_id)
        
                    if file_buffer:
                        try:
                            # --- CHANGE 2: USE sc.read_h5ad FOR .h5ad FILES ---
                            # This is the correct way to read AnnData's file format.
                            st.write("Reading .h5ad file...")
                            adata = sc.read_h5ad(filename=file_buffer)
                            st.success(f"✅ Successfully loaded AnnData object from `{selected_filename}`!")

                            st.header("UMAP Configuration")
                            col1, col2 = st.columns(2)
                            with col1:
                                n_neighbors = st.slider('Number of Neighbors (n_neighbors)', 2, 100, 15, 1)
                                min_dist = st.slider('Minimum Distance (min_dist)', 0.0, 1.0, 0.1, 0.01)
                            with col2:
                                n_components = st.slider('Number of Components (n_components)', 1, 10, 2, 1)
                                metric = st.selectbox('Metric', ['euclidean', 'manhattan', 'cosine', 'haversine'])
                            
                            st.header("Data Preview (Cell Metadata)")
                            # We display the observation metadata (adata.obs) which is a DataFrame.
                            st.dataframe(adata.obs.head())
                            
                            scaler = StandardScaler()
                            adata.X = scaler.fit_transform(adata.X)
                            
                            reducer = umap.UMAP(
                                n_neighbors=n_neighbors, min_dist=min_dist, 
                                n_components=n_components, metric=metric, random_state=42
                            )
                            embedding = reducer.fit_transform(adata.X)
                            
                            st.header("UMAP Projection")
                            if n_components >= 2:
                                umap_df = pd.DataFrame(embedding[:, :2], columns=['UMAP1', 'UMAP2'])
                                chart = alt.Chart(umap_df).mark_circle(size=60).encode(
                                    x=alt.X('UMAP1', scale=alt.Scale(zero=False)),
                                    y=alt.Y('UMAP2', scale=alt.Scale(zero=False)),
                                    tooltip=['UMAP1', 'UMAP2']
                                ).properties(width=700, height=500).interactive()
                                st.altair_chart(chart, use_container_width=True)
                            else:
                                st.write("Set 'Number of Components' to 2 or more for a 2D plot.")

                        except Exception as e:
                            st.error(f"Failed to process the data from {selected_filename}. Error: {e}")
