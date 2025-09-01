import streamlit as st
import pandas as pd
import io
import altair as alt
import umap.umap_ as umap
import scanpy as sc
from sklearn.preprocessing import StandardScaler
import json # <--- ADD THIS IMPORT

# --- Google Drive Authentication & Download Functions ---

from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

@st.cache_resource(ttl=600)
def authenticate_gdrive():
    """
    Uses Streamlit secrets to create credentials and build the 
    Google Drive service object. Caches the result for 10 minutes.
    """
    st.write("Authenticating with Google Drive...")
    try:
        # --- THIS PART IS UPDATED TO MATCH YOUR SECRETS ---
        # 1. Get the JSON string from your secrets
        json_credentials = st.secrets["google_credentials"]["service_account_json"]
        # 2. Parse the JSON string into a Python dictionary
        dict_credentials = json.loads(json_credentials)
        # 3. Use the dictionary to build the credentials
        creds = Credentials.from_service_account_info(
            dict_credentials,
            scopes=['https://www.googleapis.com/auth/drive.readonly']
        )
        # --- END OF UPDATED SECTION ---

        service = build('drive', 'v3', credentials=creds)
        st.write("Authentication successful.")
        return service
    except Exception as e:
        st.error(f"Google Drive authentication failed: {e}")
        return None

@st.cache_data(ttl=600)
def download_file_from_drive(_service, file_id):
    """
    Downloads a file from Google Drive into an in-memory bytes buffer.
    Caches the downloaded data for 10 minutes.
    The underscore in _service tells st.cache_data to ignore this argument 
    when checking for cache validity.
    """
    st.write(f"Downloading file with ID: {file_id}...")
    try:
        request = _service.files().get_media(fileId=file_id)
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
    file_id = '19ImkTu4j75f9CpHw-l0p2yQ0aW-a5t2z' 
    
    file_buffer = download_file_from_drive(drive_service, file_id)
    
    if file_buffer:
        try:
            adata_df = pd.read_csv(file_buffer) 
            st.success("✅ Successfully loaded data from Google Drive!")

            st.header("UMAP Configuration")
            
            col1, col2 = st.columns(2)
            with col1:
                n_neighbors = st.slider('Number of Neighbors (n_neighbors)', 2, 100, 15, 1)
                min_dist = st.slider('Minimum Distance (min_dist)', 0.0, 1.0, 0.1, 0.01)
            with col2:
                n_components = st.slider('Number of Components (n_components)', 1, 10, 2, 1)
                metric = st.selectbox('Metric', ['euclidean', 'manhattan', 'cosine', 'haversine'])
            
            st.header("Data Preview")
            st.dataframe(adata_df.head())
            
            adata = sc.AnnData(adata_df)
            
            scaler = StandardScaler()
            adata.X = scaler.fit_transform(adata.X)
            
            reducer = umap.UMAP(
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                n_components=n_components,
                metric=metric,
                random_state=42
            )
            embedding = reducer.fit_transform(adata.X)
            
            st.header("UMAP Projection")
            if n_components >= 2:
                umap_df = pd.DataFrame(embedding[:, :2], columns=['UMAP1', 'UMAP2'])
                chart = alt.Chart(umap_df).mark_circle(size=60).encode(
                    x=alt.X('UMAP1', scale=alt.Scale(zero=False)),
                    y=alt.Y('UMAP2', scale=alt.Scale(zero=False)),
                    tooltip=['UMAP1', 'UMAP2']
                ).properties(
                    width=700,
                    height=500
                ).interactive()
                st.altair_chart(chart, use_container_width=True)
            else:
                st.write("Set 'Number of Components' to 2 or more for a 2D plot.")

        except Exception as e:
            st.error(f"Failed to process the data. Error: {e}")
