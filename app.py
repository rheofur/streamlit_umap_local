import streamlit as st
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os

# --- Page Configuration ---
st.set_page_config(
    page_title="Interactive UMAP Plotter",
    page_icon="🧬",
    layout="wide"
)

# --- App Title and Description ---
st.title("🔬 Interactive Single-Cell UMAP Plotter")
st.write("""
Upload your processed single-cell data in `.h5ad` format to generate an interactive UMAP plot.
Choose any categorical observation from your data to color the cells.
""")

# --- Sidebar for User Input ---
with st.sidebar:
    st.header("⚙️ Settings")

    # File Uploader
    uploaded_file = st.file_uploader("Upload your .h5ad file", type=["h5ad"])

    # Placeholder for AnnData object
    adata = None

    if uploaded_file is not None:
        try:
            # To read the file, we need to write it to a temporary location
            with open("temp.h5ad", "wb") as f:
                f.write(uploaded_file.getbuffer())
            
            adata = sc.read_h5ad("temp.h5ad")
            st.success("File successfully loaded!")
            
            # Clean up the temporary file
            os.remove("temp.h5ad")

        except Exception as e:
            st.error(f"Error loading file: {e}")
            st.stop() # Stop execution if file loading fails
    else:
        st.info("Awaiting file upload.")
        st.stop() # Stop execution until a file is uploaded


# --- Main Panel for Plotting ---
if adata:
    st.header("📊 UMAP Visualization")

    # Get categorical columns from adata.obs
    categorical_obs = adata.obs.select_dtypes(include=['category', 'object']).columns.tolist()

    if not categorical_obs:
        st.warning("No categorical data found in `adata.obs` to color the plot.")
        st.stop()

    # Dropdown to select the coloring category
    color_by = st.selectbox(
        "Color UMAP plot by:",
        options=categorical_obs,
        index=0,  # Default to the first option
        help="Select a column from your observation data (`adata.obs`) to color the cells."
    )

    # Button to generate the plot
    if st.button(f"Generate UMAP colored by {color_by}"):
        with st.spinner("Generating plot..."):
            
            # --- Plotting Logic ---
            # Create a matplotlib figure and axes
            fig, ax = plt.subplots(figsize=(8, 6))

            # Generate the UMAP plot using scanpy, passing the axes
            sc.pl.umap(adata, color=color_by, ax=ax, show=False, legend_loc='on data')
            
            # Display the plot in Streamlit
            st.pyplot(fig)

            # --- Dataframe Preview ---
            st.write("### Data Preview (`adata.obs`)")
            st.dataframe(adata.obs[[color_by]].head())
