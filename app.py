import streamlit as st
import pandas as pd
import io
import scanpy as sc
import matplotlib.pyplot as plt

# --- Page Configuration ---
st.set_page_config(layout="wide")

# --- Main Application Logic ---
def main_app():
    """
    Contains the core logic of the Streamlit app.
    """
    st.title("🔬 Interactive UMAP Visualizer for Local Data")

    # --- Sidebar for Controls ---
    st.sidebar.header("Controls")
    
    # --- Local File Uploader ---
    uploaded_file = st.sidebar.file_uploader(
        "1. Upload your pre-processed data (.h5ad)",
        type=['h5ad']
    )

    if uploaded_file is not None:
        with st.spinner(f"Loading `{uploaded_file.name}`..."):
            # When a new file is uploaded, it replaces the old data.
            adata = sc.read_h5ad(filename=uploaded_file)
            st.session_state.original_adata = adata
            st.session_state.adata = adata # This is the one we will filter and plot
            st.sidebar.success(f"✅ Successfully loaded `{uploaded_file.name}`!")

    # --- Filtering and Visualization Logic (largely unchanged) ---
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

    # --- Main Visualization Area ---
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
main_app()
