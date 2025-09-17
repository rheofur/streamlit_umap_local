# Interactive UMAP Visualizer
This is a Streamlit web application for interactively visualizing pre-processed single-cell data. You can upload your own data in the .h5ad format, filter it by cell metadata, and generate UMAP plots colored by metadata or gene expression.

# Requirements
- Python 3.8+
- pip

# Installation
## Download the code:
Download the application files (app.py, requirements.txt, etc.) and place them in a new folder on your computer.

## Navigate to the app directory:
Open your terminal or command prompt and use the cd command to go into the folder where you saved the files.
```
cd path/to/your/app-folder
```

## Install the required Python packages:
Run the following command to install all the necessary libraries from the requirements.txt file. We highly recommend using clean conda or python virtual environment.
```
# (Optional)
# conda create -n streamlit_umap python=3.8
# conda activate streamlit_ump
pip install -r requirements.txt
```

# How to Run the App
## Launch the Streamlit app:
In your terminal (while in the app's folder), run the following command:

```
streamlit run app.py --server.maxUploadSize 2000 # This sets maximum file size=2000MB=2GB
```
Of note, above "server.maxUploadSize" option sets maximum file size handled by upload. Since h5ad file will be loaded on RAM, adjust value appropriate for your file size (default maximum file size would be ~200MB)

## View in your browser:
The app should automatically open in a new tab in your web browser. If it doesn't, your terminal will provide a local URL (usually http://localhost:8501) that you can visit.

## Notes
Once the application is running, you can:

Upload your data: Use the file uploader in the sidebar to select and upload your .h5ad file.

Filter cells: Choose a metadata column and select the values you want to keep.

Generate plots: Select metadata columns or genes to color your UMAP plots.

Download plots: Click the "Download Plot" button to save your visualizations as PNG images.

Note: The application expects the uploaded .h5ad file to contain UMAP coordinates in adata.obsm['X_umap'].
