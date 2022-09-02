# psm_utils online

psm_utils online is a Streamlit-based web server, built on top of the psm_utils
Python package. It allows you to easily get proteomics peptide-spectrum match
(PSM) statistics for any supported PSM file type, and to convert search engine
results from one PSM file format into another.


## Development setup

Install with pip in editable mode, using the `dev` and `online` optional
dependencies:

```sh
pip install --editable .[dev,online]
```

Start a local Streamlit web server:

```sh
streamlit run .\online\Home.py
```