# PROTOCOL
- Author: Mauricio Morales
- Rules
    - Record desision making in this `PROTOCOL` file.
    - Reports must be in `.qmd` or `.ipynb` files.
    - Libraries must be handeled using [uv](https://docs.astral.sh/uv/) and [uvr](https://github.com/nbafrank/uvr).
    - For plotting purposes, use `tiff` and `pdf` formats. Keep `pdf` format in `GitHub`.
    - Focus `OOP` workflows.

## Jun 14, 2026

- Adapting project layout for downstream analysis. 
- Configurating `uv` and `uvr` options for the necesary libraries in this project (e.g. [Seurat](https://satijalab.org/seurat/) and [Scanpy](https://scanpy.scverse.org/en/stable/)).

## Jun 18, 2026

- Selecting dependencies for the project.
- Integrating `renv` library for enviroment managment and sharing options.

## Jun 19, 2026

- Setting enviroment to work with Docker-renv (exploration): [Using renv with Docker](https://rstudio.github.io/renv/articles/docker.html).
    - Not necessary for now, at least until the project is done.
- Seurat integration for [Single Cell Analysis: Seurat](https://satijalab.org/seurat/).