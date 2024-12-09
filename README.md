[![DOI](https://zenodo.org/badge/754549889.svg)](https://zenodo.org/doi/10.5281/zenodo.10635384)

# Normalization Strategies for Lipidome Data in Cell Line Panels

Code to generate figures for "Normalization Strategies for Lipidome Data in Cell Line Panels", by Leegwater, H., Zhang, Z., Zhang, X., Hankemeier, T., Harms, A., Zweemer, A., Le Dévédec, S. and Kindt, A. (2024). Journal of Chemometrics e3636. <https://doi.org/10.1002/cem.3636> .

## Usage

Within the code folder, you will find two R markdown files and a functions folder. The R markdown files can be used to rerun all code and to create all figures. Functions that one might want to reuse can be found in the functions folder.

### Data

Metabolomics data have been deposited to the EMBL-EBI MetaboLights ([`Yurekten et al., 2024`](https://doi.org/10.1093/nar/gkad1045)) with the identifier MTBLS9493and is accessible at https://www.ebi.ac.uk/metabolights/MTBLS9493. 

When this dataset status is changed to public, we will add it and metadata to the data folder in this repository. For now, you can take a look at the html reports to see what the data could look like.

### Figures

Figures are generated reproducibly in R using [`renv`](https://rstudio.github.io/renv/index.html):

1.  Download/clone this repository

2.  Open the project file (`.Rproj`) in RStudio

3.  Run

    ``` r
    renv::restore()
    ```

    to install R package dependencies.

4.  Open `code/Calculations_to_get_all_data_for_all_figures.Rmd` and choose *Run* \> *Run All*. You may need to set a custom directory for the data, since this is not yet part of this repository.

5.  Open `code/All_figures.Rmd` and choose *Run* \> *Run All*. Figures will appear in the specified `output_dir` folder.

## Acknowledgments

Thanks to [@burgerga](https://www.github.com/burgerga) for suggestions on archiving this repository.
