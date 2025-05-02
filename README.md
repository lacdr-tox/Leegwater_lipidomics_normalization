[![DOI](https://zenodo.org/badge/754549889.svg)](https://zenodo.org/doi/10.5281/zenodo.10635384)

# Lipidomics in breast cancer cell lines

This repository contains code to analyze data and generate figures for the following two manuscripts:

```         
"Leegwater, H., Zhang, Z., Zhang, X., Hankemeier, T., Harms, A. C., Zweemer, A. J. M., Le Dévédec, S. E., & Kindt, A. (2025). Normalization Strategies for Lipidome Data in Cell Line Panels. Journal of Chemometrics, 39(1), e3636. https://doi.org/10.1002/cem.3636
```

and soon code will be added for:

```         
"Distinct lipidomic profiles in breast cancer cell lines relate to proliferation and EMT phenotypes", by Leegwater, H. _et al._ (submitted).
```

If you already have questions about this new manuscript and the code, send us a message! We plan to release the code when the manuscript is accepted.

## Usage

Within the `code` folder, you will find R markdown files (Normalization manuscript) and Quarto files (Biological interpretation manuscript) and a `functions` folder. The R markdown/Quarto files can be used to rerun all code and to create all figures. Functions that one might want to reuse can be found in the functions folder.

Note that both projects were run with different versions of R and R packages. If you run into trouble with versions, the `docs` folder has the reports of the version used for manuscript submission, with all session info recorded.

### Data

Metabolomics data have been deposited to the EMBL-EBI MetaboLights ([`Yurekten et al., 2024`](https://doi.org/10.1093/nar/gkad1045)) with the identifier MTBLS9493 and is accessible at <https://www.ebi.ac.uk/metabolights/MTBLS9493>.

We hope to set all data to public soon, which will be when both papers are published. Then, we will create a new folder in this GitHub directory with the data in a separate `data` folder. For now, you can take a look at the HTML reports to see what the data and results could look like.

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

Thanks to [\@burgerga](https://www.github.com/burgerga) for suggestions on publishing this repository.
