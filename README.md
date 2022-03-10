# [IRscope+](https://irscope.shinyapps.io/irscope-main/) - Updated Inverted Repeat scoping tool

![Website](https://img.shields.io/website?url=https%3A%2F%2Firscope.shinyapps.io%2Firscope-main%2F)

Link to the [web application](https://irscope.shinyapps.io/irscope-main/).

## Table of Contents

- [About the program and how to use it](#about-the-program-and-how-to-use-it)
  * [How the project is divided](#how-the-project-is-divided)
- [Questions](#questions)
- [Reference](#reference)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)

## About the program and how to use it

IIRscope is a tool for visualizing the genes and indel mismatches on the boundaries of the junction sites of the chloroplast genome. 

You can either use the web application at [https://irscope.shinyapps.io/irscope-main/]() or download the code from the repository and run it locally. In the latter case you will need to use [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/). You need to open the file [`app.R`](https://github.com/AmiryousefiLab/IIRscope/blob/main/app.R) and press *Run App*. A GUI equal to the one at the web will appear and you can test everything manually.

### How the project is divided

There are two main parts:

- **Related to the Shiny app**: `app` file to set up the UI and server for the web application.
- **IRscope R package**: in *IRscope/* folder. It wraps all the functions used by the program to find the IR, mismatches and for the plotting. It contains several files:
    - `imports`: libraries used by the package.
    - `read_gb_file`: functions to fetch or read the gbfiles.
    - `parse_gb_file`: given the gb file or dogma file extracts the needed information.
    - `process_data`: main functions to calculate all the information needed for the plotting (`IRs` and `IRsD`).
    - `gene_info`: gets the information about the genes in the needed format.
    - `detect_ir`: functions to detect the IR and indel mismatches.
    - `plot_genome`: main functions to plot (`IRs2` and `IRs2D`) and the needed auxiliar functions.

## Questions

You can access the *FAQ* section from the *Questions?* section in the [web app](https://irscope.shinyapps.io/irscope-main/). If you don't find what you're looking for there, you can post an issue [here](https://github.com/AmiryousefiLab/IIRscope/issues).

## Reference 

Publication at Bioinformatics <https://academic.oup.com/bioinformatics/article/34/17/3030/4961430>.

## Contact

- Ali Amiryousefi - ali.amiryousefi@helsinki.fi
- Carmen Diez - carmen.diezmenendez@estudiante.uam.es

## Acknowledgements

* [Chloroplot](https://github.com/shuyuzheng/Chloroplot/): visualizing organelle genomes.
* University of Helsinki.
* [Img Shields](https://shields.io): web working symbol.
* [Freepik icons](https://www.flaticon.com/authors/freepik): SVG icons in *Coloring* section from [Flaticon](https://www.flaticon.com/).
