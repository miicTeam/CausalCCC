## CausalCCC

 This repository contains the source code for CausalCCC, a method based on MIIC (Multivariate Information based Inductive Causation) that learns a large class of causal or non-causal graphical models from purely observational data while including the effects of unobserved latent variables. Cell-Cell Communication (CCC) methods recover Ligand-Receptor interactions between single cell populations, but do not usually integrate intracellular causal pathways upstream and downstream of these L-R interactions. CausalCCC fills this gap by providing an integrative online server to reconstruct gene-gene interaction pathways across interacting cell types from single-cell transcriptomic data.  

### Abstract
CausalCCC is an interactive web server, which integrates causal discovery and cell-cell communication (CCC) methods to reconstruct gene interaction networks across interacting cell types, based on single-cell transcriptomic data and ligand-receptor pairs of interest. CausalCCC graphical outputs provide a unique integration of upstream and downstream signaling pathways across interacting cell types, which cannot be obtained by intracellular pathway analysis or CCC methods alone. CausalCCC is available as an R package or a [web server](https://miic.curie.fr/causalCCC.php)

### Input 

Users need to provide the raw count expression matrices for the sender and receiver cells as a text table (csv, tsv, txt, etc) with columns corresponding to selected genes and rows to single cells. A list of ligand-receptor (L-R) pairs of interest should also be provided to integrate CCC analysis. If needed, CausalCCC offers a wrapper function to seamlessly prepare all input files from a single-cell object (Seurat or Anndata) and a selection of CCC methods (see below). Tutorials and guidance are provided in the Quick Start and Advanced Mode sections. No preprocessed or large single-cell object need to be uploaded. 

### Output

The output is a causal network which elucidates how intracellular signaling upstream of ligand expression in sender cells can lead to receptor activation and downstream signaling in receiver cells. The CausalCCC server provides extensive, interactive data visualization tools within the browser, allowing users to explore their CausalCCC networks in detail. For guidance on interpreting network results, we refer to the “How to interpret a network” section in the tutorial web page. 

### Processing method 

CausalCCC methodology integrates a robust and scalable causal network reconstruction method, MIIC (Multivariate Information-based Inductive Causation), validated in several peer-reviewed publications (REFS), with user-defined or internally computed ligand-receptor pairs from multiple CCC methodologies, such as CellphoneDBv5, NATMI, iTALK, Log2FC (Liana+), SCA, and Connectome (DB2020) (REFS). If needed, CausalCCC input files can be readily obtained with a single wrapper compatible with Seurat objects in R and Anndata objects in Python. CausalCCC includes a demo dataset within the workbench page and offers comprehensive tutorials and support materials. Compared to an R package, CausalCCC offers enhanced computational power along with interactive network visualization, delivering a clear advantage for researchers using single-cell objects. 

## Prerequisites
MIIC contains R and C++ sources.
- To compile from source, a compiler with support for c++14 language features is required.
- MIIC imports the following R packages: ppcor, scales, stats, Rcpp

## Installation

You need both the MIIC R package version catering for CausalCCC and the CausalCCC package containing pre-processing automatic functions. From GitHub only (for now):
```R
# install.packages("devtools")
devtools::install_github("miicTeam/miic_R_package@causalccc", force = T)
devtools::install_github("miicTeam/CausalCCC", force = T)

```

## Quick start

Please see our detailed documentation in the vignettes or directly on the webserver [tutorials](https://miic.curie.fr/tutorial_causalCCC.php)

## Documentation
You can find the documentation pages in the "man" folder, or use R functions `help()` and `?`. Please see our [wiki](https://...) and [tutorials](https://www....)

For python lovers, please refer to our preprocessing vignette [here]. 

## Authors
- Louise Dupuis, Orianne Debeaupuis
- Franck Simon
- Hervé Isambert

## License
GPL-2 | GPL-3

## Cite this work


## References

Simon F., Comes M. C., Tocci T., Dupuis L., Cabeli V., Lagrange N., Mencattini A., Parrini M. C., Martinelli E., Isambert H.; [CausalXtract: a flexible pipeline to extract causal effects from live-cell time-lapse imaging data; eLife, reviewed preprint](https://www.biorxiv.org/content/10.1101/2024.02.06.579177v1.abstract).

Ribeiro-Dantas M. D. C., Li H., Cabeli V., Dupuis L., Simon F., Hettal L., Hamy A. S., Isambert H.; [Learning interpretable causal networks from very large datasets, application to 400,000 medical records of breast cancer patients; iScience, 2024](https://arxiv.org/abs/2303.06423).

Cabeli V., Li H., Ribeiro-Dantas M., Simon F., Isambert H.; [Reliable causal discovery based on mutual information supremum principle for finite dataset; Why21 at NeurIPS, 2021](https://why21.causalai.net/papers/WHY21_24.pdf).

Cabeli V., Verny L., Sella N., Uguzzoni G., Verny M., Isambert H.; Learning clinical networks from medical records based on information estimates in mixed-type data; PLoS computational biology., 2020. [doi:10.1371/journal.pcbi.1007866](https://doi.org/10.1371/journal.pcbi.1007866) | [code](https://github.com/vcabeli/miic_PLoS)

Li H., Cabeli V., Sella N., Isambert H.; [Constraint-based causal structure learning with consistent separating sets; In Advances in Neural Information Processing Systems 2019.](https://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets) | [code](https://github.com/honghaoli42/consistent_pcalg)

Verny L., Sella N., Affeldt S., Singh PP., Isambert H.; Learning causal networks with latent variables from multivariate information in genomic data;  PLoS Comput. Biol., 2017. [doi:10.1371/journal.pcbi.1005662](https://doi.org/10.1371/journal.pcbi.1005662)
