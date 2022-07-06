## BUSCOv5 - Benchmarking sets of Universal Single-Copy Orthologs.

For full documentation please consult the user guide: https://busco.ezlab.org/busco_userguide.html

Main changes in v5:

- Metaeuk is used as default gene predictor for eukaryote pipeline. Augustus is maintained and can be used optionally instead of Metaeuk.
- Introduction of batch mode: input argument can be a folder containing input files
- The folder structure has changed, so if doing a manual installation, make sure to completely remove any previous versions of BUSCO before installing v5.

***
### Installation

#### Conda
Conda installation instructions are in the userguide here:
https://busco.ezlab.org/busco_userguide.html#conda-package

#### Docker
BUSCO is available through DockerHub - instructions here:
https://busco.ezlab.org/busco_userguide.html#docker-image

#### Manual installation
Manual installation is possible, though it is important to validate each of the dependencies before running BUSCO. 
More details in the user guide: https://busco.ezlab.org/busco_userguide.html#manual-installation

***
### Troubleshooting
To get help with BUSCO use: ``busco -h`` and ``python3 scripts/generate_plot.py -h``

Report problems on the BUSCO issue board at https://gitlab.com/ezlab/busco/issues

***
### How to cite BUSCO

**If you have used BUSCO v4 or v5 in your analyses, please cite:**

Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, [BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes](https://academic.oup.com/mbe/article/38/10/4647/6329644). Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654

Additional protocols and applications are described in:  
Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). [BUSCO: Assessing genomic data quality and beyond](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.323). Current Protocols, 1, e323. doi: 10.1002/cpz1.323


**Previous publications:**   
BUSCO: Assessing Genome Assembly and Annotation Completeness. Mathieu Seppey, Mosè Manni, Evgeny M. Zdobnov In: Kollmar M. (eds) Gene Prediction. Methods in Molecular Biology, vol 1962. Humana, New York, NY. 2019 doi.org/10.1007/978-1-4939-9173-0_14

BUSCO applications from quality assessments to gene prediction and phylogenomics. Robert M. Waterhouse, Mathieu Seppey, Felipe A. Simão, Mosè Manni, Panagiotis Ioannidis, Guennadi Klioutchnikov, Evgenia V. Kriventseva, and Evgeny M. Zdobnov Mol Biol Evol, published online Dec 6, 2017 doi: 10.1093/molbev/msx319

BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov Bioinformatics, published online June 9, 2015 doi: 10.1093/bioinformatics/btv351

Copyright (c) 2016-2022, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.
