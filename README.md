**BUSCOv4 - Benchmarking sets of Universal Single-Copy Orthologs.**

Main changes in v4:

- Automated selection of lineages issued from https://www.orthodb.org/ release 10

- Automated download of all necessary files and datasets to conduct a run

- Use prodigal for non-eukaryotic genomes


To install, clone the repository and enter ``sudo python setup.py install`` or ``python setup.py install --user``

Or use the offical Docker container: https://hub.docker.com/r/ezlabgva/busco/tags

To get help, ``busco -h`` and ``python scripts/generate_plot.py -h``

Do not forget to create a ``config.ini`` file in the ``config/`` subfolder (the script busco_configurator.py can help you filling it). You can set the ``BUSCO_CONFIG_FILE`` 
environment variable to define a custom path (including the filename) to the ``config.ini`` file, 
useful for switching between configurations or in a multi-user environment.

**!!!** Don't use "odb9" datasets with BUSCOv4. If you need to reproduce previous analyses, use BUSCOv3 (https://gitlab.com/ezlab/busco/-/tags/3.0.2)

See also the user guide: https://busco.ezlab.org/busco_userguide.html
and report problems on the BUSCO issue board at https://gitlab.com/ezlab/busco/issues

**How to cite BUSCO**

*BUSCO: Assessing Genome Assembly and Annotation Completeness.*
Mathieu Seppey, Mosè Manni, Evgeny M. Zdobnov
In: Kollmar M. (eds) Gene Prediction. Methods in Molecular Biology, vol 1962. Humana, New York, NY. 2019
doi.org/10.1007/978-1-4939-9173-0_14

*BUSCO applications from quality assessments to gene prediction and phylogenomics.*
Robert M. Waterhouse, Mathieu Seppey, Felipe A. Simão, Mosè Manni, Panagiotis Ioannidis, Guennadi Klioutchnikov, Evgenia V. Kriventseva, and Evgeny M. Zdobnov
*Mol Biol Evol*, published online Dec 6, 2017 
doi: 10.1093/molbev/msx319 

*BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.*
Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov
*Bioinformatics*, published online June 9, 2015 
doi: 10.1093/bioinformatics/btv351

Copyright (c) 2016-2019, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.
