**BUSCOv5.beta - Benchmarking sets of Universal Single-Copy Orthologs.**

BUSCO v5.beta is available for testing. We welcome feedback from users.
Report problems on the BUSCO issue board at https://gitlab.com/ezlab/busco/issues.

Main changes in v5:

- Replace Augustus gene predictor with Metaeuk (https://github.com/soedinglab/metaeuk)

**BUSCOv4 - Benchmarking sets of Universal Single-Copy Orthologs.**

Last stable version: https://gitlab.com/ezlab/busco/-/tags/4.1.2

For full documentation please consult the user guide: https://busco.ezlab.org/busco_userguide.html



To install, clone the repository and enter ``sudo python3 setup.py install`` or ``python3 setup.py install --user``

More details in the user guide: https://busco.ezlab.org/busco_userguide.html#manual-installation

Do not forget to edit the ``config/config.ini`` file to match your environment. The script `scripts/busco_configurator.py` can help with this. 
You can set the ``BUSCO_CONFIG_FILE`` environment variable to define the path (including the filename) to that ``config.ini`` file. 

```
export BUSCO_CONFIG_FILE="/path/to/myconfig.ini"
```
Alternatively you can pass the config file path as a command line argument using ``--config /path/to/config.ini``.


If you have trouble installing one of the many third-party tools, try the official Docker container: https://hub.docker.com/r/ezlabgva/busco/tags

To get help with BUSCO use: ``busco -h`` and ``python3 scripts/generate_plot.py -h``

**!!!** Do not use "odb9" datasets with BUSCOv4. If you need to reproduce previous analyses, use BUSCOv3 (https://gitlab.com/ezlab/busco/-/tags/3.0.2)

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

Copyright (c) 2016-2020, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.
