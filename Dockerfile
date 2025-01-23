ARG BUSCO_VER="5.8.2"

FROM ubuntu:focal AS app

ARG BUSCO_VER
ARG BBMAP_VER="39.13"
ARG BLAST_VER="2.16.0"
ARG MINIPROT_VER="0.13"
ARG SEPP_VER="4.5.5"
ARG METAEUK_VER="7-bba0d80"
ARG DEBIAN_FRONTEND=noninteractive
ARG BUSCO_COMMIT=1590c5edd221266cf2252f6d8ad01eec567e9680

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="BUSCO"
LABEL software.version="${BUSCO_VER}"
LABEL description="Assessing genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs"
LABEL website="https://busco.ezlab.org/"
LABEL license="https://gitlab.com/ezlab/busco/-/raw/master/LICENSE"
LABEL maintainer="Kutluhan Incekara"
LABEL maintainer.email="kutluhan.incekara@ct.gov"

# install dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    apt-transport-https \
    ca-certificates \
    gnupg \
    curl \
    wget \
    python3-pip \
    python3-pandas \
    python3-setuptools\
    python3-requests \
    hmmer \
    prodigal \
    augustus \
    r-cran-ggplot2 \
    gcc-x86-64-linux-gnu \
    openjdk-8-jre-headless \
    libjenkins-json-java \
    libgoogle-gson-java \
    libjson-java \
    lbzip2 \
    && rm -rf /var/lib/apt/lists/* && apt-get autoclean \
    && ln -s /usr/bin/python3 /usr/bin/python

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    rm -rf /var/lib/apt/lists/*

# install other necessary tools
# BioPython (python3-biopython installs 1.73. It causes python error in this version)
RUN pip install --no-cache-dir biopython
# blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz &&\
    tar -xvf ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz && rm ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz
# sepp 
RUN wget https://github.com/smirarab/sepp/archive/refs/tags/v${SEPP_VER}.tar.gz &&\
    tar -xvf v${SEPP_VER}.tar.gz && rm v${SEPP_VER}.tar.gz &&\
    cd sepp-${SEPP_VER} &&\
    python setup.py config -c && python setup.py install
# bbtools
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VER}.tar.gz &&\
    tar -xvf BBMap_${BBMAP_VER}.tar.gz && rm BBMap_${BBMAP_VER}.tar.gz &&\    
    mv /bbmap/* /usr/local/bin/
# metaeuk
RUN wget https://github.com/soedinglab/metaeuk/releases/download/${METAEUK_VER}/metaeuk-linux-sse41.tar.gz &&\
    tar -xvf metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz &&\
    mv  /metaeuk/bin/* /usr/local/bin/
# miniprot
RUN wget https://github.com/lh3/miniprot/releases/download/v${MINIPROT_VER}/miniprot-${MINIPROT_VER}_x64-linux.tar.bz2 &&\
    tar -C /usr/local/bin/ --strip-components=1 --no-same-owner -xvf miniprot-${MINIPROT_VER}_x64-linux.tar.bz2 miniprot-${MINIPROT_VER}_x64-linux/miniprot &&\
    rm miniprot-${MINIPROT_VER}_x64-linux.tar.bz2

# and finally busco
RUN wget https://github.com/sage-wright/busco/archive/${BUSCO_COMMIT}.zip && \
    unzip ${BUSCO_COMMIT}.zip && \
    mv busco-${BUSCO_COMMIT} busco && \
    cd busco && \
    python -m pip install google-cloud-storage google-auth && \
    python -m pip install .

ENV AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config/" \
    PATH="${PATH}:/ncbi-blast-${BLAST_VER}+/bin:/usr/share/augustus/scripts:/busco/scripts" \
    LC_ALL=C

WORKDIR /data

CMD busco -h && generate_plot.py -h

## Tests ##
FROM app AS test

ARG BUSCO_VER

RUN busco -h && generate_plot.py -h