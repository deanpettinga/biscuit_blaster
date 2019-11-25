FROM continuumio/miniconda3:latest

# install unix dependencies
RUN apt-get update -y && \
    apt-get install git -y && \
    apt-get install build-essential -y && \
    apt-get install zlib1g-dev -y && \
    apt-get install libncurses5-dev -y && \
    apt-get install wget -y && \
    
    # intall biscuit
    cd / && \
    wget "https://github.com/zwdzwd/biscuit/releases/download/v0.3.8.20180515/biscuit_0_3_8_x86_64" && \
    cp biscuit_0_3_8_x86_64 /usr/bin/biscuit && \
    chmod 755 /usr/bin/biscuit && \
    
    # install samtools,samblaster,biscuit with conda
    conda create -n env -c bioconda -c conda-forge -c defaults python=3.6 \
	samblaster=0.1.24 \
	samtools=1.9 \
	pigz=2.3.4

ENV PATH /opt/conda/envs/env/bin:$PATH
VOLUME [ "/data" ]
