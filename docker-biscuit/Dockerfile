FROM continuumio/miniconda3:latest

# download the latest available precompiled binary version of BISCUIT as
RUN apt-get update && apt-get install -y \
		curl &&\
		curl -OL $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest | grep browser_download_url | grep linux_amd64 | cut -d '"' -f 4) &&\
		chmod +x biscuit_* &&\
		mv biscuit_* /usr/local/bin/ &&\
		ln -s /usr/local/bin/biscuit_* /usr/local/bin/biscuit &&\
		conda create -n env -c bioconda -c conda-forge -c defaults python=3.6 \
		samblaster=0.1.24 \
		samtools=1.9 \
		pigz=2.3.4

ENV PATH /opt/conda/envs/env/bin:$PATH
