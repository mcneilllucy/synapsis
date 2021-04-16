FROM nfcore/base

LABEL authors="dmccarthy@svi.edu.au" \
    maintainer="Davis McCarthy <dmccarthy@svi.edu.au>" \
    description="Docker image containing all requirements for AAAA_2019_Project-Template"

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        git \
        libbz2-dev \
        zlib1g-dev \
        && rm -rf /var/lib/apt/lists/*

COPY environment.yml /
RUN conda env create -f /environment.yml python=3.6 && conda clean -a
ENV PATH /opt/conda/envs/aaaa-2021-Project-Template/bin:$PATH
