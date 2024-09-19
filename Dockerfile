# A basic Dockerfile for building a Docker image for the Garfield tool

# Base image
FROM ubuntu:20.04

# maintainer
LABEL maintainer="your_email@example.com"

# disable interactive functions
ENV DEBIAN_FRONTEND=noninteractive

# system update and install necessary packages
RUN apt-get update && \
    apt-get install -y \
    perl \
    r-base \
    wget \
    curl \
    unzip \
    build-essential \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*



# newest: https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20240818.zip
# install plink
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20240818.zip && \
    unzip plink_linux_x86_64_20240818.zip -d /usr/local/bin/ && \
    rm plink_linux_x86_64_20240818.zip

# set cxx14 flag for `ranger` R-package
RUN mkdir -p ~/.R && \
echo "CXX = g++ -std=gnu++14" >> ~/.R/Makevars

# install R-packages
RUN R -e "install.packages('ranger', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('genio', repos='https://cloud.r-project.org/')"

# install logicFS from BiocManager
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org/'); BiocManager::install('logicFS')"

# copy the codes
COPY . /app

# set the perl library path
ENV PERL5LIB=/app/lib

# set working directory
WORKDIR /app

# run the installation script
RUN perl INSTALL.pl

# set the entrypoint
ENTRYPOINT ["perl", "Garfield"]
