# FROM alpine:3.10
FROM python:3

RUN mkdir -p /src

ADD ./docker-build/ViennaRNA-2.4.12.tar.gz /src
ADD ./docker-build/mfold-3.6.tar.gz /src

# RUN apk add python3 \
#       gcc g++ make zlib \
#       gfortran libgfortran 
# RUN yum install -y python3 gcc \
#       gcc-c++ kernel-devel make \
#       gcc-gfortran libgfortran \
#       which git groff

RUN apt-get update && \
      apt-get install -y --no-install-recommends \
      gfortran groff

WORKDIR /src

RUN cd mfold-3.6 && \
      ./configure && \
      make && \
      make install

RUN cd /src/ViennaRNA-2.4.12 && \
      ./configure --disable-debug --without-perl && \
      make && \
      make install && \
      make check

RUN mkdir -p tests
COPY ./samples/Trunc_Test.fa tests/trunc.fa

RUN apt-get update && \
      apt-get \
            install -y --no-install-recommends \
            texlive-font-utils ghostscript

RUN python -m venv /aptamer && \
      . /aptamer/bin/activate && \
      pip install --upgrade pip && \
      pip install numpy \
            scipy \
            biopython \
            python-Levenshtein

RUN mkdir -p /aptamer/lib/aptamer
COPY ./aptamer /aptamer/lib/aptamer
COPY ./create-graph /aptamer/bin/create-graph
COPY ./predict-structures /aptamer/bin/predict-structures

RUN mkdir -p /files && mv /src/tests/trunc.fa /files/test.fa

VOLUME ["/data"]
VOLUME ["/files"]

WORKDIR "/data"

ENV VIRTUAL_ENV /aptamer
ENV PYTHONPATH /usr/local/lib/python3.7/site-packages
ENV PATH "/aptamer/bin:$PATH"
ENV SCRIPT_PATH /aptamer/lib/aptamer
ENV DATA_DIR /data
SHELL ["/bin/bash"]
CMD ["bash"]

