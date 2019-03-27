FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    gfortran \
    libopenmpi-dev \
    openmpi-bin \
    python3.6 \
    python3-pip

RUN pip3 install \
    scipy \
    numpy \
    pandas \
    xarray \
    toolz \
    dask \
    future \
    matplotlib \
    pytest

RUN pwd

RUN ls

RUN cd /root/lib/hypre/src && ./configure && make install

RUN cd ../../../ && pip install -e .
