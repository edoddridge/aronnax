FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    gfortran \
    libopenmpi-dev \
    openmpi-bin

RUN pip install \
    scipy \
    numpy \
    pandas \
    xarray \
    toolz \
    dask \
    future \
    matplotlib \
    pytest

RUN cd lib/hypre/src && ./configure && make install

RUN cd ../../../ && pip install -e .
