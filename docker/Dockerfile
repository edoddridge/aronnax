FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    gfortran \
    libopenmpi-dev \
    openmpi-bin \
    git \
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

RUN git clone --recursive https://github.com/edoddridge/aronnax.git

RUN cd aronnax/lib/hypre/src && ./configure && make install

RUN cd ../../../ && pip3 install -e .
