FROM ubuntu:latest

ENV GFDL_WORK /tmp
ENV GFDL_BASE /isca
ENV GFDL_DATA /data
ENV GFDL_ENV ubuntu_conda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# # ignore missing hardware needed lfor openMPI speedup
# ENV ["OMPI_MCA_btl", "^openib"]
# # avoid mpi vader error [d3f8787e619d:05992] Read -1, expected 8192, errno = 1
# # https://github.com/open-mpi/ompi/issues/4948
# ENV  ["OMPI_MCA_btl_vader_single_copy_mechanism", "none"]

RUN apt-get update && apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive \
	apt-get install -y \
	build-essential \
	curl \
	gfortran \
	git \
    libnetcdf-dev \
    libpnetcdf-dev \
    libnetcdff-dev \
    libhdf5-openmpi-dev \
    python3 \
    python3-pip \
    tcl \
    tcl-dev \
    wget

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

# support for running as a local user so file permissions are correct
RUN curl -o /usr/local/bin/gosu -SL "https://github.com/tianon/gosu/releases/download/1.10/gosu-$(dpkg --print-architecture)" \
    && chmod +x /usr/local/bin/gosu

# creating a new user group called isca_build
RUN groupadd -g 9004 isca_build
# creating a new user that's part of the isca_build group called isca
RUN useradd -ms /bin/bash -u 9002 -g isca_build isca

#copy the current directory's code to the /isca directory in the docker image.
COPY . /isca
RUN ls /isca/.git
RUN chown -R isca /isca

RUN conda env create -f /isca/ci/environment-py3.9.yml
# RUN conda create --name isca_env_test xarray

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "isca_env", "/bin/bash", "-c"]

# Activate the environment, and make sure it's activated:
RUN echo "Make sure xarray is installed:"
RUN python -c "import xarray"

RUN pip install -e /isca/src/extra/python

RUN mkdir -p /data && chown -R isca /data
VOLUME /data
VOLUME /isca

WORKDIR /isca

CMD python -m pytest
# COPY docker-entrypoint.sh /usr/local/bin/entrypoint.sh
# RUN chmod +x /usr/local/bin/entrypoint.sh
# ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

#RUN python3 -c "import isca; cb = isca.IscaCodeBase.from_directory('/isca'); cb.compile()"
#CMD python3 /isca/exp/held_suarez.py --compile --up-to -i 3 -n 2
