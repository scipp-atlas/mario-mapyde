ARG BUILDER_IMAGE=python:3.8-slim
FROM ${BUILDER_IMAGE} as builder

USER root
WORKDIR /usr/local

SHELL [ "/bin/bash", "-c" ]

RUN apt-get -qq -y update && \
    apt-get -qq -y install \
      gcc \
      g++ \
      gfortran \
      make \
      vim \
      zlibc \
      zlib1g-dev \
      libbz2-dev \
      rsync \
      bash-completion \
      python3-dev \
      less \
      wget && \
    apt-get -y autoclean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

# Install LHAPDF
ARG LHAPDF_VERSION=6.3.0
RUN mkdir /code && \
    cd /code && \
    wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-${LHAPDF_VERSION}.tar.gz -O LHAPDF-${LHAPDF_VERSION}.tar.gz && \
    tar xvfz LHAPDF-${LHAPDF_VERSION}.tar.gz && \
    cd LHAPDF-${LHAPDF_VERSION} && \
    ./configure --help && \
    export CXX=$(which g++) && \
    export PYTHON=$(which python) && \
    ./configure \
      --prefix=/usr/local && \
    make -j$(($(nproc) - 1)) && \
    make install && \
    rm -rf /code

# Install MadGraph5_aMC@NLO for Python 3
ARG MG_VERSION=2.7.3
RUN cd /usr/local && \
    wget -q https://launchpad.net/mg5amcnlo/2.0/2.7.x/+download/MG5_aMC_v${MG_VERSION}.py3.tar.gz && \
    tar xzf MG5_aMC_v${MG_VERSION}.py3.tar.gz && \
    rm MG5_aMC_v${MG_VERSION}.py3.tar.gz


# Enable tab completion by uncommenting it from /etc/bash.bashrc
# The relevant lines are those below the phrase "enable bash completion in interactive shells"
RUN export SED_RANGE="$(($(sed -n '\|enable bash completion in interactive shells|=' /etc/bash.bashrc)+1)),$(($(sed -n '\|enable bash completion in interactive shells|=' /etc/bash.bashrc)+7))" && \
    sed -i -e "${SED_RANGE}"' s/^#//' /etc/bash.bashrc && \
    unset SED_RANGE

# Use C.UTF-8 locale to avoid issues with ASCII encoding
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

ENV HOME /home/docker
WORKDIR ${HOME}/data
RUN cp /root/.profile ${HOME}/.profile && \
    cp /root/.bashrc ${HOME}/.bashrc && \
    echo "" >> ${HOME}/.bashrc && \
    echo 'export PATH=${HOME}/.local/bin:$PATH' >> ${HOME}/.bashrc && \
    echo 'export PATH=/usr/local/MG5_aMC_v2_7_3_py3/bin:$PATH' >> ${HOME}/.bashrc && \
    echo 'alias ls="ls -ltrh --color"' >> ${HOME}/.bashrc && \
    python -m pip install --upgrade --no-cache-dir pip setuptools wheel && \
    python -m pip install --no-cache-dir six numpy

#ENV USER docker
#USER docker
ENV PYTHONPATH=/usr/local/lib:$PYTHONPATH
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
ENV PATH ${HOME}/.local/bin:$PATH
ENV PATH /usr/local/MG5_aMC_v2_7_3_py3/bin:$PATH

RUN echo "install pythia8" | mg5_aMC && \
    echo "install mg5amc_py8_interface" | mg5_aMC

ENTRYPOINT ["/bin/bash", "-l", "-c"]
CMD ["/bin/bash"]
