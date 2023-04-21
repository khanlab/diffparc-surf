# Stage: python
FROM nvidia/cuda:10.1-cudnn7-runtime-ubuntu18.04 as python
WORKDIR /root
ENV DEBIAN_FRONTEND="noninteractive" \
    PYENV_VER="v2.3.17" \
    PYTHON_VER="3.8.12" \
    PYENV_ROOT="/root/.pyenv" \
    PATH="/root/.pyenv/shims:/root/.pyenv/bin:${PATH}"
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
    build-essential \
    curl \
    git \
    libbz2-dev \
    libffi-dev \ 
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \ 
    libssl-dev \
    llvm \
    python-openssl \
    tk-dev \
    wget \
    xz-utils \
    zlib1g-dev \
    software-properties-common \
    && git clone -b ${PYENV_VER} --depth=1 https://github.com/pyenv/pyenv.git .pyenv \
    && pyenv install -v ${PYTHON_VER} \
    && pyenv global ${PYTHON_VER} \
    && python -m pip install --upgrade pip

# Stage: builder
# Notes: g++ (snakebids), tcsh (freesurfer), libtiff5 (mrtrix3)
FROM python as builder
RUN mkdir -p /opt \
    && apt-get update -qq \ 
    && apt-get install -y -q --no-install-recommends \
    libeigen3-dev \
    zlib1g-dev \
    libqt5opengl5-dev \
    libqt5svg5-dev \
    libgl1-mesa-dev \
    libncurses5 \
    libxt6 \
    libtiff5 \
    parallel \
    rsync \
    tcsh \
    g++ \
    unzip \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Stage: diffparc (python wheel)
FROM builder as diffparc
COPY . /opt/diffparc-surf
RUN cd /opt/diffparc-surf \
    && pip install --prefer-binary --no-cache-dir \
    poetry==1.4.0 \
    && poetry build -f wheel

# Stage: ants
FROM builder as ants
ARG ANTS_VER=2.4.3
RUN wget https://github.com/ANTsX/ANTs/releases/download/v${ANTS_VER}/ants-${ANTS_VER}-centos7-X64-gcc.zip -O ants.zip \
    && unzip -qq ants.zip -d /opt \
    && mv /opt/ants-${ANTS_VER} /opt/ants \
    && rm ants.zip

# Stage: itksnap (built with Ubuntu16.04 - glibc 2.23)
FROM khanlab/itksnap:main as itksnap
RUN cp -R /opt/itksnap/ /opt/itksnap-mini/ \
    && cd /opt/itksnap-mini/bin \
    && rm c2d itksnap* 

# Stage: mrtrix (nogui)
FROM builder as mrtrix
ARG MRTRIX_VER=3.0.4
ARG MRTRIX_CONFIG_FLAGS="-nogui"
ARG MRTRIX_BUILD_FLAGS="" 
RUN git clone --depth 1 https://github.com/MRtrix3/mrtrix3.git /opt/mrtrix3 \
    && cd /opt/mrtrix3 \
    && git fetch --tags \
    && git checkout ${MRTRIX_VER} \
    && ./configure ${MRTRIX_CONFIG_FLAGS} \
    && ./build ${MRTRIX_BUILD_FLAGS} \
    && mkdir -p mrtrix3_runtime \
    && cp -R bin core lib src mrtrix3_runtime \ 
    && rm -r tmp

# Stage: niftyreg
FROM builder as niftyreg
ARG NIFTYREG_VER=1.3.9
RUN wget https://sourceforge.net/projects/niftyreg/files/nifty_reg-${NIFTYREG_VER}/NiftyReg-${NIFTYREG_VER}-Linux-x86_64-Release.tar.gz/download -O niftyreg.tar.gz \
    && tar -xf niftyreg.tar.gz -C /opt \
    && mv /opt/NiftyReg-${NIFTYREG_VER}-Linux-x86_64-Release /opt/niftyreg \
    && cd /opt/niftyreg/bin \ 
    && ls . | grep -xv "reg_aladin" | xargs rm \
    && rm /root/niftyreg.tar.gz

# Stage: synthstrip
FROM freesurfer/synthstrip:1.3 as synthstrip

# Stage: synthseg
FROM builder as synthseg
RUN git clone https://github.com/BBillot/Synthseg.git /opt/SynthSeg \
    && curl https://www.dropbox.com/s/i62tzl821mqi3vd/synths_models_freesurfer_7.3.2.zip -Lo synthseg_models_freesurfer_7.3.2.zip \
    && unzip -qqo synthseg_models_freesurfer_7.3.2.zip -d /opt/SynthSeg/models \
    && rm synthseg_models_freesurfer_7.3.2.zip 

# Stage: workbench
FROM builder as workbench
ARG WB_VER=1.5.0
RUN wget https://humanconnectome.org/storage/app/media/workbench/workbench-linux64-v${WB_VER}.zip -O workbench.zip \
    && unzip -qq workbench.zip -d /opt \
    && rm -r /opt/workbench/plugins_linux64 \ 
    && rm workbench.zip 

# Stage: runtime
FROM builder as runtime
COPY --from=diffparc /opt/diffparc-surf/dist/*.whl /opt/diffparc-surf/
COPY --from=ants \ 
    # Commands to copy
    /opt/ants/bin/antsApplyTransforms \
    /opt/ants/bin/antsRegistration \ 
    /opt/ants/bin/Atropos \ 
    /opt/ants/bin/N4BiasFieldCorrection \ 
    # Target destination
    /opt/ants/bin/
COPY --from=itksnap /opt/itksnap-mini /opt/itksnap/
COPY --from=mrtrix /opt/mrtrix3/mrtrix3_runtime /opt/mrtrix3/
COPY --from=niftyreg /opt/niftyreg /opt/niftyreg/
COPY --from=synthseg /opt/SynthSeg /opt/SynthSeg/
COPY --from=synthstrip /freesurfer /opt/freesurfer/
COPY --from=workbench /opt/workbench /opt/workbench/
# Setup environments
ENV OS=Linux \
    ANTSPATH=/opt/ants/bin \
    FREESURFER_HOME=/opt/freesurfer \ 
    LD_LIBRARY_PATH=/opt/itksnap/lib:/opt/niftyreg/lib:/opt/workbench/libs_linux64:/opt/workbench/libs_linux64_software_opengl:/usr/local/cuda/lib64:/usr/local/cuda/lib64:${LD_LIBRARY_PATH} \
    PATH=/opt/ants:/opt/ants/bin:/opt/freesurfer:/opt/itksnap/bin:/opt/mrtrix3/bin:/opt/niftyreg/bin/:/opt/workbench/bin_linux64:/usr/local/cuda/bin:${PATH}
RUN WHEEL=`ls /opt/diffparc-surf | grep whl` \
    && pip install /opt/diffparc-surf/${WHEEL} \
    && rm -r /opt/diffparc-surf \
    && apt-get purge -y -q curl g++ unzip wget \
    && apt-get --purge -y -qq autoremove
ENTRYPOINT ["diffparc"]