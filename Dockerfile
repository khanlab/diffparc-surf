# Stage: requirements
# Notes: g++ (snakebids), tcsh (freesurfer), libtiff5 (mrtrix3)
FROM python:3.8-slim-bullseye as requirements
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
    curl \
    git \
    g++ \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Stage: build
# Python dependencies
FROM requirements as build 
COPY . /opt/diffparc-surf
RUN cd /opt/diffparc-surf \
    && pip install --prefer-binary --no-cache-dir \
    poetry==1.4.0 \
    && poetry build -f wheel

# Stage: ants
FROM requirements as ants
ARG ANTS_VER=2.4.3
RUN wget https://github.com/ANTsX/ANTs/releases/download/v${ANTS_VER}/ants-${ANTS_VER}-centos7-X64-gcc.zip -O ants.zip \
    && unzip -qq ants.zip -d /opt \
    && mv /opt/ants-${ANTS_VER} /opt/ants \
    && rm ants.zip

# Stage: c3d
FROM requirements as c3d
ARG c3d_VER=1.0.0 
RUN wget https://sourceforge.net/projects/c3d/files/c3d/${c3d_VER}/c3d-${c3d_VER}-Linux-x86_64.tar.gz/download -O c3d.tar.gz \
    && tar -xf c3d.tar.gz -C /opt \
    && mv /opt/c3d-${c3d_VER}-Linux-x86_64 /opt/c3d \ 
    && rm /opt/c3d/bin/c2d /opt/c3d/bin/c3d_gui \
    && rm c3d.tar.gz 

# Stage: greedy
# Experimental builds
FROM requirements as greedy
ARG GREEDY_VER=1.2.0
RUN wget https://sourceforge.net/projects/greedy-reg/files/Experimental/greedy-${GREEDY_VER}-Linux-gcc64.tar.gz/download -O greedy.tar.gz \
    && tar -xf greedy.tar.gz -C /opt \ 
    && mv /opt/greedy-${GREEDY_VER}-Linux-gcc64 /opt/greedy \
    && rm greedy.tar.gz 

# Stage: mrtrix (nogui)
FROM requirements as mrtrix
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
FROM requirements as niftyreg
ARG NIFTYREG_VER=1.3.9
RUN wget https://sourceforge.net/projects/niftyreg/files/nifty_reg-${NIFTYREG_VER}/NiftyReg-${NIFTYREG_VER}-Linux-x86_64-Release.tar.gz/download -O niftyreg.tar.gz \
    && tar -xf niftyreg.tar.gz -C /opt \
    && mv /opt/NiftyReg-${NIFTYREG_VER}-Linux-x86_64-Release /opt/niftyreg \
    && cd /opt/niftyreg/bin \ 
    && ls . | grep -xv "reg_aladin" | xargs rm \
    && cd / \
    && rm niftyreg.tar.gz

# # Stage: synthstrip
FROM freesurfer/synthstrip:1.3 as synthstrip

# Stage: workbench
FROM requirements as workbench
ARG WB_VER=1.5.0
RUN wget https://humanconnectome.org/storage/app/media/workbench/workbench-linux64-v${WB_VER}.zip -O workbench.zip \
    && unzip -qq workbench.zip -d /opt \
    && rm -r /opt/workbench/plugins_linux64 \ 
    && rm workbench.zip 

# Stage: runtime
FROM requirements as runtime
COPY --from=build /opt/diffparc-surf/dist/*.whl /opt/diffparc-surf/
COPY --from=ants \ 
    # Commands to copy
    /opt/ants/bin/antsApplyTransforms \
    /opt/ants/bin/antsRegistration \ 
    /opt/ants/bin/Atropos \ 
    /opt/ants/bin/N4BiasFieldCorrection \ 
    # Target destination
    /opt/ants/bin/
COPY --from=c3d /opt/c3d /opt/c3d/
COPY --from=greedy /opt/greedy/bin/greedy /opt/greedy/bin/
COPY --from=mrtrix /opt/mrtrix3/mrtrix3_runtime /opt/mrtrix3/
COPY --from=niftyreg /opt/niftyreg /opt/niftyreg/
COPY --from=synthstrip /freesurfer /opt/freesurfer/
COPY --from=workbench /opt/workbench /opt/workbench/
RUN WHEEL=`ls /opt/diffparc-surf | grep whl` \
    && pip install /opt/diffparc-surf/${WHEEL} \
    && rm -r /opt/diffparc-surf \
    && apt-get purge -y -q curl g++ unzip wget \
    && apt-get --purge -y -qq autoremove
# Setup environments
ENV OS=Linux \
    ANTSPATH=/opt/ants/bin \
    FREESURFER_HOME=/freesurfer \ 
    LD_LIBRARY_PATH=/opt/niftyreg/lib:/opt/workbench/libs_linux64:/opt/workbench/libs_linux64_software_opengl:${LD_LIBRARY_PATH} \
    PATH=/opt/ants:/opt/ants/bin:/opt/c3d/bin:/opt/freesurfer:/opt/greedy/bin:/opt/mrtrix3/bin:/opt/niftyreg/bin/:/opt/workbench/bin_linux64:${PATH}

# FROM snakemake/snakemake:stable

# COPY . /src
# RUN pip install snakebids
# ENTRYPOINT [ "/src/diffparc/run.py" ]
