# Stage: requirements
# Notes: g++ (snakebids), tcsh (freesurfer), libtiff5 (mrtrix3)
FROM python:3.9-slim-bullseye as requirements
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
    && cp -R bin core lib src mrtrix3_runtime 

# Stage: ants
FROM requirements as ants
ARG ANTS_VER=2.4.3
RUN wget https://github.com/ANTsX/ANTs/releases/download/v${ANTS_VER}/ants-${ANTS_VER}-centos7-X64-gcc.zip -O ants.zip \
    && unzip -qq ants.zip -d /opt \
    && mv /opt/ants-${ANTS_VER} /opt/ants \
    && rm ants.zip

# Stage: runtime
FROM requirements as runtime
COPY --from=ants \ 
    # Commands to copy
    /opt/ants/bin/antsApplyTransforms \
    /opt/ants/bin/antsRegistration \ 
    /opt/ants/bin/Atropos \ 
    /opt/ants/bin/N4BiasFieldCorrection \ 
    # Target destination
    /opt/ants/bin/
COPY --from=mrtrix /opt/mrtrix3/mrtrix3_runtime /opt/mrtrix3/
# Setup environments
ENV OS=Linux \
    ANTSPATH=/opt/ants/bin \
    PATH=/opt/mrtrix3/bin:/opt/ants:/opt/ants/bin:$PATH

# FROM snakemake/snakemake:stable

# COPY . /src
# RUN pip install snakebids
# ENTRYPOINT [ "/src/diffparc/run.py" ]
