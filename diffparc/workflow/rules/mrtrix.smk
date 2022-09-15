# ----------- MRTRIX PREPROC BEGIN ------------#
rule nii2mif:
    input:
        dwi=bids(
            root="results",
            suffix="dwi.nii.gz",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        bval=bids(
            root="results",
            suffix="dwi.bval",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        bvec=bids(
            root="results",
            suffix="dwi.bvec",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        mask=bids(
            root="results",
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **config["subj_wildcards"]
        ),
    output:
        dwi=bids(
            root="results",
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
        mask=bids(
            root="results",
            datatype="dwi",
            suffix="mask.mif",
            **config["subj_wildcards"],
        ),
    threads: 4
    resources:
        mem_mb=16000,
    group:
        "subj1"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrconvert {input.dwi} {output.dwi} -fslgrad {input.bvec} {input.bval} -nthreads {threads} && "
        "mrconvert {input.mask} {output.mask} -nthreads {threads}"


rule dwi2response_msmt:
    # Dhollander, T.; Mito, R.; Raffelt, D. & Connelly, A. Improved white matter response function estimation for 3-tissue constrained spherical deconvolution. Proc Intl Soc Mag Reson Med, 2019, 555
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        bvec=bids(
            root="results",
            suffix="dwi.bvec",
            desc="preproc",
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        bval=bids(
            root="results",
            suffix="dwi.bval",
            desc="preproc",
            datatype="dwi",
            **config["subj_wildcards"]
        ),
    output:
        wm_rf=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="wm",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
        gm_rf=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="gm",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
        csf_rf=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="csf",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj1"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2response dhollander {input.dwi} {output.wm_rf} {output.gm_rf} {output.csf_rf} -fslgrad {input.bvec} {input.bval} -nthreads {threads} -mask {input.mask}"


rule dwi2fod_msmt:
    # Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        wm_rf=rules.dwi2response_msmt.output.wm_rf,
        gm_rf=rules.dwi2response_msmt.output.gm_rf,
        csf_rf=rules.dwi2response_msmt.output.csf_rf,
    output:
        wm_fod=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="wm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="gm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="csf",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2fod -nthreads {threads} -mask {input.mask} msmt_csd {input.dwi} {input.wm_rf} {output.wm_fod} {input.gm_rf} {output.gm_fod} {input.csf_rf} {output.csf_fod} "


rule mtnormalise:
    # Raffelt, D.; Dhollander, T.; Tournier, J.-D.; Tabbara, R.; Smith, R. E.; Pierre, E. & Connelly, A. Bias Field Correction and Intensity Normalisation for Quantitative Analysis of Apparent Fibre Density. In Proc. ISMRM, 2017, 26, 3541
    # Dhollander, T.; Tabbara, R.; Rosnarho-Tornstrand, J.; Tournier, J.-D.; Raffelt, D. & Connelly, A. Multi-tissue log-domain intensity and inhomogeneity normalisation for quantitative apparent fibre density. In Proc. ISMRM, 2021, 29, 2472
    input:
        wm_fod=rules.dwi2fod_msmt.output.wm_fod,
        gm_fod=rules.dwi2fod_msmt.output.gm_fod,
        csf_fod=rules.dwi2fod_msmt.output.csf_fod,
        mask=rules.nii2mif.output.mask,
    output:
        wm_fod=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="wmnorm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="normalized",
            suffix="gm_fod.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids(
            root="results",
            datatype="dwi",
            alg="msmt",
            desc="normalized",
            suffix="csf_fod.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mtnormalise -nthreads {threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod}"


rule dwi2response_csd:
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        bvec=bids(
            root="results",
            suffix="dwi.bvec",
            desc="preproc",
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        bval=bids(
            root="results",
            suffix="dwi.bval",
            desc="preproc",
            datatype="dwi",
            **config["subj_wildcards"]
        ),
    output:
        wm_rf=bids(
            root="results",
            datatype="dwi",
            alg="csd",
            desc="wm",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj1"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2response fa {input.dwi} {output.wm_rf}  -fslgrad {input.bvec} {input.bval} -nthreads {threads} -mask {input.mask}"


rule dwi2fod_csd:
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        wm_rf=rules.dwi2response_csd.output.wm_rf,
    output:
        wm_fod=bids(
            root="results",
            datatype="dwi",
            alg="csd",
            desc="wm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2fod -nthreads {threads} -mask {input.mask} csd {input.dwi} {input.wm_rf} {output.wm_fod}  "


rule dwi2tensor:
    input:
        rules.nii2mif.output.dwi,
    output:
        tensor=bids(
            root="results",
            datatype="dwi",
            suffix="tensor.mif",
            **config["subj_wildcards"],
        ),
    group:
        "subj1"
    threads: 8
    resources:
        mem_mb=32000,
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2tensor {input} {output}"


rule tensor2metrics:
    input:
        tensor=rules.dwi2tensor.output.tensor,
        mask=rules.nii2mif.output.mask,
    output:
        fa=bids(
            root="results",
            datatype="dwi",
            suffix="fa.mif",
            **config["subj_wildcards"],
        ),
    group:
        "subj1"
    threads: 8
    resources:
        mem_mb=32000,
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tensor2metric -fa {output.fa} -mask {input.mask} {input.tensor}"


# -------------- MRTRIX PREPROC END ----------------#


# ----------- MRTRIX TRACTOGRAPHY BEGIN ------------#
rule create_seed:
    input:
        rules.tensor2metrics.output.fa,
    params:
        threshold=0.15,
    output:
        seed=bids(
            root="results",
            datatype="dwi",
            suffix="seed.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrthreshold {input} -abs {params.threshold} {output}"


rule tckgen:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        wm_fod=rules.mtnormalise.output.wm_fod,
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        seed=rules.create_seed.output.seed,
    params:
        streamlines=5000,
        seed_strategy=lambda wildcards, input: f"-seed_image {input.seed}",
    output:
        tck=bids(
            root="results",
            datatype="dwi",
            desc="iFOD2",
            suffix="tractography.tck",
            **config["subj_wildcards"],
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckgen -nthreads {threads} -algorithm iFOD2 -mask {input.mask} {params.seed_strategy} -select {params.streamlines} {input.wm_fod} {output.tck}"