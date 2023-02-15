rule nii2mif:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bval=bids(
            root=root,
            suffix="dwi.bval",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bvec=bids(
            root=root,
            suffix="dwi.bvec",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        mask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        dwi=temp(
            bids(
                root=root,
                datatype="dwi",
                suffix="dwi.mif",
                **subj_wildcards,
            )
        ),
        mask=temp(
            bids(
                root=root,
                datatype="dwi",
                suffix="mask.mif",
                **subj_wildcards,
            )
        ),
    threads: 4
    resources:
        mem_mb=16000,
    group:
        "subj"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrconvert {input.dwi} {output.dwi} -fslgrad {input.bvec} {input.bval} -nthreads {threads} && "
        "mrconvert {input.mask} {output.mask} -nthreads {threads}"


rule dseg_nii2mif:
    input:
        "{file}_dseg.nii.gz",
    output:
        temp("{file}_dseg.mif"),
    container:
        config["singularity"]["mrtrix"]
    group:
        "subj"
    shell:
        "mrconvert {input} {output} -nthreads {threads}"


rule dwi2response_msmt:
    # Dhollander, T.; Mito, R.; Raffelt, D. & Connelly, A. Improved white matter response function estimation for 3-tissue constrained spherical deconvolution. Proc Intl Soc Mag Reson Med, 2019, 555
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
    output:
        wm_rf=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="wm",
            suffix="response.txt",
            **subj_wildcards,
        ),
        gm_rf=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="gm",
            suffix="response.txt",
            **subj_wildcards,
        ),
        csf_rf=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="csf",
            suffix="response.txt",
            **subj_wildcards,
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2response dhollander {input.dwi} {output.wm_rf} {output.gm_rf} {output.csf_rf}  -nthreads {threads} -mask {input.mask}"


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
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="wm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
        gm_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="gm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
        csf_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="csf",
            suffix="fod.mif",
            **subj_wildcards,
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj"
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
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="wmnorm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
        gm_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="normalized",
            suffix="gm_fod.mif",
            **subj_wildcards,
        ),
        csf_fod=bids(
            root=root,
            datatype="dwi",
            alg="msmt",
            desc="normalized",
            suffix="csf_fod.mif",
            **subj_wildcards,
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mtnormalise -nthreads {threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod}"


rule dwi2response_csd:
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
    output:
        wm_rf=bids(
            root=root,
            datatype="dwi",
            alg="csd",
            desc="wm",
            suffix="response.txt",
            **subj_wildcards,
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2response fa {input.dwi} {output.wm_rf} -nthreads {threads} -mask {input.mask}"


rule dwi2fod_csd:
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        wm_rf=rules.dwi2response_csd.output.wm_rf,
    output:
        wm_fod=bids(
            root=root,
            datatype="dwi",
            alg="csd",
            desc="wm",
            suffix="fod.mif",
            **subj_wildcards,
        ),
    threads: 8
    resources:
        mem_mb=32000,
    group:
        "subj"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2fod -nthreads {threads} -mask {input.mask} csd {input.dwi} {input.wm_rf} {output.wm_fod}  "


rule dwi2tensor:
    input:
        rules.nii2mif.output.dwi,
    output:
        tensor=bids(
            root=root,
            datatype="dwi",
            suffix="tensor.mif",
            **subj_wildcards,
        ),
    group:
        "subj"
    threads: 8
    resources:
        mem_mb=32000,
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2tensor {input} {output}"


rule dwi_to_tensor:
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
    output:
        tensor=bids(
            root=root,
            datatype="dwi",
            suffix="tensor.mif",
            **subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
    group:
        "subj"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2tensor -mask {input.mask} {input.dwi} {output.tensor} -nthreads {threads}"


rule tensor_to_metrics:
    input:
        tensor=rules.dwi_to_tensor.output.tensor,
        mask=rules.nii2mif.output.mask,
    output:
        fa=bids(
            root=root,
            datatype="dwi",
            suffix="FA.nii.gz",
            **subj_wildcards,
        ),
        md=bids(
            root=root,
            datatype="dwi",
            suffix="MD.nii.gz",
            **subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
    group:
        "subj"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tensor2metric -mask {input.mask} -fa {output.fa} -adc {output.md} {input.tensor}"


def get_fod_for_tracking(wildcards):
    if config["fod_algorithm"] == "csd":
        return (
            bids(
                root=root,
                datatype="dwi",
                alg="csd",
                desc="wm",
                suffix="fod.mif",
                **subj_wildcards,
            ),
        )
    elif config["fod_algorithm"] == "msmt_csd":
        return (
            bids(
                root=root,
                datatype="dwi",
                desc="wmnorm",
                alg="msmt",
                suffix="fod.mif",
                **subj_wildcards,
            ),
        )


rule resample_metric_to_aux_dseg:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            desc="{dseg_method}",
            suffix="dseg.nii.gz",
            **subj_wildcards
        ),
        metric=bids(
            root=root,
            datatype="dwi",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
    output:
        metric=bids(
            root=root,
            datatype="dwi",
            resliced="{dseg_method}",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.dseg} {input.metric} -reslice-identity -o {output.metric}"
