# obtain VBM-like features using the synthseg labels to obtain tissue density, modulated by the jacobian
wildcard_constraints:
    tissue="[a-zA-Z0-9]+",
    fwhm="[0-9]+",


# using the synthseg tissue labels is one approach, other one would be to use a GMM or MRF-based segmentation
rule synthseg_to_tissue:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="synthseg",
            suffix="dseg.nii.gz"
        ),
        label_tsv=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["aux_dseg"]["synthseg"]["label_tsv"],
        ),
    params:
        tissue_lut=config["vbm"]["tissue_lut"],
    output:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            method="synthseg",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/seg_to_tissue.py"


rule extract_tissue_density:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            method="{method}",
            suffix="dseg.nii.gz"
        ),
    params:
        smoothing_fwhm="1x1x1mm",  #nominal amount of smoothing prior to warping+modulation (will be smoothed afterwards)
        label_num=lambda wildcards: config["vbm"]["tissue_lut"][wildcards.tissue],
    output:
        density=temp(
            bids(
                root=root,
                datatype="vbm",
                **subj_wildcards,
                method="{method}",
                desc="{tissue}",
                suffix="density.nii.gz"
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.dseg} -retain-labels {params.label_num} "
        " -binarize -smooth {params.smoothing_fwhm} "
        " -o {output.density}"


rule compose_warps_to_template:
    input:
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            from_="subject",
            to=config["template"],
            **subj_wildcards
        ),
        affine_xfm_ras=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="ras",
            **subj_wildcards
        ),
        ref=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_desc-masked_T1w.nii.gz",
    output:
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            desc="composeaffine",
            from_="subject",
            to=config["template"],
            **subj_wildcards
        ),
    threads: 8
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "greedy -d 3 -threads {threads} -rf {input.ref} -rc {output.warp} -r {input.warp} {input.affine_xfm_ras}"


rule warp_tissue_density_to_template:
    input:
        density=rules.extract_tissue_density.output.density,
        warp=rules.compose_warps_to_template.output.warp,
        ref=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_desc-masked_T1w.nii.gz",
    output:
        density=bids(
            root=root,
            datatype="vbm",
            **subj_wildcards,
            space=config["template"],
            method="{method}",
            desc="{tissue}",
            suffix="density.nii.gz"
        ),
    threads: 8
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "greedy -d 3 -threads {threads} -rf {input.ref} -rm {input.density} {output.density} -r {input.warp}"


rule calc_warp_det_jac:
    input:
        warp=rules.compose_warps_to_template.output.warp,
    output:
        detjac=bids(
            root=root,
            datatype="warps",
            suffix="detjac.nii.gz",
            desc="composeaffine",
            from_="subject",
            to=config["template"],
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["ants"]
    shell:
        "CreateJacobianDeterminantImage 3 {input.warp} {output.detjac}"


rule modulate_tissue_density:
    input:
        density=rules.warp_tissue_density_to_template.output.density,
        detjac=rules.calc_warp_det_jac.output.detjac,
    output:
        density=bids(
            root=root,
            datatype="vbm",
            **subj_wildcards,
            space=config["template"],
            method="{method}",
            desc="{tissue}",
            suffix="moddensity.nii.gz"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        #note: order of images for divide seems to be different from c3d docs
        "c3d  {input.detjac} {input.density} -divide -o {output.density}"


rule smooth_density_map:
    input:
        density=rules.modulate_tissue_density.output.density,
    params:
        smoothing_fwhm="{fwhm}x{fwhm}x{fwhm}mm",
    output:
        density=bids(
            root=root,
            datatype="vbm",
            **subj_wildcards,
            space=config["template"],
            method="{method}",
            fwhm="{fwhm}mm",
            desc="{tissue}",
            suffix="moddensity.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.density} "
        " -smooth {params.smoothing_fwhm} "
        " -o {output.density}"


rule transform_dti_metric_to_template:
    input:
        metric=bids(
            root=root,
            datatype="dwi",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            desc="composeaffine",
            from_="subject",
            to=config["template"],
            **subj_wildcards
        ),
        ref=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_desc-masked_T1w.nii.gz",
    output:
        metric=bids(
            root=root,
            datatype="vbm",
            space=config["template"],
            desc="dti",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
    threads: 8
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "greedy -d 3 -threads {threads} -rf {input.ref} -rm {input.metric} {output.metric} -r {input.warp}"


rule smooth_dti_metric:
    input:
        metric=bids(
            root=root,
            datatype="vbm",
            space=config["template"],
            desc="dti",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
    params:
        smoothing_fwhm="{fwhm}x{fwhm}x{fwhm}mm",
    output:
        metric=bids(
            root=root,
            datatype="vbm",
            space=config["template"],
            desc="dti",
            fwhm="{fwhm}mm",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input} "
        " -smooth {params.smoothing_fwhm} "
        " -o {output}"
