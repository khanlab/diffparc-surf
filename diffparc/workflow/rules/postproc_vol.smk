
rule nlin_transform_conn_to_template:
    input:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            from_="subject",
            to=config["template"],
            **subj_wildcards
        ),
        affine_xfm_itk=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="itk",
            **subj_wildcards
        ),
        ref=os.path.join(workflow.basedir, "..", config["template_t1w"]),
    output:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            space=config["template"],
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    threads: 8
    resources:
        mem_mb=8000,
    log:
        bids(
            root="logs",
            hemi="{hemi}",
            space=config["template"],
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="transformconntotemplate.log",
            **subj_wildcards,
        ),
    group:
        "subj"
    shell:
        #using nearestneighbor to avoid bluring with background -- background set as -1
        "antsApplyTransforms -d 3 -e 3  --interpolation NearestNeighbor -i {input.conn_nii}  -o {output.conn_nii}  -r {input.ref} -t {input.warp} -t {input.affine_xfm_itk} &> {log} "


rule linear_transform_conn_to_template:
    input:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
        affine_xfm_itk=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="itk",
            **subj_wildcards
        ),
        ref=os.path.join(workflow.basedir, "..", config["template_t1w"]),
    output:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            space=config["template"],
            warp="linear",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    threads: 8
    resources:
        mem_mb=8000,
    log:
        bids(
            root="logs",
            hemi="{hemi}",
            space=config["template"],
            warp="linear",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="transformconntotemplate.log",
            **subj_wildcards,
        ),
    group:
        "subj"
    shell:
        #using nearestneighbor to avoid bluring with background -- background set as -1
        "antsApplyTransforms -d 3 -e 3  --interpolation NearestNeighbor -i {input.conn_nii}  -o {output.conn_nii}  -r {input.ref} -t {input.affine_xfm_itk} &> {log} "


rule maxprob_conn_native:
    input:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    output:
        maxprob_nii=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            segtype="maxprob",
            suffix="dseg.nii.gz",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c4d {input} -slice w 0:-1 -vote -o {output} "


rule maxprob_conn_linMNI:
    """ generate maxprob connectivity, adding outside striatum, and inside striatum (at a particular "streamline count" threshold) to identify unlabelled regions """
    input:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            space=config["template"],
            warp="linear",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    output:
        conn_nii=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            space=config["template"],
            warp="linear",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            segtype="maxprob",
            suffix="dseg.nii.gz",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c4d {input} -slice w 0:-1 -vote -o {output} "


rule binarize_trim_template_seed:
    input:
        seed=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
        resample_res=lambda wildcards: config["resample_seed_res"],
    output:
        seed_thr=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_hemi-{hemi}_label-{seed}_desc-cropped_mask.nii.gz",
    log:
        bids(
            root="logs",
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            suffix="binarizetrimtemplateseed.log"
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d {input} -threshold {params.threshold} inf 1 0 -resample-mm {params.resample_res} -trim 0vox -type uchar -o {output} &> {log}"


rule nlin_transform_conn_to_crop_template:
    """ this uses the template seed tightly cropped for efficiency 
    when working with a large number of targets"""
    input:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            from_="subject",
            to=config["template"],
            **subj_wildcards
        ),
        affine_xfm_itk=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="itk",
            **subj_wildcards
        ),
        ref=rules.binarize_trim_template_seed.output.seed_thr,
    output:
        conn_nii=temp(
            bids(
                root=root,
                datatype="tracts",
                hemi="{hemi}",
                warp="nlincrop",
                space=config["template"],
                desc="{targets}",
                label="{seed}",
                seedspervoxel="{seedspervoxel}",
                method="{method}",
                suffix="unmaskedconn.nii.gz",
                **subj_wildcards,
            )
        ),
    container:
        config["singularity"]["diffparc"]
    threads: 8
    resources:
        mem_mb=8000,
    log:
        bids(
            root="logs",
            hemi="{hemi}",
            space=config["template"],
            warp="nlincrop",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="transformconntotemplate.log",
            **subj_wildcards,
        ),
    group:
        "subj"
    shell:
        #linear interpolation since we have dilated the seed mask (so no risk of blending with background)
        "antsApplyTransforms -d 3 -e 3  --interpolation Linear -i {input.conn_nii}  -o {output.conn_nii}  -r {input.ref} -t {input.warp} -t {input.affine_xfm_itk} &> {log}"


rule multiply_nlin_conn_with_template_mask:
    input:
        data=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            warp="nlincrop",
            space=config["template"],
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="unmaskedconn.nii.gz",
            **subj_wildcards,
        ),
        mask=rules.binarize_trim_template_seed.output.seed_thr,
    output:
        data=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            warp="nlincrop",
            space=config["template"],
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/mask_4d_with_3d.py"
