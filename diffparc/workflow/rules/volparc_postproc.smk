
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
        config["singularity"]["ants"]
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
        config["singularity"]["ants"]
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
        config["singularity"]["itksnap"]
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
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c4d {input} -slice w 0:-1 -vote -o {output} "
