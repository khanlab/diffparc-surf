rule binarize_trim_subject_seed:
    input:
        seed=get_subject_seed_probseg,  #grabs either shapeinject or atlasreg probseg
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
        resample_res=lambda wildcards: config["seeds"][wildcards.seed]["vol_seed_res"],
    output:
        seed_thr=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    log:
        bids(
            root="logs",
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            suffix="binarizetrimsubjectseed.log"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -threshold 0.5 inf 1 0 -resample-mm {params.resample_res} -trim 0vox -type uchar -o {output} &> {log}"


rule create_voxel_seed_images:
    input:
        seed=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    output:
        voxseeds_dir=temp(
            directory(
                bids(
                    root=config["tmp_dir"],
                    **subj_wildcards,
                    hemi="{hemi}",
                    label="{seed}",
                    datatype="dwi",
                    suffix="voxseeds"
                )
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/create_voxel_seed_images.py"


rule track_from_voxels:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        wm_fod=get_fod_for_tracking,
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **subj_wildcards,
        ),
        mask=bids(
            root=root,
            datatype="dwi",
            suffix="mask.mif",
            **subj_wildcards,
        ),
        vox_seeds_dir=bids(
            root=config["tmp_dir"],
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            datatype="dwi",
            suffix="voxseeds"
        ),
    params:
        seeds_per_voxel="{seedpervox}",
    output:
        tck_dir=temp(
            directory(
                bids(
                    root=config["tmp_dir"],
                    datatype="dwi",
                    hemi="{hemi}",
                    label="{seed}",
                    seedpervox="{seedpervox}",
                    suffix="voxtracts",
                    **subj_wildcards,
                )
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "mkdir -p {output.tck_dir} && "
        "parallel --bar --jobs {threads} "
        "tckgen -quiet -nthreads 0  -algorithm iFOD2 -mask {input.mask} "
        " {input.wm_fod} {output.tck_dir}/vox_{{1}}.tck "
        " -seed_random_per_voxel {input.vox_seeds_dir}/seed_{{1}}.nii {params.seeds_per_voxel} "
        " ::: `ls {input.vox_seeds_dir} | grep -Po '(?<=seed_)[0-9]+'`"


rule connectivity_from_voxels:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        tck_dir=bids(
            root=config["tmp_dir"],
            datatype="dwi",
            hemi="{hemi}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="voxtracts",
            **subj_wildcards,
        ),
        targets=bids(
            root=root,
            **subj_wildcards,
            desc="{targets}",
            datatype="anat",
            suffix="dseg.mif"
        ),
    output:
        conn_dir=temp(
            directory(
                bids(
                    root=config["tmp_dir"],
                    datatype="dwi",
                    hemi="{hemi}",
                    desc="{targets}",
                    label="{seed}",
                    seedpervox="{seedpervox}",
                    suffix="voxconn",
                    **subj_wildcards,
                )
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "mkdir -p {output.conn_dir} && "
        "parallel --eta --jobs {threads} "
        "tck2connectome -nthreads 0 -quiet {input.tck_dir}/vox_{{1}}.tck {input.targets} {output.conn_dir}/conn_{{1}}.csv -vector"
        " ::: `ls {input.tck_dir} | grep -Po '(?<=vox_)[0-9]+'`"


rule gen_conn_csv:
    input:
        conn_dir=bids(
            root=config["tmp_dir"],
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="voxconn",
            **subj_wildcards,
        ),
    params:
        header_line=lambda wildcards: ",".join(
            config["targets"][wildcards.targets]["labels"]
        ),
    output:
        conn_csv=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/gather_csv_files.py"


rule conn_csv_to_image:
    input:
        conn_csv=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
        seed_nii=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    output:
        conn_nii=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/conn_csv_to_image.py"


rule nlin_transform_conn_to_template:
    input:
        conn_nii=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
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
            datatype="dwi",
            hemi="{hemi}",
            space=config["template"],
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
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
            seedpervox="{seedpervox}",
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
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
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
            datatype="dwi",
            hemi="{hemi}",
            space=config["template"],
            warp="linear",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
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
            seedpervox="{seedpervox}",
            suffix="transformconntotemplate.log",
            **subj_wildcards,
        ),
    group:
        "subj"
    shell:
        #using nearestneighbor to avoid bluring with background -- background set as -1
        "antsApplyTransforms -d 3 -e 3  --interpolation NearestNeighbor -i {input.conn_nii}  -o {output.conn_nii}  -r {input.ref} -t {input.affine_xfm_itk} &> {log} "


rule maxprob_conn_native:
    """ generate maxprob connectivity, adding outside striatum, and inside striatum (at a particular "streamline count" threshold) to identify unlabelled regions """
    input:
        conn_nii=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    output:
        conn_nii=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
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
            datatype="dwi",
            hemi="{hemi}",
            space=config["template"],
            warp="linear",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    output:
        conn_nii=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            space=config["template"],
            warp="linear",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
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
