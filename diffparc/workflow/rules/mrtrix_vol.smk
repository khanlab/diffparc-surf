rule binarize_trim_subject_seed:
    input:
        seed=get_subject_seed_probseg,  #grabs either shapeinject or atlasreg probseg
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
        resample_res=lambda wildcards: config["resample_seed_res"],
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
        "c3d {input} -threshold {params.threshold} inf 1 0 -resample-mm {params.resample_res} -trim 0vox -type uchar -o {output} &> {log}"


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
                    datatype="tracts",
                    method="mrtrix",
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
            datatype="tracts",
            method="mrtrix",
            suffix="voxseeds"
        ),
    params:
        seedspervoxel="{seedspervoxel}",
        mrtrix_rng_seed=config["mrtrix_rng_seed"],
    output:
        tck_dir=temp(
            directory(
                bids(
                    root=config["tmp_dir"],
                    datatype="tracts",
                    hemi="{hemi}",
                    label="{seed}",
                    seedspervoxel="{seedspervoxel}",
                    method="mrtrix",
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
    benchmark:
        bids(
            root="benchmarks",
            hemi="{hemi}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="mrtrix",
            suffix="voltrack.tsv",
            **subj_wildcards,
        )
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "mkdir -p {output.tck_dir} && "
        "parallel --eta --jobs {threads} "
        "MRTRIX_RNG_SEED={params.mrtrix_rng_seed} tckgen -quiet -nthreads 0  -algorithm iFOD2 -mask {input.mask} "
        " {input.wm_fod} {output.tck_dir}/vox_{{1}}.tck "
        " -seed_random_per_voxel {input.vox_seeds_dir}/seed_{{1}}.nii {params.seedspervoxel} "
        " ::: `ls {input.vox_seeds_dir} | grep -Po '(?<=seed_)[0-9]+'`"


rule connectivity_from_voxels:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        tck_dir=bids(
            root=config["tmp_dir"],
            datatype="tracts",
            hemi="{hemi}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="mrtrix",
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
                    datatype="tracts",
                    hemi="{hemi}",
                    desc="{targets}",
                    label="{seed}",
                    seedspervoxel="{seedspervoxel}",
                    method="mrtrix",
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
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="mrtrix",
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
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="mrtrix",
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
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="mrtrix",
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
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="mrtrix",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/conn_csv_to_image.py"
