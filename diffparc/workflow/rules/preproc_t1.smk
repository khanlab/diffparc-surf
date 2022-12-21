rule import_t1:
    """ This currently just grabs the first T1w image """
    input:
        lambda wildcards: expand(
            input_path["T1w"],
            zip,
            **snakebids.filter_list(input_zip_lists["T1w"], wildcards)
        )[0],
    output:
        bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule synthstrip_t1:
    input:
        t1=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
    output:
        mask=temp(
            bids(
                root=root,
                datatype="anat",
                **subj_wildcards,
                desc="nofixhdrbrain",
                suffix="mask.nii.gz"
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["synthstrip"]
    threads: 8
    shell:
        "python3 /freesurfer/mri_synthstrip -i {input.t1} -m {output.mask} --no-csf"


rule fixheader_synthstrip:
    input:
        t1=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="nofixhdrbrain",
            suffix="mask.nii.gz"
        ),
    output:
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="brain",
            suffix="mask.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.t1} {input.mask} -copy-transform -o {output.mask}"


rule n4_t1_withmask:
    input:
        t1=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="brain",
            suffix="mask.nii.gz"
        ),
    output:
        t1=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
    threads: 8
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -d 3 -i {input.t1} -x {input.mask} -o {output}"


rule mask_subject_t1w:
    input:
        t1=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="brain",
            suffix="mask.nii.gz"
        ),
    output:
        t1=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="T1w.nii.gz",
            desc="masked"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -multiply -o {output}"
