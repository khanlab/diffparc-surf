# rules for normalizing dti metrics to help deal with site effects


rule zscore_norm:
    input:
        metric=bids(
            root=root, datatype="dwi", suffix="{metric}.nii.gz", **subj_wildcards
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
        metric=bids(
            root=root,
            datatype="dwi",
            desc="zscore",
            suffix="norm{metric}.nii.gz",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/normalize_by_zscore.py"
