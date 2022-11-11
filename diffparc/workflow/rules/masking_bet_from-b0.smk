# input: bids(root='work',suffix='b0.nii.gz',desc='dwiref',datatype='dwi',**config['subj_wildcards'])
# output: bids(root='work',desc='brain',method='bet_from-b0',suffix='mask.nii.gz',datatype='dwi',**config['subj_wildcards'])


rule import_avg_b0:
    input:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="dwiref",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="dwiref",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


# n4
rule n4_avg_b0:
    input:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="dwiref",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="n4",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "N4BiasFieldCorrection -i {input} -o {output}"


# rescale intensities, clip off first/last 5% of intensities, then rescale to 0-2000
rule rescale_avg_b0:
    input:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="n4",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="rescale",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d -verbose {input} -clip 5% 95% -stretch 0% 99% 0 2000 -o {output}"


rule bet_avg_b0:
    input:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="rescale",
            datatype="dwi",
            **subj_wildcards
        ),
    params:
        bet_frac=config["b0_bet_frac"],
    output:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="bet",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["prepdwi"]  #fsl
    group:
        "subj"
    shell:
        "bet {input} {output} -f {params.bet_frac}"


rule binarize_avg_b0:
    input:
        bids(
            root="work",
            suffix="b0.nii.gz",
            desc="bet",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        bids(
            root="work",
            suffix="mask.nii.gz",
            desc="brain",
            method="bet_from-b0",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -binarize  -o {output}"
