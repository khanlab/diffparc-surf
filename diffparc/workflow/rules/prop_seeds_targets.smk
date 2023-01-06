rule create_upsampled_cropped_seed_ref:
    input:
        ref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    params:
        resample_res=config["resample_seed_res"],
    output:
        ref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -resample-mm {params.resample_res}  -o {output}"


rule transform_seed_to_subject:
    input:
        seed=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
        ref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **subj_wildcards
        ),
        inv_warp=bids(
            root=root,
            datatype="warps",
            suffix="invwarp.nii.gz",
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
    output:
        seed=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            datatype="anat",
            suffix="probseg.nii.gz"
        ),
    envmodules:
        "ants",
    container:
        config["singularity"]["ants"]
    log:
        bids(
            root="logs",
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            suffix="transformseedtosubject.log"
        ),
    group:
        "subj"
    threads: 8
    shell:
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.seed} -o {output} -r {input.ref} -t [{input.affine_xfm_itk},1] {input.inv_warp}  &> {log}"


rule transform_targets_to_subject:
    input:
        targets=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["targets"][wildcards.targets]["template_dseg"],
        ),
        ref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **subj_wildcards
        ),
        inv_warp=bids(
            root=root,
            datatype="warps",
            suffix="invwarp.nii.gz",
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
    output:
        targets=bids(
            root=root,
            **subj_wildcards,
            desc="{targets}",
            datatype="anat",
            suffix="dseg.nii.gz"
        ),
    envmodules:
        "ants",
    container:
        config["singularity"]["ants"]
    log:
        bids(
            root="logs",
            **subj_wildcards,
            desc="{targets}",
            suffix="transformtargetstosubjects.log"
        ),
    group:
        "subj"
    threads: 8
    shell:
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.targets} -o {output} -r {input.ref} -t [{input.affine_xfm_itk},1] {input.inv_warp}  &> {log}"


rule reslice_shapeinject_to_ref:
    input:
        probseg=bids(
            root=root,
            datatype="anat",
            suffix="probseg.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="shapeinject",
            **subj_wildcards
        ),
        ref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        probseg=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="probseg.nii.gz",
                hemi="{hemi}",
                label="{seed}",
                desc="shapeinjectupsampled",
                **subj_wildcards
            )
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input.ref} {input.probseg} -reslice-identity -o {output.probseg}"


def get_subject_seed_upsampled_probseg(wildcards):
    if config["seeds"][wildcards.seed]["use_synthseg"]:
        return (
            bids(
                root=root,
                datatype="anat",
                suffix="probseg.nii.gz",
                hemi="{hemi}",
                label="{seed}",
                desc="shapeinjectupsampled",
                **subj_wildcards
            ),
        )
    else:
        return (
            bids(
                root=root,
                **subj_wildcards,
                hemi="{hemi}",
                label="{seed}",
                datatype="anat",
                suffix="probseg.nii.gz"
            ),
        )
