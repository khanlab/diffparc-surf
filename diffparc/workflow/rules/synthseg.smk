rule run_synthseg:
    input:
        t1=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
    output:
        dseg=temp(
            bids(
                root=root,
                datatype="anat",
                **subj_wildcards,
                desc="synthsegnoresample",
                suffix="dseg.nii.gz"
            )
        ),
    container:
        config["singularity"]["synthseg"]
    threads: 8
    group:
        "subj"
    shell:
        "python /SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --cpu --threads {threads}"


rule run_synthseg_withcortparc:
    input:
        t1=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
    output:
        dseg=temp(
            bids(
                root=root,
                datatype="anat",
                **subj_wildcards,
                desc="synthsegcortparcnoresample",
                suffix="dseg.nii.gz"
            )
        ),
    container:
        config["singularity"]["synthseg"]
    threads: 8
    group:
        "subj"
    shell:
        "python /SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --cpu --threads {threads} --parc"


rule reslice_synthseg_to_t1:
    input:
        t1=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="{dseg}noresample",
            suffix="dseg.nii.gz"
        ),
    output:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="{dseg,synthseg|synthsegcortparc}",
            suffix="dseg.nii.gz"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d -interpolation NearestNeighbor {input.ref} {input.dseg} -reslice-identity -o {output.dseg}"


rule run_synthseg_template:
    input:
        t1=os.path.join(workflow.basedir, "..", config["template_t1w"]),
    output:
        dseg=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_desc-synthseg_dseg.nii.gz",
    container:
        config["singularity"]["synthseg"]
    threads: 8
    group:
        "subj"
    shell:
        "python /SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --cpu --threads {threads}"


rule extract_synthseg_label:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            desc="synthseg",
            suffix="dseg.nii.gz",
            **subj_wildcards,
        ),
    params:
        labels=lambda wildcards: config["synthseg_labels"][wildcards.seed][
            wildcards.hemi
        ],
        smoothing="1x1x1mm",  #sigma
    output:
        probseg=temp(
            bids(
                root=root,
                hemi="{hemi}",
                label="{seed}",
                desc="synthseg",
                datatype="anat",
                suffix="probseg.nii.gz",
                **subj_wildcards,
            )
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -retain-labels {params.labels} -binarize -smooth {params.smoothing} -o {output}"


rule template_shape_injection:
    input:
        flo=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
        ref=bids(
            root=root,
            hemi="{hemi}",
            label="{seed}",
            desc="synthseg",
            datatype="anat",
            suffix="probseg.nii.gz",
            **subj_wildcards,
        ),
    params:
        input_fixed_moving=lambda wildcards, input: f"-i {input.ref} {input.flo}",
        input_moving_warped=lambda wildcards, input, output: f"-rm {input.flo} {output.warped_flo}",
        affine_iterations="100x50x10",
        fluid_iterations="100x50x10",  #default 100x50x10
        gradient_sigma="1.732vox",  #default 1.732vox
        warp_sigma="0.707vox",  #default 0.707vox
        timestep="1.0",  #default 1.0
    output:
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="shapeinject",
            from_=config["template"],
            to="subject",
            **subj_wildcards
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            suffix="invwarp.nii.gz",
            desc="shapeinject",
            hemi="{hemi}",
            label="{seed}",
            from_=config["template"],
            to="subject",
            **subj_wildcards
        ),
        warped_flo=bids(
            root=root,
            datatype="anat",
            suffix="probseg.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="shapeinject",
            **subj_wildcards
        ),
    threads: 8
    resources:
        mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time=60,  # 1 hrs
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    log:
        bids(
            root="logs",
            suffix="shapeinject.log",
            hemi="{hemi}",
            label="{seed}",
            **subj_wildcards
        ),
    shadow:
        "minimal"
    shell:
        "greedy -d 3 -threads {threads} -a -moments 2 -det 1 -m SSD {params.input_fixed_moving} -o affine.txt -n {params.affine_iterations} &> {log} && "
        "greedy -d 3 -threads {threads} -it affine.txt -m SSD  {params.input_fixed_moving} -o {output.warp} -oinv {output.invwarp} -n {params.fluid_iterations} -s {params.gradient_sigma} {params.warp_sigma} -e {params.timestep} &>> {log} && "
        "greedy -d 3 -threads {threads} -rf {input.ref} {params.input_moving_warped} -r {output.warp} affine.txt  &>> {log}"


def get_subject_seed_probseg(wildcards):
    if config["seeds"][wildcards.seed]["use_synthseg"]:
        return (
            bids(
                root=root,
                datatype="anat",
                suffix="probseg.nii.gz",
                hemi="{hemi}",
                label="{seed}",
                desc="shapeinject",
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


rule linear_transform_aux_dseg_to_template:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="{dseg_method}",
            suffix="dseg.nii.gz"
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
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="{dseg_method}",
            space=config["template"],
            warp="linear",
            suffix="dseg.nii.gz"
        ),
    container:
        config["singularity"]["ants"]
    threads: 8
    resources:
        mem_mb=8000,
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.dseg}  -o {output.dseg}  -r {input.ref} -t {input.affine_xfm_itk} "
