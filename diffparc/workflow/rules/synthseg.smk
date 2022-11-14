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
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="synthseg",
            suffix="dseg.nii.gz"
        ),
    container:
        config["singularity"]["synthseg"]
    threads: 8
    group:
        "subj"
    shell:
        # currently only have models for v1 in the container (onedrive link wasn't working!)
        "python /SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --parc --cpu --threads {threads}"


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
        probseg=bids(
            root=root,
            hemi="{hemi}",
            space="individual",
            label="{seed}",
            from_="synthseg",
            datatype="anat",
            suffix="probseg.nii.gz",
            **subj_wildcards,
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
            space="individual",
            label="{seed}",
            from_="synthseg",
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
            datatype="anat",
            suffix="warp.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="shapeinject",
            from_="{template}",
            to="subject",
            **subj_wildcards
        ),
        invwarp=bids(
            root=root,
            datatype="anat",
            suffix="invwarp.nii.gz",
            desc="shapeinject",
            hemi="{hemi}",
            label="{seed}",
            from_="{template}",
            to="subject",
            **subj_wildcards
        ),
        warped_flo=bids(
            root=root,
            datatype="anat",
            suffix="probseg.nii.gz",
            space="individual",
            hemi="{hemi}",
            label="{seed}",
            from_="{template}",
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
            template="{template}",
            **subj_wildcards
        ),
    shadow:
        "minimal"
    shell:
        "greedy -d 3 -threads {threads} -a -moments 2 -det 1 -m SSD {params.input_fixed_moving} -o affine.txt -n {params.affine_iterations} &> {log} && "
        "greedy -d 3 -threads {threads} -it affine.txt -m SSD  {params.input_fixed_moving} -o {output.warp} -oinv {output.invwarp} -n {params.fluid_iterations} -s {params.gradient_sigma} {params.warp_sigma} -e {params.timestep} &>> {log} && "
        "greedy -d 3 -threads {threads} -rf {input.ref} {params.input_moving_warped} -r {output.warp} affine.txt  &>> {log}"


def get_cmd_synthseg_targets(wildcards, input, output):
    cmd = [f"c3d {input.dseg} -popas IN"]
    for target in config["synthseg_targets"].keys():
        # for each target, we push the input image, retain labels, binarize, rescale, then keep that on the stack for accumulation
        in_labels = " ".join([str(i) for i in config["synthseg_targets"][target]["in"]])
        out_label = config["synthseg_targets"][target]["out"]

        cmd.append(f"-push IN -retain-labels {in_labels} -binarize -scale {out_label}")

    cmd.append(f"-accum -add -endaccum -o {output.dseg}")

    return " ".join(cmd)


rule synthseg_to_targets:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="synthseg",
            suffix="dseg.nii.gz"
        ),
    params:
        cmd=get_cmd_synthseg_targets,
    output:
        dseg=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            desc="{targets}",
            from_="synthseg",
            datatype="anat",
            suffix="dseg.nii.gz"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "{params.cmd}"


rule warp_template_target_dseg:
    input:
        dseg=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["targets"][wildcards.targets]["template_dseg"],
        ),
        ref=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        xfm=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="xfm.txt",
            from_="subject",
            to="{template}",
            desc="affine",
            type_="itk"
        ),
    output:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="dseg.nii.gz",
            from_="{template}",
            targets="{targets}",
            reg="affine",
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.dseg} -o {output.dseg} -r {input.ref} "
        " -t [{input.xfm},1] "


# to map template labels to synthseg cortex, we can
# get the cortical GM map from synthseg (>=1000), then
# relabel according to nearest template label

# to do this efficiently, we can map the background voxels from
# the template label to the nearest foreground voxel first,
# then, just mask that..


rule nearest_label_synthseg_targets:
    input:
        synthseg_dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="synthseg",
            suffix="dseg.nii.gz"
        ),
        template_dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="dseg.nii.gz",
            from_=config["template"],
            targets="{targets}",
            reg="affine",
        ),
    output:
        dseg=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            desc="{targets}",
            from_="synthsegnearest",
            datatype="anat",
            suffix="dseg.nii.gz"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input.template_dseg} -replace 0 inf -split -foreach  -sdt -scale -1 -endfor -merge -popas LBL "
        " {input.synthseg_dseg} -threshold 1000 inf 1 0 -push LBL -multiply -o {output.dseg}"
