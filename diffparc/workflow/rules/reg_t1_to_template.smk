


rule affine_to_template:
    input:
        flo=bids(
            root=root,
            suffix="T1w.nii.gz",
            desc="n4",
            datatype="anat",
            **subj_wildcards
        ),
        ref=os.path.join(workflow.basedir, "..", config["template_t1w"]),
    output:
        warped_subj=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="T1w.nii.gz",
            space="{template}",
            desc="affine"
        ),
        xfm_ras=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="xfm.txt",
            from_="subject",
            to="{template}",
            desc="affine",
            type_="ras"
        ),
    container:
        config["singularity"]["prepdwi"]  #niftyreg
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped_subj} -aff {output.xfm_ras}"


rule convert_template_xfm_ras2itk:
    input:
        bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="xfm.txt",
            from_="subject",
            to="{template}",
            desc="{desc}",
            type_="ras"
        ),
    output:
        bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="xfm.txt",
            from_="subject",
            to="{template}",
            desc="{desc}",
            type_="itk"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d_affine_tool {input}  -oitk {output}"


rule warp_brainmask_from_template_affine:
    input:
        mask=os.path.join(workflow.basedir, "..", config["template_mask"]),
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
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="mask.nii.gz",
            from_="{template}",
            reg="affine",
            desc="brain"
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.mask} -o {output.mask} -r {input.ref} "
        " -t [{input.xfm},1] "


rule warp_tissue_probseg_from_template_affine:
    input:
        probseg=os.path.join(workflow.basedir, "..", config["template_tissue_probseg"]),
        ref=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        xfm=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="xfm.txt",
            from_="subject",
            to="{template}",
            desc="{desc}",
            type_="itk"
        ),
    output:
        probseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="probseg.nii.gz",
            label="{tissue}",
            from_="{template}",
            reg="{desc}"
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.probseg} -o {output.probseg} -r {input.ref} "
        " -t [{input.xfm},1]"


rule n4_t1_withmask:
    input:
        t1=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="mask.nii.gz",
            from_="{template}".format(template=config["template"]),
            reg="affine",
            desc="brain"
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


rule mask_template_t1w:
    input:
        t1=os.path.join(workflow.basedir, "..", config["template_t1w"]),
        mask=os.path.join(workflow.basedir, "..", config["template_mask"]),
    output:
        t1=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template="{template}"
        )
        + "_desc-masked_T1w.nii.gz",
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -multiply -o {output}"


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
            suffix="mask.nii.gz",
            from_="atropos3seg",
            desc="brain"
        ),
    output:
        t1=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="T1w.nii.gz",
            from_="atropos3seg",
            desc="masked"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -multiply -o {output}"


rule greedy_affine_init:
    input:
        flo=[
            bids(
                root=root,
                datatype="anat",
                **subj_wildcards,
                suffix="T1w.nii.gz",
                from_="atropos3seg",
                desc="masked"
            ),
            expand(
                bids(
                    root=root,
                    datatype="anat",
                    **subj_wildcards,
                    suffix="probseg.nii.gz",
                    label="{tissue}",
                    desc="atropos3seg"
                ),
                tissue=config["tissue_labels"],
                allow_missing=True,
            ),
        ],
        ref=[
            get_template_prefix(
                root=root, subj_wildcards=subj_wildcards, template="{template}"
            )
            + "_desc-masked_T1w.nii.gz",
            expand(
                os.path.join(
                    workflow.basedir, "..", config["template_tissue_probseg"]
                ),
                tissue=config["tissue_labels"],
                allow_missing=True,
            ),
        ],
        init_xfm=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="xfm.txt",
            from_="subject",
            to="{template}",
            desc="affine",
            type_="itk"
        ),
    params:
        input_fixed_moving=lambda wildcards, input: [
            f"-i {fixed} {moving}" for fixed, moving in zip(input.ref, input.flo)
        ],
        input_moving_warped=lambda wildcards, input, output: [
            f"-rm {moving} {warped}"
            for moving, warped in zip(input.flo, output.warped_flo)
        ],
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
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
        invwarp=bids(
            root=root,
            datatype="anat",
            suffix="invwarp.nii.gz",
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
        warped_flo=[
            bids(
                root=root,
                datatype="anat",
                suffix="T1w.nii.gz",
                space="{template}",
                desc="greedy",
                **subj_wildcards
            )
        ],
        affine=bids(
            root=root,
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="itk",
            **subj_wildcards
        ),
        affine_xfm_ras=bids(
            root=root,
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="ras",
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
        bids(root="logs", suffix="greedy.log", template="{template}", **subj_wildcards),
    shell:
        #affine first
        "greedy -d 3 -threads {threads} -a -m NCC 2x2x2 {params.input_fixed_moving} -o {output.affine_xfm_ras} -ia-image-centers -n {params.affine_iterations} &> {log} && "

        "greedy -d 3 -threads {threads} -m NCC 2x2x2 {params.input_fixed_moving} -it {output.affine_xfm_ras} -o {output.warp} -oinv {output.invwarp} -n {params.fluid_iterations} -s {params.gradient_sigma} {params.warp_sigma} -e {params.timestep} &>> {log} && "

        "c3d_affine_tool {output.affine_xfm_ras} -oitk {output.affine} &>> {log} && "

        "greedy -d 3 -threads {threads} -rf {input.ref[0]} {params.input_moving_warped} -r {output.warp} {output.affine_xfm_ras} &>> {log}"
        #then deformable:
        #then convert affine to itk format that ants uses
        #and finally warp the moving image


"""
rule ants_syn_affine_init:
    input:
        flo=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="T1w.nii.gz",
            from_="atropos3seg",
            desc="masked"
        ),

        ref=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template="{template}"
        )+"_desc-masked_T1w.nii.gz",

        init_xfm=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="xfm.txt",
            from_="subject",
            to="{template}",
            desc="affine",
            type_="itk"
        ),
    params:
        out_prefix=bids(
            root=root,
            datatype="anat",
            suffix="",
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
        base_opts="--write-composite-transform -d {dim} --float 1 ".format(
            dim=config["ants"]["dim"]
        ),
        intensity_opts=config["ants"]["intensity_opts"],
        init_transform=lambda wildcards, input: "-r {xfm}".format(xfm=input.init_xfm),
        linear_multires="-c [{reg_iterations},1e-6,10] -f {shrink_factors} -s {smoothing_factors}".format(
            reg_iterations=config["ants"]["linear"]["reg_iterations"],
            shrink_factors=config["ants"]["linear"]["shrink_factors"],
            smoothing_factors=config["ants"]["linear"]["smoothing_factors"],
        ),
        linear_metric=lambda wildcards, input: "-m MI[{template},{target},1,32,Regular,0.25]".format(
            template=input.ref, target=input.flo
        ),
        deform_model="-t {deform_model}".format(
            deform_model=config["ants"]["deform"]["transform_model"]
        ),
        deform_multires="-c [{reg_iterations},1e-9,10] -f {shrink_factors} -s {smoothing_factors}".format(
            reg_iterations=config["ants"]["deform"]["reg_iterations"],
            shrink_factors=config["ants"]["deform"]["shrink_factors"],
            smoothing_factors=config["ants"]["deform"]["smoothing_factors"],
        ),
        deform_metric=lambda wildcards, input: "-m {metric}[{template},{target},1,4]".format(
            metric=config["ants"]["deform"]["sim_metric"],
            template=input.ref,
            target=input.flo,
        ),
    output:
        out_composite=bids(
            root=root,
            datatype="anat",
            suffix="Composite.h5",
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
        out_inv_composite=bids(
            root=root,
            datatype="anat",
            suffix="InverseComposite.h5",
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
        warped_flo=bids(
            root=root,
            datatype="anat",
            suffix="T1w.nii.gz",
            space="{template}",
            desc="SyN",
            **subj_wildcards
        ),
    threads: 8
    resources:
        mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time=60,  # 1 hrs
    container:
        config["singularity"]["prepdwi"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsRegistration {params.base_opts} {params.intensity_opts} "
        "{params.init_transform} "
        "{params.deform_model} {params.deform_metric} {params.deform_multires} "
        "-o [{params.out_prefix},{output.warped_flo}]"
"""


rule warp_dseg_from_template:
    input:
        dseg=lambda wildcards: workflow.source_path(
            os.path.join("..", "..", config["template_atlas_dseg"])
        ).format(**wildcards),
        ref=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        inv_warp=bids(
            root=root,
            datatype="anat",
            suffix="invwarp.nii.gz",
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
        affine_xfm_ras=bids(
            root=root,
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="ras",
            **subj_wildcards
        ),
    output:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="dseg.nii.gz",
            atlas="{atlas}",
            from_="{template}",
            reg="SyN"
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.dseg} -o {output.dseg} -r {input.ref} "
        " -t [{input.affine_xfm_ras},1] {input.inv_warp} "
        #use inverse xfm (going from template to subject)


rule warp_tissue_probseg_from_template:
    input:
        probseg=lambda wildcards: workflow.source_path(
            os.path.join("..", "..", config["template_tissue_probseg"])
        ).format(**wildcards),
        ref=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        inv_warp=bids(
            root=root,
            datatype="anat",
            suffix="invwarp.nii.gz",
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
        affine_xfm_ras=bids(
            root=root,
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="ras",
            **subj_wildcards
        ),
    output:
        probseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="probseg.nii.gz",
            label="{tissue}",
            from_="{template}",
            reg="SyN"
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.probseg} -o {output.probseg} -r {input.ref} "
        " -t [{input.affine_xfm_ras},1] {input.inv_warp} "
        #use inverse xfm (going from template to subject)


rule warp_brainmask_from_template:
    input:
        mask=lambda wildcards: workflow.source_path(
            os.path.join("..", "..", config["template_mask"])
        ).format(**wildcards),
        ref=bids(root=root, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        inv_composite=bids(
            root=root,
            datatype="anat",
            suffix="InverseComposite.h5",
            from_="subject",
            to="{template}",
            **subj_wildcards
        ),
    output:
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="mask.nii.gz",
            from_="{template}",
            reg="SyN",
            desc="brain"
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.mask} -o {output.mask} -r {input.ref} "
        " -t {input.inv_composite} "
        #use inverse xfm (going from template to subject)


rule dilate_brainmask:
    input:
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="mask.nii.gz",
            from_="{template}",
            reg="{desc}",
            desc="brain"
        ),
    params:
        dil_opt=" ".join(
            ["-dilate 1 3x3x3vox" for i in range(config["n_init_mask_dilate"])]
        ),
    output:
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="mask.nii.gz",
            from_="{template}",
            reg="{desc}",
            desc="braindilated"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} {params.dil_opt} -o {output}"


# dilate labels N times to provide more of a fudge factor when assigning GM labels
rule dilate_atlas_labels:
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="dseg.nii.gz",
            atlas="{atlas}",
            from_="{template}"
        ),
    params:
        dil_opt=" ".join(
            ["-dilate 1 3x3x3vox" for i in range(config["n_atlas_dilate"])]
        ),
    output:
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="dseg.nii.gz",
            atlas="{atlas}",
            from_="{template}",
            desc="dilated"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} {params.dil_opt} -o {output}"


rule resample_mask_to_dwi:
    input:
        mask=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            suffix="mask.nii.gz",
            from_="{template}",
            reg="SyN",
            desc="brain"
        ),
        ref=bids(
            root=root,
            desc="topup",
            datatype="dwi",
            method="jac",
            **subj_wildcards,
            suffix="b0.nii.gz"
        ),
    params:
        interpolation="NearestNeighbor",
    output:
        mask=bids(
            root=root,
            **subj_wildcards,
            desc="brain",
            suffix="mask.nii.gz",
            method="template",
            from_="{template}",
            reg="SyN",
            datatype="dwi"
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 --input-image-type 0 --input {input.mask} --reference-image {input.ref}  --interpolation {params.interpolation} --output {output.mask} --verbose"
