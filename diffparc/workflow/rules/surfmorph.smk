# to replicate earlier surface morphometry pipeline, we want to:
# 1. rigid reg target probseg to template probseg
# 2. run lddmm/greedy/ants to get mapping from template to target
# 3. apply transform to the template surface
# 4. calculate displacement based on target minus template
# 5. normalize by the average displacement vector in local neighborhood
#   - (was 10mm radius cube before, now we use 8mm fwhm surf smoothing - can tune this hyperparameter)
# 6. calculate inward/outward displacement using dot product with surface normal

# note: with synthseg -> template-shape-injection, this could potentially be optimized a bit more..


rule gen_template_surface:
    input:
        nii=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
        decimate_percent=config["surface_decimate_percent"],
    output:
        surf_gii=temp(
            get_template_prefix(
                root=root, subj_wildcards=subj_wildcards, template="{template}"
            )
            + "_hemi-{hemi}_desc-nostruct_{seed}.surf.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["pyvista"]
    script:
        "../scripts/gen_isosurface.py"


rule set_surface_structure:
    input:
        surf_gii=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template="{template}"
        )
        + "_hemi-{hemi}_desc-nostruct_{seed}.surf.gii",
    params:
        structure=lambda wildcards: config["hemi_to_structure"][wildcards.hemi],
    output:
        surf_gii=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template="{template}"
        )
        + "_hemi-{hemi}_{seed}.surf.gii",
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output} && "
        "wb_command -set-structure {output} {params.structure} -surface-type ANATOMICAL"


def get_subject_seg_for_shapereg(wildcards):
    if config["use_synthseg_seed"]:
        return (
            bids(
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
        )
    else:
        return (
            bids(
                root=root,
                **subj_wildcards,
                hemi="{hemi}",
                space="individual",
                label="{seed}",
                from_="{template}",
                datatype="anat",
                suffix="probseg.nii.gz"
            ),
        )


rule rigid_shape_reg:
    """ rigidly register subj to template shape using moment-invariants """
    input:
        template=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
        target=get_subject_seg_for_shapereg,
    params:
        general_opts="-d 3",
        rigid_opts="-m SSD -moments 2 -det 1",
    output:
        xfm_ras=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="ras",
            label="{seed}",
            datatype="surf",
            **subj_wildcards
        ),
        xfm_itk=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="itk",
            label="{seed}",
            datatype="surf",
            **subj_wildcards
        ),
        warped_target=bids(
            root=root,
            **subj_wildcards,
            desc="rigid",
            hemi="{hemi}",
            label="{seed}",
            space="{template}",
            datatype="surf",
            suffix="probseg.nii.gz"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} {params.rigid_opts}  -i {input.template} {input.target} -o {output.xfm_ras} &&  "
        "greedy -threads {threads} {params.general_opts} -rf {input.template} -rm {input.target} {output.warped_target}  -r {output.xfm_ras} && "
        "c3d_affine_tool {output.xfm_ras} -oitk {output.xfm_itk}"


rule fluid_shape_reg:
    """ fluid registration to get the displacements """
    input:
        template=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
        target=bids(
            root=root,
            **subj_wildcards,
            desc="rigid",
            hemi="{hemi}",
            label="{seed}",
            space="{template}",
            datatype="surf",
            suffix="probseg.nii.gz"
        ),
    output:
        warp=bids(
            root=root,
            datatype="surf",
            suffix="warp.nii.gz",
            hemi="{hemi}",
            from_="subject",
            to="{template}",
            label="{seed}",
            **subj_wildcards
        ),
        warped=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="fluid",
            label="{seed}",
            space="{template}",
            datatype="surf",
            suffix="probseg.nii.gz"
        ),
    threads: 8
    resources:
        mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time=60,  # 1 hrs
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "greedy -d 3 -threads {threads} -m SSD -i {input.template} {input.target} -o {output.warp} -n 100x50x10 && "
        "greedy -d 3 -threads {threads} -rf {input.template} -rm {input.target} {output.warped} -r {output.warp} "


rule convert_warpfield:
    input:
        warp=bids(
            root=root,
            datatype="surf",
            suffix="warp.nii.gz",
            hemi="{hemi}",
            from_="subject",
            to="{template}",
            label="{seed}",
            **subj_wildcards
        ),
    output:
        warp=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            suffix="surfwarp.nii.gz",
            to_="subject",
            from_="{template}",
            label="{seed}",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule deform_surface:
    input:
        surf_gii=rules.set_surface_structure.output,
        warp=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            suffix="surfwarp.nii.gz",
            to_="subject",
            from_="{template}",
            label="{seed}",
            **subj_wildcards
        ),
    output:
        surf_warped=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            desc="fluid",
            from_="{template}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-apply-warpfield {input.surf_gii} {input.warp} {output.surf_warped}"


rule compute_displacement_metrics:
    input:
        ref_surf=rules.set_surface_structure.output,
        comp_surf=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            desc="fluid",
            from_="{template}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    output:
        scalar=temp(
            bids(
                root=root,
                hemi="{hemi}",
                **subj_wildcards,
                from_="{template}",
                datatype="surf",
                label="{seed}",
                suffix="scalardisp.shape.gii"
            )
        ),
        vector=temp(
            bids(
                root=root,
                hemi="{hemi}",
                **subj_wildcards,
                from_="{template}",
                datatype="surf",
                label="{seed}",
                suffix="vectordisp.shape.gii"
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command  -surface-to-surface-3d-distance {input.comp_surf} {input.ref_surf} {output.scalar} -vectors {output.vector}"


rule calc_template_surf_normals:
    input:
        ref_surf=rules.set_surface_structure.output,
    output:
        normals=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template="{template}"
        )
        + "_hemi-{hemi}_label-{seed}_normals.shape.gii",
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-normals {input} {output}"


rule smooth_displacement_vectors:
    input:
        ref_surf=rules.set_surface_structure.output,
        vector=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    params:
        fwhm=config["surfdisp_normalization_fwhm_mm"],
    output:
        vector=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="smoothed",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -metric-smoothing {input.ref_surf} {input.vector} {params.fwhm} {output.vector} -fwhm"


rule normalize_displacement_by_smoothed:
    """subtracts the displacement by the smoothed displacement to provide a spatially-local
    displacement estimate"""
    input:
        vector=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
        smoothed=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="smoothed",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    output:
        normalized=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="normalized",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -metric-math 'vec - vecsmooth' {output.normalized}"
        " -var vec {input.vector} "
        " -var vecsmooth {input.smoothed}  "


rule calc_inout_displacement:
    """uses dot-product with surface normal to get inward/outward displacement"""
    input:
        vec=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            desc="normalized",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
        norm=rules.calc_template_surf_normals.output,
    output:
        inout=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="inout.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -metric-math '(v1*n1)+(v2*n2)+(v3*n3)' {output.inout}"
        " -var v1 {input.vec} -column 1  "
        " -var v2 {input.vec} -column 2  "
        " -var v3 {input.vec} -column 3  "
        " -var n1 {input.norm} -column 1  "
        " -var n2 {input.norm} -column 2  "
        " -var n3 {input.norm} -column 3  "


rule create_cifti_metric_dscalar:
    input:
        left_metric=bids(
            root=root,
            **subj_wildcards,
            hemi="L",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.shape.gii"
        ),
        right_metric=bids(
            root=root,
            **subj_wildcards,
            hemi="R",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.shape.gii"
        ),
    output:
        cifti_dscalar=bids(
            root=root,
            **subj_wildcards,
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.dscalar.nii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-create-dense-scalar {output} -left-metric {input.left_metric} -right-metric {input.right_metric}"


def get_cifti_dscalar_subjects(wildcards):
    filename = (
        bids(
            root=root,
            **subj_wildcards,
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.dscalar.nii"
        ),
    )

    if "session" in subj_wildcards:
        return sorted(
            expand(
                filename,
                subject=inputs.subjects,
                session=inputs.sessions,
                **wildcards,
                allow_missing=True
            )
        )
    else:
        return sorted(
            expand(filename, subject=inputs.subjects, **wildcards, allow_missing=True)
        )


rule merge_dscalar_metrics_over_subjects:
    input:
        cifti_dscalar_subjects=get_cifti_dscalar_subjects,
    params:
        merge_opt=lambda wildcards, input: " ".join(
            [f"-cifti {cifti}" for cifti in input.cifti_dscalar_subjects]
        ),
    output:
        cifti_dscalar_group=bids(
            root=root,
            subject="group",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.dscalar.nii",
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "group"
    shell:
        "wb_command -cifti-merge {output} {params.merge_opt}"
