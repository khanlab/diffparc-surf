rule upsample_template_probseg:
    input:
        nii=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
    params:
        resample=lambda wildcards: config["seeds"][wildcards.seed]["probseg_resample"],
    output:
        nii=temp(
            get_template_prefix(
                root=root, subj_wildcards=subj_wildcards, template=config["template"]
            )
            + "_hemi-{hemi}_desc-upsampled_label-{seed}_probseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input}  -resample {params.resample} -o {output}"


rule gen_template_surface:
    input:
        nii=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_hemi-{hemi}_desc-upsampled_label-{seed}_probseg.nii.gz",
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
        decimate_percent=lambda wildcards: config["seeds"][wildcards.seed][
            "surface_decimate_percent"
        ],
    output:
        surf_gii=temp(
            get_template_prefix(
                root=root, subj_wildcards=subj_wildcards, template=config["template"]
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
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_hemi-{hemi}_desc-nostruct_{seed}.surf.gii",
    params:
        structure=lambda wildcards: config["hemi_to_structure"][wildcards.hemi],
    output:
        surf_gii=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_hemi-{hemi}_{seed}.surf.gii",
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output} && "
        "wb_command -set-structure {output} {params.structure} -surface-type ANATOMICAL"


rule convert_surf_gii_to_vtk_polydata:
    """ generic rule for converting surf gii into vtk polydata, so we can calculate enclosed volume of the surface"""
    input:
        surf_gii="{file}.surf.gii",
    output:
        surf_vtk=temp("{file}.surf.vtk"),
    group:
        "subj"
    container:
        config["singularity"]["pyvista"]
    script:
        "../scripts/convert_surf_gii_to_vtk.py"


rule calc_surf_volume:
    """calcs enclosed volume from surface, writes single line to a text file"""
    input:
        surf_vtk="{file}.surf.vtk",
    output:
        txt=temp("{file}.surfvolume.txt"),
    group:
        "subj"
    container:
        config["singularity"]["pyvista"]
    script:
        "../scripts/calculate_surf_volume_vtk.py"


rule convert_affine_to_world:
    input:
        xfm_itk=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="itk",
            **subj_wildcards
        ),
    output:
        xfm_world=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_=config["template"],
            to_="subject",
            desc="world",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -convert-affine -from-itk {input} -to-world {output} -inverse"


rule linear_transform_surf_to_template:
    """this is used to obtain the volmni (ie affine-normalized) surfaces
    that are somewhat adjusted for head size"""
    input:
        surf=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
        xfm_world=bids(
            root=root,
            datatype="warps",
            suffix="affine.txt",
            from_=config["template"],
            to_="subject",
            desc="world",
            **subj_wildcards
        ),
    output:
        surf=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="surf",
            space=config["template"],
            warp="linear",
            suffix="{seed}.surf.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-apply-affine {input.surf} {input.xfm_world} {output.surf}"
