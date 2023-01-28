
rule calc_surface_area_metric:
    """this uses the template surface area"""
    input:
        template_surf=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_hemi-{hemi}_{seed}.surf.gii",
    output:
        metric=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="surf",
            label="{seed}",
            suffix="surfarea.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-vertex-areas {input} {output}"


rule calc_surface_area_ratio_metric:
    """ log ratio of subject surface area over template surface area,
    as a measure of expansion (+ve) /contraction (-ve) of the surface """
    input:
        template_surf=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_hemi-{hemi}_{seed}.surf.gii",
        subject_surf=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    output:
        metric=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="surf",
            label="{seed}",
            suffix="surfarearatio.shape.gii"
        ),
    group:
        "subj"
    shadow:
        "minimal"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-vertex-areas {input.template_surf} template_surfarea.shape.gii && "
        "wb_command -surface-vertex-areas {input.subject_surf} subject_surfarea.shape.gii && "
        "wb_command -metric-math 'log(SUBJECT/TEMPLATE)' {output} "
        "  -var TEMPLATE template_surfarea.shape.gii -var SUBJECT subject_surfarea.shape.gii "
