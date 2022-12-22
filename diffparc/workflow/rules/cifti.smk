
rule create_cifti_metric_dscalar:
    input:
        left_metric=bids(
            root=root,
            **subj_wildcards,
            hemi="L",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.shape.gii"
        ),
        right_metric=bids(
            root=root,
            **subj_wildcards,
            hemi="R",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.shape.gii"
        ),
    output:
        cifti_dscalar=bids(
            root=root,
            **subj_wildcards,
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


rule merge_dscalar_metrics_over_subjects:
    input:
        cifti_dscalar_subjects=expand(
            bids(
                root=root,
                **subj_wildcards,
                datatype="surf",
                label="{seed}",
                suffix="{metric}.dscalar.nii"
            ),
            zip,
            **subj_zip_list,
            allow_missing=True
        ),
    params:
        merge_opt=lambda wildcards, input: " ".join(
            [f"-cifti {cifti}" for cifti in input.cifti_dscalar_subjects]
        ),
    output:
        cifti_dscalar_group=bids(
            root=root,
            subject="group",
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
