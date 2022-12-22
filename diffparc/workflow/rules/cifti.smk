
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


def get_cmd_spec_file(wildcards, input, output):
    specfile = output.spec_file

    cmds = list()
    for infile in input.L:
        cmds.append(
            " ".join(
                [
                    "wb_command",
                    "-add-to-spec-file",
                    specfile,
                    config["hemi_to_structure"]["L"],
                    infile,
                ]
            )
        )
    for infile in input.R:
        cmds.append(
            " ".join(
                [
                    "wb_command",
                    "-add-to-spec-file",
                    specfile,
                    config["hemi_to_structure"]["R"],
                    infile,
                ]
            )
        )
    for infile in input.LR:
        cmds.append(
            " ".join(
                [
                    "wb_command",
                    "-add-to-spec-file",
                    specfile,
                    config["hemi_to_structure"]["R"],
                    infile,
                ]
            )
        )

    for infile in input.other:
        cmds.append(
            " ".join(
                [
                    "wb_command",
                    "-add-to-spec-file",
                    specfile,
                    "OTHER",
                    infile,
                ]
            )
        )

    return " && ".join(cmds)


def get_inputs_spec_file(wildcards):
    inputs_dict = {"L": [], "R": [], "LR": [], "other": []}

    dti_metrics = list(set(["FA", "MD"]).intersection(set(config["surface_metrics"])))

    for hemi in config["hemispheres"]:
        inputs_dict[hemi].extend(
            expand(
                bids(
                    root=root,
                    datatype="surf",
                    hemi=hemi,
                    suffix="{seed}.surf.gii",
                    **subj_wildcards,
                ),
                **wildcards,
            )
        )

    morph_suffixes = [
        f"{metric}.dscalar.nii"
        for metric in list(
            set(["surfarea", "inout"]).intersection(set(config["surface_metrics"]))
        )
    ]

    inputs_dict["LR"].extend(
        expand(
            bids(
                root=root,
                datatype="surf",
                label="{seed}",
                suffix="{suffix}",
                **subj_wildcards,
            ),
            suffix=morph_suffixes,
            **wildcards,
        )
    )

    conn_suffixes = ["maxprob.dlabel.nii", "conn.dscalar.nii"]
    conn_suffixes.extend([f"{dti}.dscalar.nii" for dti in dti_metrics])

    inputs_dict["LR"].extend(
        expand(
            bids(
                root=root,
                datatype="surf",
                desc="{targets}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                method="{method}",
                suffix="{suffix}",
                **subj_wildcards,
            ),
            targets=config["seeds"][wildcards.seed]["targets"],
            seedspervertex=config["seeds"][wildcards.seed]["seeds_per_vertex"],
            method=config["methods"],
            suffix=conn_suffixes,
            **wildcards,
        )
    )

    inputs_dict["other"].extend(
        expand(
            bids(
                root=root,
                datatype="anat",
                desc="preproc",
                suffix="T1w.nii.gz",
                **subj_wildcards,
            ),
            **wildcards,
        )
    )

    inputs_dict["other"].extend(
        expand(
            bids(
                root=root,
                datatype="dwi",
                suffix="{dti}.nii.gz",
                **subj_wildcards,
            ),
            dti=dti_metrics,
            **wildcards,
        )
    )

    return inputs_dict


rule create_spec:
    input:
        unpack(get_inputs_spec_file),
    params:
        cmd=get_cmd_spec_file,
    output:
        spec_file=bids(
            root=root, datatype="surf", suffix="{seed}.spec", **subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"
