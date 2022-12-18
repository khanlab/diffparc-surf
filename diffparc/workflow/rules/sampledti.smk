"""
    These rules sample dti metrics along the streamlines seeded from each 
    vertex using tcksample, aggregating along each streamline, then across 
    all streamlines at the vertex.

"""


rule tcksample_from_vertices:
    """ sample dti metric from streamlines connected to each vertex """
    input:
        tck_dir=bids(
            root=config["tmp_dir"],
            datatype="surf",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="vertextracts",
            **subj_wildcards,
        ),
        metric=bids(
            root=root,
            datatype="dwi",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
    params:
        stat_tck="-stat_tck {statsalong}",  #whether to use min, max, mean, median
    output:
        sampledti_dir=temp(
            directory(
                bids(
                    root=config["tmp_dir"],
                    datatype="surf",
                    hemi="{hemi}",
                    desc="{targets}",
                    label="{seed}",
                    seedspervertex="{seedspervertex}",
                    method="mrtrix",
                    suffix="sampledti",
                    metric="{metric}",
                    statsalong="{statsalong}",
                    **subj_wildcards,
                )
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "mkdir -p {output.sampledti_dir} && "
        "parallel --eta --jobs {threads} "
        "tcksample -nthreads 0 -quiet {input.tck_dir}/vertex_{{1}}.tck "
        " {input.metric} {output.sampledti_dir}/sample_{{1}}.txt "
        " {params.stat_tck} "
        " ::: `ls {input.tck_dir} | grep -Po '(?<=vertex_)[0-9]+'`"


rule sampledti_to_metric:
    """converts the tcksample txt files to a gifti metric"""
    input:
        sampledti_dir=bids(
            root=config["tmp_dir"],
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="mrtrix",
            suffix="sampledti",
            metric="{metric}",
            statsalong="{statsalong}",
            **subj_wildcards,
        ),
    output:
        gii_metric=temp(
            bids(
                root=root,
                datatype="surf",
                hemi="{hemi}",
                desc="{targets}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                method="mrtrix",
                metric="{metric}",
                statsalong="{statsalong}",
                statsacross="{statsacross}",
                suffix="nostructsampledti.shape.gii",
                **subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/tcksample_to_gifti_metric.py"


rule set_structure_sampledti_metric:
    input:
        gii_metric=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="mrtrix",
            metric="{metric}",
            statsalong="{statsalong}",
            statsacross="{statsacross}",
            suffix="nostructsampledti.shape.gii",
            **subj_wildcards,
        ),
    params:
        structure=lambda wildcards: config["hemi_to_structure"][wildcards.hemi],
    output:
        gii_metric=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="mrtrix",
            metric="{metric}",
            statsalong="{statsalong}",
            statsacross="{statsacross}",
            suffix="sampledti.shape.gii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output} && "
        "wb_command -set-structure {output} {params.structure}"


def get_gifti_metrics_sampledti(wildcards):
    files = dict()
    files["left_metric"] = bids(
        root=root,
        datatype="surf",
        hemi="L",
        desc="{targets}",
        label="{seed}",
        seedspervertex="{seedspervertex}",
        method="mrtrix",
        metric="{metric}",
        statsalong="{statsalong}",
        statsacross="{statsacross}",
        suffix="sampledti.shape.gii",
        **subj_wildcards,
    ).format(
        **wildcards,
        statsalong=config["stat_along_tcks"],
        statsacross=config["stat_across_tcks"],
    )
    files["right_metric"] = bids(
        root=root,
        datatype="surf",
        hemi="R",
        desc="{targets}",
        label="{seed}",
        seedspervertex="{seedspervertex}",
        method="mrtrix",
        metric="{metric}",
        statsalong="{statsalong}",
        statsacross="{statsacross}",
        suffix="sampledti.shape.gii",
        **subj_wildcards,
    ).format(
        **wildcards,
        statsalong=config["stat_along_tcks"],
        statsacross=config["stat_across_tcks"],
    )
    return files


rule create_cifti_sampledti_dscalar:
    input:
        unpack(get_gifti_metrics_sampledti),
    output:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            seedspervertex="{seedspervertex}",
            method="mrtrix",
            label="{seed}",
            suffix="{metric}.dscalar.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-create-dense-scalar {output} -left-metric {input.left_metric} -right-metric {input.right_metric}"
