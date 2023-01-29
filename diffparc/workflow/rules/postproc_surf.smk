
rule conn_csv_to_metric:
    input:
        csv=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
    params:
        normalize_percentile=lambda wildcards: config["seeds"][wildcards.seed][
            "normalize_percentile"
        ],
        seeds_per_vertex=lambda wildcards: float(wildcards.seedspervertex),
    output:
        gii_metric=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="nostructconn.shape.gii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/conn_csv_to_gifti_metric.py"


rule set_structure_conn_metric:
    input:
        gii_metric=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="nostructconn.shape.gii",
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
            method="{method}",
            suffix="conn.shape.gii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output} && "
        "wb_command -set-structure {output} {params.structure}"


rule smooth_conn_metric:
    """ this smooths the connectivity over the surface, to improve SNR 
    if we are using fewer streamlines """
    input:
        metric=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.shape.gii",
            **subj_wildcards,
        ),
        surf=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    output:
        metric=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            fwhm="{fwhm}mm",
            suffix="conn.shape.gii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -metric-smoothing {input.surf} {input.metric} {wildcards.fwhm} {output.metric} -fwhm"


rule create_cifti_conn_dscalar:
    input:
        left_metric=lambda wildcards: bids(
            root=root,
            datatype="surf",
            hemi="L",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            fwhm=config["seeds"][wildcards.seed]["fwhm_conn"],
            suffix="conn.shape.gii",
            **subj_wildcards,
        ).format(**wildcards),
        right_metric=lambda wildcards: bids(
            root=root,
            datatype="surf",
            hemi="R",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            fwhm=config["seeds"][wildcards.seed]["fwhm_conn"],
            suffix="conn.shape.gii",
            **subj_wildcards,
        ).format(**wildcards),
    output:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.dscalar.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-create-dense-scalar {output} -left-metric {input.left_metric} -right-metric {input.right_metric}"


rule create_cifti_conn_dscalar_maxprob:
    input:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.dscalar.nii",
            **subj_wildcards,
        ),
    output:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.dscalar.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-reduce {input} INDEXMAX {output}"


rule create_cifti_sumconn_dscalar:
    """ sum up connectivity at a voxel across targets, will be used to threshold"""
    input:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.dscalar.nii",
            **subj_wildcards,
        ),
    output:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="sumconn.dscalar.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-reduce {input} SUM {output}"


rule mask_maxprob_by_sumconn_threshold:
    """ use fraction of connectivity to threshold """
    input:
        maxprob=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.dscalar.nii",
            **subj_wildcards,
        ),
        sumconn=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="sumconn.dscalar.nii",
            **subj_wildcards,
        ),
    params:
        threshold=lambda wildcards: float(
            config["seeds"][wildcards.seed]["streamline_threshold_percent"]
        )
        / float(wildcards.seedspervertex),
    output:
        masked=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maskedmaxprob.dscalar.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-math '(SUMCONN > {params.threshold}) * MAXPROB' {output.masked} -var SUMCONN {input.sumconn} -var MAXPROB {input.maxprob}"


rule create_cifti_maxprob_dlabel:
    input:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maskedmaxprob.dscalar.nii",
            **subj_wildcards,
        ),
        label_list_txt=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["targets"][wildcards.targets]["label_list_txt"],
        ),
    output:
        cifti_dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-label-import {input.cifti_dscalar} {input.label_list_txt} {output.cifti_dlabel}"


rule split_cifti_maxprob_dlabel:
    """split the cifti dlabel into metric gii files"""
    input:
        cifti_dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    output:
        left_label_gii=bids(
            root=root,
            datatype="surf",
            hemi="L",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.label.gii",
            **subj_wildcards,
        ),
        right_label_gii=bids(
            root=root,
            datatype="surf",
            hemi="R",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.label.gii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-separate {input.cifti_dlabel} COLUMN "
        " -label CORTEX_LEFT {output.left_label_gii} "
        " -label CORTEX_RIGHT {output.right_label_gii}"


rule parcellate_cifti_metric:
    input:
        cifti_dscalar=bids(
            root=root,
            **subj_wildcards,
            datatype="surf",
            label="{seed}",
            method="{method}",
            suffix="{metric}.dscalar.nii"
        ),
        cifti_dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    output:
        cifti_pscalar=bids(
            root=root,
            **subj_wildcards,
            datatype="surf",
            label="{seed}",
            parcel="{targets}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="{metric}.pscalar.nii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-parcellate {input.cifti_dscalar} {input.cifti_dlabel} COLUMN "
        " {output}"
