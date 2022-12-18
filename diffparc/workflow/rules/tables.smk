def get_dscalar_nii(wildcards):
    metric = wildcards.metric
    if metric == "surfarea" or metric == "inout":
        dscalar = bids(
            root=root,
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.dscalar.nii",
            **subj_wildcards,
        )
    elif metric == "FA" or metric == "MD":
        dscalar = bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            label="{seed}",
            suffix="{metric}.dscalar.nii",
            **subj_wildcards,
        )
    return dscalar.format(**wildcards)


rule write_surf_metrics_legacy_csv:
    """ for backwards compatiblity with old diffparc - 
    separate file for each metric, using identical column names (parcels), 
    and index column as "subj", formatted as sub-{subject}_ses-{session} """
    input:
        dscalar=get_dscalar_nii,
        dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    params:
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            from_="{template}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            form="legacy",
            suffix="{metric}.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["pythondeps"]
    group:
        "subj"
    script:
        "../scripts/write_surf_metrics_legacy.py"


rule write_indepconn_metric_csv:
    """ this reads in the raw connectivity
    before normalization, and calculates the
    mean number of streamlines to each target
    from non-zero values"""
    input:
        csv_left=bids(
            root=root,
            datatype="surf",
            hemi="L",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
        csv_right=bids(
            root=root,
            datatype="surf",
            hemi="R",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
    params:
        target_labels=lambda wildcards: config["targets"][wildcards.targets]["labels"],
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            from_="{template}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            form="legacy",
            suffix="indepconn.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["pythondeps"]
    group:
        "subj"
    script:
        "../scripts/write_indepconn_metric.py"


"""
rule write_surf_metrics_long_csv:
    input:
        dscalars=expand(
            bids(
                root=root,
                from_="{template}",
                datatype="surf",
                label="{seed}",
                suffix="{metric}.dscalar.nii",
                **subj_wildcards,
            ),
            metric=config["surface_metrics"],
            allow_missing=True,
        ),
        dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method='{method}',
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    params:
        metrics=config["surface_metrics"],
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            from_="{template}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method='{method}',
            form="long",
            suffix="surfmetrics.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["pythondeps"]
    group:
        "subj"
    script:
        "../scripts/write_surf_metrics_long.py"


rule write_surf_metrics_wide_csv:
    input:
        dscalars=expand(
            bids(
                root=root,
                from_="{template}",
                datatype="surf",
                label="{seed}",
                suffix="{metric}.dscalar.nii",
                **subj_wildcards,
            ),
            metric=config["surface_metrics"],
            allow_missing=True,
        ),
        dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method='{method}',
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    params:
        metrics=config["surface_metrics"],
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            from_="{template}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method='{method}',
            form="wide",
            suffix="surfmetrics.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["pythondeps"]
    group:
        "subj"
    script:
        "../scripts/write_surf_metrics_wide.py"

"""


rule concat_subj_csv:
    input:
        csvs=expand(
            bids(
                root=root,
                datatype="tabular",
                from_="{template}",
                desc="{targets}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                form="{form}",
                suffix="{suffix}.csv",
                **subj_wildcards
            ),
            zip,
            **subj_zip_list,
            allow_missing=True
        ),
        #loop over subjects and sessions 
    output:
        csv=bids(
            root=root,
            subject="group",
            datatype="tabular",
            from_="{template}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            form="{form}",
            suffix="{suffix}.csv",
        ),
    container:
        config["singularity"]["pythondeps"]
    group:
        "agg"
    script:
        "../scripts/concat_csv.py"
