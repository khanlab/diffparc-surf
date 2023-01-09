def get_dscalar_nii(wildcards):
    metric = wildcards.metric
    if metric == "surfarea" or metric == "inout":
        dscalar = bids(
            root=root,
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


def get_maxprob_dsegs(wildcards):
    out_dict = {}
    if wildcards.metric == "vol":
        for hemi in config["hemispheres"]:
            out_dict[f"maxprob_{hemi}"] = bids(
                root=root,
                datatype="anat",
                hemi=hemi,
                desc="{targets}",
                label="{seed}",
                seedspervoxel="{seedspervoxel}",
                method="{method}",
                segtype="maxprob",

                suffix="dseg.nii.gz",
                **subj_wildcards,
            ).format(**wildcards)
    elif wildcards.metric == "volmni":
        for hemi in config["hemispheres"]:
            out_dict[f"maxprob_{hemi}"] = bids(
                root=root,
                datatype="anat",
                hemi=hemi,
                space=config["template"],
                warp="linear",
                desc="{targets}",
                label="{seed}",
                seedspervoxel="{seedspervoxel}",
                method="{method}",
                segtype="maxprob",
                suffix="dseg.nii.gz",
                **subj_wildcards,
            ).format(**wildcards)
    return out_dict


rule write_vol_metrics_legacy_csv:
    """ for backwards compatiblity with old diffparc - 
    separate file for each metric, using identical column names (parcels), 
    and index column as "subj", formatted as sub-{subject}_ses-{session} """
    input:
        unpack(get_maxprob_dsegs),
    params:
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
        parc_list=lambda wildcards: config["targets"][wildcards.targets]["labels"],
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            form="legacy",
            suffix="{metric,vol|volmni}.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["pythondeps"]
    group:
        "subj"
    script:
        "../scripts/write_vol_metrics_legacy.py"


rule write_indepconn_metric_csv:
    """ this reads in the raw connectivity
    before normalization, and calculates the
    mean number of streamlines to each target
    from non-zero values"""
    input:
        csv_left=bids(
            root=root,
            datatype="tracts",
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
            datatype="tracts",
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


rule concat_subj_csv:
    input:
        csvs=expand(
            bids(
                root=root,
                datatype="tabular",
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
