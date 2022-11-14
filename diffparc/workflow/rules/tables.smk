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
