rule write_surf_metrics_tsv:
    input:
        dscalars = expand(bids(
                root=root,
                from_="{template}",
                datatype="surf",
                label="{seed}",
                suffix="{metric}.dscalar.nii",
                **subj_wildcards,
                ),
            metric=config['surface_metrics'], allow_missing=True),

        dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    output:
        dlabel=bids(
            root=root,
            datatype="surf",
            from_="{template}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="surfmetrics.tsv",
            **subj_wildcards,
        ),       
    container: config['singularity']['pythondeps']
    group: 'subj'
    script:
        '../scripts/write_surf_metrics_tsv.py'


