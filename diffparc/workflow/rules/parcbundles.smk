"""  

     These rules create the streamline bundles that start from 
    a striatum parcellation, and reach the corresponding target, 
      e.g. a bundle with streamlines from the Caudal_motor striatal region 
        that end up reaching the Caudal_motor cortical region.
     
    This is done by first finding the specific streamline tck files
    that correspond to the vertices labelled by that region (create_parc_tcklist),
    then obtaining the subset of those streamlines that reach the cortical target
    region (create_parc_bundle). 
    
    Note: this can be zero streamlines -- if this is the case, the final bundle
    tck file is just a zero-sized file.. 

    One thing we can do with this is compute a volumetric tract density image 
    (create_parc_tdi), which we can threshold and use to mask with e.g. an FA 
    or MD map..

"""


def get_tck_filename(wildcards):
    return bids(
        root=config["tmp_dir"],
        datatype="tracts",
        hemi="{hemi}",
        label="{seed}",
        seedspervertex="{seedspervertex}",
        suffix="vertextracts/vertex_{{index:05d}}.tck",
        **subj_wildcards,
    ).format(**wildcards)


rule create_parc_tcklist:
    input:
        label_gii=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="mrtrix",
            suffix="maxprob.label.gii",
            **subj_wildcards,
        ),
    params:
        tck_filename=get_tck_filename,
        #gets label number from the parc name, with 1-based indexing
        label_num=lambda wildcards: config["targets"][wildcards.targets][
            "labels"
        ].index(wildcards.parc)
        + 1,
    output:
        tcklist=temp(
            bids(
                root=root,
                datatype="tracts",
                hemi="{hemi}",
                desc="{targets}",
                parc="{parc}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                suffix="tcklist.txt",
                **subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/create_parc_tcklist.py"


rule create_parc_bundle:
    """ create parc bundle from list of streamlines connected to the region.
    if there are no streamlines, then we simply touch the file - the next
    rule will check for a zero-sized file"""
    input:
        tcklist=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="tcklist.txt",
            **subj_wildcards,
        ),
        tck_dir=bids(
            root=config["tmp_dir"],
            datatype="tracts",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="vertextracts",
            **subj_wildcards,
        ),
        mask=bids(
            root=root,
            **subj_wildcards,
            targets="{targets}",
            desc="{parc}",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    output:
        bundle=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="bundle.tck",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "if [ `cat {input.tcklist} | wc -l` == 0 ]; "
        "then "
        "  touch {output.bundle}; "
        "else "
        "  tckedit `cat {input.tcklist}` {output.bundle} -include {input.mask}; "
        "fi"


rule sample_dti_from_parcbundle:
    input:
        bundle=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="bundle.tck",
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
        sample_txt=temp(
            bids(
                root=root,
                datatype="tracts",
                hemi="{hemi}",
                desc="{targets}",
                parc="{parc}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                method="mrtrix",
                suffix="sampledti.txt",
                metric="{metric}",
                statsalong="{statsalong}",
                **subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "tcksample -nthreads 0 -quiet {input.bundle} "
        " {input.metric} {output.sample_txt} "
        " {params.stat_tck} "


rule sampledti_to_metric:
    """converts the tcksample txt files to a gifti metric, setting all parc 
    vertices to the dti aggregate value"""
    input:
        sample_txts=lambda wildcards: expand(
            bids(
                root=root,
                datatype="tracts",
                hemi="{hemi}",
                desc="{targets}",
                parc="{parc}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                method="mrtrix",
                suffix="sampledti.txt",
                metric="{metric}",
                statsalong="{statsalong}",
                **subj_wildcards,
            ),
            parc=config["targets"][wildcards.targets]["labels"],
            **wildcards,
            allow_missing=True,
        ),
        label_gii=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.label.gii",
            **subj_wildcards,
        ),
    params:
        parcs=lambda wildcards: config["targets"][wildcards.targets]["labels"],
        sample_txt_file=lambda wildcards: bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="mrtrix",
            suffix="sampledti.txt",
            metric="{metric}",
            statsalong="{statsalong}",
            **subj_wildcards
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
                method="{method,mrtrix}",
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
        "../scripts/tcksample_parcs_to_gifti_metric.py"


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
            suffix="{metric,FA|MD}.dscalar.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-create-dense-scalar {output} -left-metric {input.left_metric} -right-metric {input.right_metric}"


rule create_parc_tdi:
    """tract density image for the parcel. if there are no streamlines,
    then we are given a zero-sized file (touched in previous rule),
    and if so, we create a zero-valued image as the tract density"""
    input:
        bundle=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="bundle.tck",
            **subj_wildcards,
        ),
        ref=bids(
            root=root,
            suffix="T1w.nii.gz",
            desc="preproc",
            datatype="anat",
            **subj_wildcards
        ),
    output:
        tdi=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="tdi.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "if [ -s {input.bundle} ]; "
        "then "
        "   tckmap {input.bundle} {output.tdi} -template {input.ref}; "
        "else "
        "   mrcalc {input.ref} 0 -mult {output.tdi}; "
        "fi"


rule threshold_tdi:
    """ threshold the tdi image using percentile of non-zero voxels.
    Note: we pass the bundle tck file to check if it is zero-sized"""
    input:
        bundle=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="bundle.tck",
            **subj_wildcards,
        ),
        tdi=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="tdi.nii.gz",
            **subj_wildcards,
        ),
    output:
        mask=temp(
            bids(
                root=root,
                datatype="tracts",
                hemi="{hemi}",
                desc="{targets}",
                parc="{parc}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                suffix="tdimask.nii.gz",
                **subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "if [ -s {input.bundle} ]; "
        "then "
        "   c3d {input.tdi} -pim ForegroundQuantile -threshold 90% +Inf 1 0 -o {output.mask}; "
        "else "
        "   c3d {input.tdi} -scale 0 -o {output.mask}; "
        "fi"
        #threshold
        #if no streamlines, just zero it out..


rule concat_all_streamlines:
    input:
        tck_dir=bids(
            root=config["tmp_dir"],
            datatype="tracts",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="vertextracts",
            **subj_wildcards,
        ),
    output:
        bundle=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="bundle.tck",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "tckedit `ls {input}/*.tck` {output.bundle}"
