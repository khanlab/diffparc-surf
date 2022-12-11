"""  
    Note: This rule file is not currently used..


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
        datatype="surf",
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
            suffix="maxprob.label.gii",
            **subj_wildcards,
        ),
    params:
        tck_filename=get_tck_filename,
        label_num=lambda wildcards: config["targets"][wildcards.targets][  #gets label number from the parc name, with 1-based indexing
            "labels"
        ].index(wildcards.parc)
        + 1,
    output:
        tcklist=temp(
            bids(
                root=root,
                datatype="surf",
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


rule extract_target_mask:
    input:
        dseg=get_dseg_targets_nii,
    params:
        label_num=lambda wildcards: config["targets"][wildcards.targets][
            "labels"
        ].index(wildcards.parc)
        + 1,
        #start at 1
    output:
        mask=temp(
            bids(
                root=root,
                **subj_wildcards,
                space="individual",
                desc="{targets}",
                parc="{parc}",
                datatype="anat",
                suffix="mask.nii.gz"
            )
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input.dseg} -retain-labels {params.label_num} -binarize -o {output.mask}"


rule create_parc_bundle:
    """ create parc bundle from list of streamlines connected to the region.
    if there are no streamlines, then we simply touch the file - the next
    rule will check for a zero-sized file"""
    input:
        tcklist=bids(
            root=root,
            datatype="surf",
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
            datatype="surf",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="vertextracts",
            **subj_wildcards,
        ),
        mask=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            desc="{targets}",
            parc="{parc}",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    output:
        bundle=bids(
            root=root,
            datatype="surf",
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


rule create_parc_tdi:
    """tract density image for the parcel. if there are no streamlines,
    then we are given a zero-sized file (touched in previous rule),
    and if so, we create a zero-valued image as the tract density"""
    input:
        bundle=bids(
            root=root,
            datatype="surf",
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
            datatype="surf",
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
            datatype="surf",
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
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="tdi.nii.gz",
            **subj_wildcards,
        ),
    output:
        mask=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            parc="{parc}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="tdimask.nii.gz",
            **subj_wildcards,
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



