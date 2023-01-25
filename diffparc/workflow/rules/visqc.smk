
rule qc_reg:
    input:
        ref=config["template_t1w"],
        flo=bids(
            root=root,
            subject="{subject}",
            suffix="T1w.nii.gz",
            space="{template}",
            desc="{desc}",
        ),
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                suffix="regqc.png",
                from_="subject",
                to="{template}",
                desc="{desc}",
            ),
            caption="../report/regqc.rst",
            category="Registration QC",
            subcategory="{desc} {template}",
        ),
        html=bids(
            root="qc",
            subject="{subject}",
            suffix="regqc.html",
            from_="subject",
            to="{template}",
            desc="{desc}",
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/vis_regqc.py"


rule qc_probseg:
    input:
        img=bids(root=root, subject="{subject}", desc="n4", suffix="T1w.nii.gz"),
        seg4d=bids(
            root=root,
            subject="{subject}",
            suffix="probseg.nii.gz",
            desc="atropos3seg",
        ),
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                suffix="probseg.png",
                desc="atropos3seg",
            ),
            caption="../report/segqc.rst",
            category="Segmentation QC",
            subcategory="3-class Tissue Segmentation",
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/vis_qc_probseg.py"


rule qc_dseg:
    input:
        img=bids(root=root, subject="{subject}", desc="n4", suffix="T1w.nii.gz"),
        seg=bids(
            root=root,
            subject="{subject}",
            suffix="dseg.nii.gz",
            atlas="{atlas}",
            from_="{template}",
            reg="SyN",
        ),
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                suffix="dseg.png",
                atlas="{atlas}",
                from_="{template}",
            ),
            caption="../report/segqc.rst",
            category="Segmentation QC",
            subcategory="{atlas} Atlas from {template}",
        ),
        html=bids(
            root="qc",
            subject="{subject}",
            suffix="dseg.html",
            atlas="{atlas}",
            from_="{template}",
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/vis_qc_dseg.py"

rule qc_structure:
    input:
        surf_mesh=expand(
            bids(
                root=root,
                datatype="surf",
                hemi="{hemi}",
                suffix="{seed}.surf.gii",
                **subj_wildcards,
            ),
            hemi=["L", "R"],
            allow_missing=True,
        ),
        surf_roi=expand(
            bids(
                root=root,
                datatype="surf",
                hemi="{hemi}",
                desc="{targets}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                method="{method}",
                suffix="maxprob.dlabel.nii",
                **subj_wildcards,
            ),
            hemi=["L", "R"],
            allow_missing=True,
        ),
        vol_nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        ),
        vol_roi=expand(
            bids(
                root=root,
                datatype="anat",
                hemi="{hemi}",
                desc="{targets}",
                label="{seed}",
                seedspervoxel="{seedspervertex}",
                method="{method}",
                segtype="maxprob",
                suffix="dseg.nii.gz",
                **subj_wildcards,
            ),
            hemi=["L", "R"],
            allow_missing=True,
        ),
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                session="{session}",
                desc="{targets}",
                method="{method}",
                seedspervertex="{seedspervertex}",
                suffix="{seed}QC.png",
            ),
            caption="../report/seg_qc.rst",
            category="Segmentation QC",
            subcategory="{seed} to {targets}",
        )
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/vis_qc_seg.py"

rule qc_synthseg:
    input:
        vol_nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        ),
        synthseg_dseg=bids(
            root=root,
            datatype="anat",
            desc="synthseg",
            suffix="dseg.nii.gz",
            **subj_wildcards,
        )
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                session="{session}",
                desc="synthseg",
                suffix="dsegQC.png"
            ),
            caption="../report/synthseg_qc.rst",
            category="Synthseg QC"
        )
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/vis_qc_synthseg.py"