rule convert_rigid_to_world:
    input:
        xfm_itk=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="itk",
            label="{seed}",
            datatype="surf",
            **subj_wildcards
        ),
    output:
        rigid_world=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="world",
            label="{seed}",
            datatype="surf",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -convert-affine -from-itk {input} -to-world {output} -inverse"


rule transform_template_surf_to_t1:
    """ transforms the template surface to the subject T1 """
    input:
        surf_gii=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            desc="fluid",
            from_="{template}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
        rigid_world=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="world",
            label="{seed}",
            datatype="surf",
            **subj_wildcards
        ),
    output:
        surf_warped=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            hemi="{hemi}",
            from_="{template}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-apply-affine {input.surf_gii} {input.rigid_world} {output.surf_warped}"


# for seeding, create a csv with vertex coords
rule create_surf_seed_csv:
    input:
        surf=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            space="individual",
            from_="{template}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    output:
        csv=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            space="individual",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="seeds.csv"
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/surf_to_seed_csv.py"


rule track_from_vertices:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        wm_fod=get_fod_for_tracking,
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **subj_wildcards,
        ),
        mask=bids(
            root=root,
            datatype="dwi",
            suffix="mask.mif",
            **subj_wildcards,
        ),
        csv=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            hemi="{hemi}",
            from_=config["template"],
            datatype="surf",
            label="{seed}",
            suffix="seeds.csv"
        ),
    params:
        radius="0.5",
        seedspervertex="{seedspervertex}",
    output:
        tck_dir=temp(
            directory(
                bids(
                    root=config["tmp_dir"],
                    datatype="surf",
                    hemi="{hemi}",
                    label="{seed}",
                    seedspervertex="{seedspervertex}",
                    suffix="vertextracts",
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
        "mkdir -p {output.tck_dir} && "
        "parallel --bar --link --jobs {threads} "
        "tckgen -quiet -nthreads 0  -algorithm iFOD2 -mask {input.mask} "
        " {input.wm_fod} {output.tck_dir}/vertex_{{1}}.tck "
        " -seed_sphere {{2}},{params.radius} -seeds {params.seedspervertex} "
        " :::  `seq --format '%05g' $(cat {input.csv} | wc -l)` ::: `cat {input.csv}` "


def get_dseg_targets(wildcards):

    if config["use_synthseg_targets"]:
        return (
            bids(
                root=root,
                **subj_wildcards,
                space="individual",
                desc="{targets}",
                from_="synthsegnearest",
                datatype="anat",
                suffix="dseg.mif"
            ),
        )

    else:
        return (
            bids(
                root=root,
                **subj_wildcards,
                space="individual",
                desc="{targets}",
                from_=config["template"],
                datatype="anat",
                suffix="dseg.mif"
            ),
        )


def get_dseg_targets_nii(wildcards):

    if config["use_synthseg_targets"]:
        return (
            bids(
                root=root,
                **subj_wildcards,
                space="individual",
                desc="{targets}",
                from_="synthsegnearest",
                datatype="anat",
                suffix="dseg.nii.gz"
            ),
        )

    else:
        return (
            bids(
                root=root,
                **subj_wildcards,
                space="individual",
                desc="{targets}",
                from_=config["template"],
                datatype="anat",
                suffix="dseg.nii.gz"
            ),
        )


rule connectivity_from_vertices:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
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
        targets=get_dseg_targets,
    output:
        conn_dir=temp(
            directory(
                bids(
                    root=config["tmp_dir"],
                    datatype="surf",
                    hemi="{hemi}",
                    desc="{targets}",
                    label="{seed}",
                    seedspervertex="{seedspervertex}",
                    suffix="vertexconn",
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
        "mkdir -p {output.conn_dir} && "
        "parallel --eta --jobs {threads} "
        "tck2connectome -nthreads 0 -quiet {input.tck_dir}/vertex_{{1}}.tck {input.targets} {output.conn_dir}/conn_{{1}}.csv -vector"
        " ::: `ls {input.tck_dir} | grep -Po '(?<=vertex_)[0-9]+'`"


rule gen_vertex_conn_csv:
    input:
        conn_dir=bids(
            root=config["tmp_dir"],
            datatype="surf",
            desc="{targets}",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="vertexconn",
            **subj_wildcards,
        ),
    params:
        header_line=lambda wildcards: ",".join(
            config["targets"][wildcards.targets]["labels"]
        ),
    output:
        conn_csv=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/gather_csv_files.py"


rule conn_csv_to_metric:
    input:
        csv=bids(
            root=root,
            datatype="surf",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.csv",
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
                suffix="nostructconn.shape.gii",
                **subj_wildcards,
            )
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


rule create_cifti_conn_dscalar:
    input:
        left_metric=bids(
            root=root,
            datatype="surf",
            hemi="L",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.shape.gii",
            **subj_wildcards,
        ),
        right_metric=bids(
            root=root,
            datatype="surf",
            hemi="R",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.shape.gii",
            **subj_wildcards,
        ),
    output:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
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
            suffix="maxprob.dscalar.nii",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-reduce {input} INDEXMAX {output}"


# need to then convert that into a label, then can use parcellate
rule create_cifti_maxprob_dlabel:
    input:
        cifti_dscalar=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="maxprob.dscalar.nii",
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
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="{metric}.dscalar.nii"
        ),
        cifti_dlabel=bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ),
    output:
        cifti_pscalar=bids(
            root=root,
            **subj_wildcards,
            from_="{template}",
            datatype="surf",
            label="{seed}",
            parcel="{targets}",
            seedspervertex="{seedspervertex}",
            suffix="{metric}.pscalar.nii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-parcellate {input.cifti_dscalar} {input.cifti_dlabel} COLUMN "
        " {output}"


rule calc_surface_area_metric:
    input:
        surf_warped=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            hemi="{hemi}",
            from_="{template}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    output:
        metric=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            from_="{template}",
            datatype="surf",
            label="{seed}",
            suffix="surfarea.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-vertex-areas {input} {output}"


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
        label_num=lambda wildcards: config["targets"][wildcards.targets][
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
