rule convert_rigid_to_world:
    input:
        xfm_itk=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_=config["template"],
            to="subj",
            desc="rigid",
            type_="itk",
            label="{seed}",
            datatype="warps",
            **subj_wildcards
        ),
    output:
        rigid_world=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_=config["template"],
            to="subj",
            desc="rigid",
            type_="world",
            label="{seed}",
            datatype="warps",
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
            datatype="morph",
            suffix="{seed}.surf.gii"
        ),
        rigid_world=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_=config["template"],
            to="subj",
            desc="rigid",
            type_="world",
            label="{seed}",
            datatype="warps",
            **subj_wildcards
        ),
    output:
        surf_warped=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
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
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
    output:
        csv=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="tracts",
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
            hemi="{hemi}",
            datatype="tracts",
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
                    datatype="tracts",
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
    benchmark:
        bids(
            root="benchmarks",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="mrtrixsurftrack.tsv",
            **subj_wildcards,
        )
    shell:
        "mkdir -p {output.tck_dir} && "
        "parallel --eta --link --jobs {threads} "
        "tckgen -quiet -nthreads 0  -algorithm iFOD2 -mask {input.mask} "
        " {input.wm_fod} {output.tck_dir}/vertex_{{1}}.tck "
        " -seed_sphere {{2}},{params.radius} -seeds {params.seedspervertex} "
        " :::  `seq --format '%05g' $(cat {input.csv} | wc -l)` ::: `cat {input.csv}` "


rule connectivity_from_vertices:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
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
        targets=bids(
            root=root,
            **subj_wildcards,
            desc="{targets}",
            datatype="anat",
            suffix="dseg.mif"
        ),
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
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="mrtrix",
            suffix="conn.csv",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/gather_csv_files.py"
