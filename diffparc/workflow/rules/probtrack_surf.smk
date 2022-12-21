# TODO: add FSL containers


rule extract_target_mask:
    input:
        targets=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            desc="{targets}",
            from_=config["template"],
            datatype="anat",
            suffix="dseg.nii.gz"
        ),
    params:
        label_num=lambda wildcards: config["targets"][wildcards.targets][
            "labels"
        ].index(wildcards.desc)
        + 1,
    output:
        bids(
            root=root,
            **subj_wildcards,
            space="individual",
            targets="{targets}",
            desc="{desc}",
            from_=config["template"],
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -retain-labels {params.label_num} -binarize -o {output}"


rule fix_sform_mask:
    input:
        brain_mask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        brain_mask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            fix="sform",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input} {output} && "
        "QFORM=`fslorient -getqform {output}` && "
        "fslorient -setsform $QFORM {output}"


rule fix_sform_target:
    input:
        bids(
            root=root,
            **subj_wildcards,
            space="individual",
            targets="{targets}",
            desc="{desc}",
            from_=config["template"],
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    output:
        bids(
            root=root,
            **subj_wildcards,
            space="individual",
            targets="{targets}",
            desc="{desc}",
            from_=config["template"],
            fix="sform",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    group:
        "subj"
    shell:
        "cp {input} {output} && "
        "QFORM=`fslorient -getqform {output}` && "
        "fslorient -setsform $QFORM {output}"


rule gen_targets_txt:
    input:
        targets=lambda wildcards: expand(
            bids(
                root=root,
                **subj_wildcards,
                space="individual",
                targets="{targets}",
                desc="{desc}",
                from_=config["template"],
                fix="sform",
                datatype="anat",
                suffix="mask.nii.gz"
            ),
            desc=config["targets"][wildcards.targets]["labels"],
            allow_missing=True,
        ),
    output:
        target_txt=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            targets="{targets}",
            from_=config["template"],
            datatype="anat",
            suffix="targets.txt"
        ),
    group:
        "subj"
    run:
        f = open(output.target_txt, "w")
        for s in input.targets:
            f.write(f"{s}\n")
        f.close()


rule run_probtrack_surface:
    input:
        bedpost_dir=bids(
            root=root,
            desc="eddy",
            suffix="diffusion.bedpostX",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        target_txt=bids(
            root=root,
            **subj_wildcards,
            space="individual",
            targets="{targets}",
            from_=config["template"],
            datatype="anat",
            suffix="targets.txt"
        ),
        surf_gii=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            space="individual",
            from_=config["template"],
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),
        brain_mask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            fix="sform",
            datatype="dwi",
            **subj_wildcards
        ),
    params:
        seeds_per_vertex="{seedspervertex}",
    output:
        out_tract_dir=directory(
            bids(
                root=root,
                **subj_wildcards,
                hemi="{hemi}",
                label="{seed}",
                desc="{targets}",
                seedspervertex="{seedspervertex}",
                datatype="surf",
                suffix="probtrack"
            )
        ),
        out_conn_txt=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            desc="{targets}",
            seedspervertex="{seedspervertex}",
            datatype="surf",
            suffix="probtrack/matrix_seeds_to_all_targets"
        ),
    group:
        "subj"
    shell:
        "probtrackx2 "
        " -x {input.surf_gii} "
        " -m {input.brain_mask} "
        " -s {input.bedpost_dir}/merged "
        " --dir={output.out_tract_dir} "
        " --targetmasks={input.target_txt} "
        " --forcedir "
        " --opd --os2t  --s2tastext "
        " --seedref={input.brain_mask}"
        " --omatrix2 "
        " --target2={input.brain_mask}"
        " --randfib=2 "
        " -V 0 "
        " -l  --onewaycondition -c 0.2 -S 2000 --steplength=0.5 "
        " -P {params.seeds_per_vertex} --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 "


rule create_conn_csv_probtrack:
    input:
        conn_txt=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            desc="{targets}",
            seedspervertex="{seedspervertex}",
            datatype="surf",
            suffix="probtrack/matrix_seeds_to_all_targets"
        ),
    params:
        col_headers=lambda wildcards: config["targets"][wildcards.targets]["labels"],
    output:
        conn_csv=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="fslprobtrack",
            datatype="surf",
            suffix="conn.csv"
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/probtrack_matrix_to_csv.py"
