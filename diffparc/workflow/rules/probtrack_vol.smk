rule binarize_upsampled_subject_seed:
    input:
        seed_nii=get_subject_seed_upsampled_probseg,
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
    output:
        seed_nii=bids(
            root=root,
            datatype="tracts",
            suffix="mask.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="seed",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -threshold 0.5 inf 1 0  -type uchar -o {output}"


rule fix_sform_seed:
    input:
        seed_nii=bids(
            root=root,
            datatype="tracts",
            suffix="mask.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="seed",
            **subj_wildcards
        ),
    output:
        seed_nii=temp(
            bids(
                root=root,
                datatype="tracts",
                suffix="mask.nii.gz",
                hemi="{hemi}",
                label="{seed}",
                desc="seed",
                fix="sform",
                **subj_wildcards
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["fsl"]
    shell:
        "cp {input} {output} && "
        "QFORM=`fslorient -getqform {output}` && "
        "fslorient -setsform $QFORM {output}"


rule create_seed_bg:
    input:
        seed_nii=bids(
            root=root,
            datatype="tracts",
            suffix="mask.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="seed",
            fix="sform",
            **subj_wildcards
        ),
    output:
        seed_nii=temp(
            bids(
                root=root,
                datatype="tracts",
                suffix="mask.nii.gz",
                hemi="{hemi}",
                label="{seed}",
                desc="seedbg",
                fix="sform",
                **subj_wildcards
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input} -replace 1 0 0 1 -o {output}"


rule run_probtrack_volume:
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
            targets="{targets}",
            datatype="tracts",
            desc="probtrack",
            suffix="targets.txt"
        ),
        seed_nii=bids(
            root=root,
            suffix="mask.nii.gz",
            datatype="tracts",
            hemi="{hemi}",
            label="{seed}",
            desc="seed",
            fix="sform",
            **subj_wildcards
        ),
        dwi_brain_mask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            fix="sform",
            datatype="dwi",
            **subj_wildcards
        ),
        seed_target_brain_mask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            fix="sform",
            datatype="dwi",
            **subj_wildcards
        ),
        targets=lambda wildcards: expand(
            bids(
                root=root,
                **subj_wildcards,
                targets="{targets}",
                desc="{desc}",
                fix="sform",
                datatype="anat",
                suffix="mask.nii.gz"
            ),
            desc=config["targets"][wildcards.targets]["labels"],
            allow_missing=True,
        ),
    params:
        seedspervoxel="{seedspervoxel}",
    output:
        out_tract_dir=directory(
            bids(
                root=root,
                **subj_wildcards,
                hemi="{hemi}",
                label="{seed}",
                desc="{targets}",
                seedspervoxel="{seedspervoxel}",
                datatype="tracts",
                suffix="probtrack"
            )
        ),
        out_conn_txt=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            desc="{targets}",
            seedspervoxel="{seedspervoxel}",
            datatype="tracts",
            suffix="probtrack/matrix_seeds_to_all_targets"
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=4000,
        time=lambda wildcards: int(0.2 * float(wildcards.seedspervoxel)),  # 15 minutes for 100 seedspervertex for undecimated striatum, so set at 0.20 minutes per seed (this will need to go up if a larger seed region is used)
    benchmark:
        bids(
            root="benchmarks",
            hemi="{hemi}",
            label="{seed}",
            desc="{targets}",
            seedspervoxel="{seedspervoxel}",
            suffix="probtrackvol.tsv",
            **subj_wildcards,
        )
    container:
        config["singularity"]["fsl"]
    shell:
        "probtrackx2 "
        " -x {input.seed_nii} "
        " -m {input.dwi_brain_mask} "
        " -s {input.bedpost_dir}/merged "
        " --dir={output.out_tract_dir} "
        " --targetmasks={input.target_txt} "
        " --forcedir "
        " --opd --os2t  --s2tastext "
        " --seedref={input.seed_target_brain_mask}"
        " --omatrix2 "
        " --target2={input.seed_target_brain_mask}"
        " --randfib=2 "
        " -V 0 "
        " -l  --onewaycondition -c 0.2 -S 2000 --steplength=0.5 "
        " -P {params.seedspervoxel} --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 "


rule merge_seed_conn_files:
    input:
        tract_dir=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            desc="{targets}",
            seedspervoxel="{seedspervoxel}",
            datatype="tracts",
            suffix="probtrack"
        ),
        seed_bg=bids(
            root=root,
            datatype="tracts",
            suffix="mask.nii.gz",
            hemi="{hemi}",
            label="{seed}",
            desc="seedbg",
            fix="sform",
            **subj_wildcards
        ),
    params:
        seed_conn_files=lambda wildcards, input: [
            f"{input.tract_dir}/seeds_to_{fname}"
            for fname in expand(
                bids(
                    include_subject_dir=False,
                    include_session_dir=False,
                    **subj_wildcards,
                    targets="{targets}",
                    desc="{desc}",
                    fix="sform",
                    suffix="mask.nii.gz",
                ),
                desc=config["targets"][wildcards.targets]["labels"],
                **wildcards,
            )
        ],
    output:
        conn_nii=bids(
            root=root,
            datatype="tracts",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="fsl",
            suffix="conn.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["fsl"]
    shell:
        "fslmerge -t {output} {input.seed_bg} {params.seed_conn_files}"
