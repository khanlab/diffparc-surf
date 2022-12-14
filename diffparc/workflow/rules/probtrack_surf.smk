rule extract_target_mask:
    input: 
        targets=get_dseg_targets_nii,
    params:
        label_num=lambda wildcards: config["targets"][wildcards.targets][
            "labels"
        ].index(wildcards.desc)+1
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
    container: config['singularity']['itksnap']
    group: 'subj'
    shell:  
        "c3d {input} -retain-labels {params.label_num} -binarize -o {output}"

rule gen_targets_txt:
    input:
        targets=lambda wildcards: expand(bids(
            root=root,
            **subj_wildcards,
            space="individual",
            targets="{targets}",
            desc="{desc}",
            from_=config["template"],
            datatype="anat",
            suffix="mask.nii.gz"
        ),desc=config['targets'][wildcards.targets]['labels'],allow_missing=True)
    output:
        target_txt = bids(
            root=root,
            **subj_wildcards,
            space="individual",
            targets="{targets}",
            from_=config["template"],
            datatype="anat",
            suffix="targets.txt"
        )
    group: 'subj'
    run:
        f = open(output.target_txt,'w')
        for s in input.targets:
            f.write(f'{s}\n')
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
        seedref=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
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
            from_="{template}",
            datatype="surf",
            suffix="{seed}.surf.gii"
        ),

    
    params:
        seeds_per_vertex=config['seeds_per_vertex']
    output:
        out_tract_dir=directory(bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            label="{seed}",
            desc='{targets}',
            from_="{template}",
            datatype="surf",
            suffix="probtrack"
        ),

)
    shell:
        "probtrackx2 "
        " -x {input.surf_gii} "
        " -m {input.bedpost_dir}/nodif_brain_mask.nii.gz "
        " -s {input.bedpost_dir}/merged "
        " --dir={output.out_tract_dir} "
        " --targetmasks={input.target_txt} "
        " --forcedir "
        " --opd --os2t  --s2tastext "
        " --seedref={input.seedref} "
        " --omatrix2 "
        " --target2={input.bedpost_dir}/nodif_brain_mask.nii.gz "
        " --randfib=2 "
        " -V 1 "
        " -l  --onewaycondition -c 0.2 -S 2000 --steplength=0.5 "
        " -P {params.seeds_per_vertex} --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 "


