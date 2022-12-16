if config["in_prepdwi_dir"]:

    def get_bedpost_dir_from_prepdwi(wildcards):
        sub_ses_dir = bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        )
        bedpost_dir = os.path.join(
            config["in_prepdwi_dir"],
            "work",
            sub_ses_dir,
            "dwi",
            "uncorrected_denoise_unring_eddy_regT1",
            "bedpost.bedpostX",
        )
        return bedpost_dir.format(**wildcards)

    rule copy_bedpost_dir_from_prepdwi:
        input:
            bedpost_dir=get_bedpost_dir_from_prepdwi,
        output:
            bedpost_dir=directory(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion.bedpostX",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
        group:
            "subj"
        shell:
            "cp -R {input.bedpost_dir} {output.bedpost_dir}"


elif config["in_snakedwi_dir"]:

    rule copy_bedpost_dir_from_snakedwi:
        input:
            bedpost_dir=bids(
                root=config["in_snakedwi_dir"],
                desc="eddy",
                suffix="diffusion.bedpostX",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
        output:
            bedpost_dir=directory(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion.bedpostX",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
        group:
            "subj"
        shell:
            "cp -R {input.bedpost_dir} {output.bedpost_dir}"


else:

    rule copy_inputs_for_bedpost:
        input:
            dwi=bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="preproc",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            bval=bids(
                root=root,
                suffix="dwi.bval",
                desc="preproc",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            bvec=bids(
                root=root,
                suffix="dwi.bvec",
                desc="preproc",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            brainmask=bids(
                root=root,
                suffix="mask.nii.gz",
                desc="brain",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
        output:
            diff_dir=directory(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
            dwi=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "data.nii.gz",
            ),
            brainmask=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "nodif_brain_mask.nii.gz",
            ),
            bval=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "bvals",
            ),
            bvec=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "bvecs",
            ),
        group:
            "subj"
        shell:
            "mkdir -p {output.diff_dir} && "
            "cp {input.dwi} {output.dwi} && "
            "cp {input.brainmask} {output.brainmask} && "
            "cp {input.bval} {output.bval} && "
            "cp {input.bvec} {output.bvec} "
            #could symlink instead??


    def get_bedpost_cmd(wildcards):
        if config.get("use_gpu_bedpost_container", False):
            return (
                f"singularity exec --nv -e {config['singularity']['fsl_604']} bedpostx_gpu"
            )
        else:
            return os.path.join(workflow.basedir, f"scripts/bedpostx-parallel")


    def get_bedpost_parallel_opt(wildcards, threads):
        if config.get("use_gpu_bedpost_container", False):
            return ""
        else:
            return f"-P {threads}"

    rule run_bedpost:
        input:
            diff_dir=bids(
                root=root,
                desc="eddy",
                suffix="diffusion",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            dwi=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "data.nii.gz",
            ),
            brainmask=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "nodif_brain_mask.nii.gz",
            ),
            bval=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "bvals",
            ),
            bvec=os.path.join(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                "bvecs",
            ),
        params:
            bedpost_cmd=get_bedpost_cmd,
            parallel_opt=get_bedpost_parallel_opt,
        output:
            bedpost_dir=directory(
                bids(
                    root=root,
                    desc="eddy",
                    suffix="diffusion.bedpostX",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
        group:
            "subj"
        threads: 32  #this needs to be set in order to avoid multiple gpus from executing
        resources:
            gpus=1,
            mem_mb=16000,
            time=360,
        shell:
            "{params.bedpost_cmd} {input.diff_dir} {params.parallel_opt} && "
            "rm -rf {output.bedpost_dir}/logs && "
            "rm -rf {input.diff_dir}"
            #remove the logs to reduce # of files  
            # remove the input dir (copy of files) 
