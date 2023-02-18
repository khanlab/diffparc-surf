rule moco_dwi:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        bval=bids(
            root=root,
            suffix="dwi.bval",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        json=bids(
            root=root,
            suffix="dwi.json",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    #        brainmask=get_dwi_mask(),
    params:
        show_eta="--eta" if config["show_parallel_eta"] else "",
    output:
        affine_dir=directory(
            bids(
                root=root,
                suffix="transforms",
                desc="moco",
                datatype="dwi",
                **subj_wildcards
            )
        ),
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="moco",
            datatype="dwi",
            **subj_wildcards
        ),
        json=bids(
            root=root, suffix="dwi.json", datatype="dwi", desc="moco", **subj_wildcards
        ),
        bval=bids(
            root=root, suffix="dwi.bval", datatype="dwi", desc="moco", **subj_wildcards
        ),
    threads: 32
    shadow:
        "minimal"
    container:
        config["singularity"]["prepdwi"]  #-- this rule needs niftyreg, c3d and mrtrix
    group:
        "subj"
    shell:
        "c4d {input.dwi} -slice w 0:-1 -oo dwi_%03d.nii && "
        "parallel {params.show_eta} --jobs {threads} "
        "reg_aladin -flo dwi_{{1}}.nii  -ref dwi_000.nii -res warped_{{1}}.nii -aff affine_xfm_ras_{{1}}.txt "
        " ::: `ls dwi_???.nii | tail -n +2 | grep -Po '(?<=dwi_)[0-9]+'` && "
        " mkdir -p {output.affine_dir} && cp affine_xfm_ras_*.txt {output.affine_dir} && "
        " echo -e '1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1' > {output.affine_dir}/affine_xfm_ras_000.txt && "
        " mrcat dwi_000.nii warped_*.nii {output.dwi}  && "
        " cp {input.bval} {output.bval} && "
        " cp {input.json} {output.json} "


rule rotate_bvecs_moco:
    input:
        affine_dir=bids(
            root=root,
            suffix="transforms",
            desc="moco",
            datatype="dwi",
            **subj_wildcards
        ),
        bvec=bids(
            root=root,
            suffix="dwi.bvec",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    output:
        bvec=bids(
            root=root, suffix="dwi.bvec", datatype="dwi", desc="moco", **subj_wildcards
        ),
    params:
        script=os.path.join(workflow.basedir, "scripts/rotate_bvecs_multi.sh"),
    container:
        config["singularity"]["prepdwi"]  #fsl
    group:
        "subj"
    shell:
        "chmod a+x {params.script} && "
        "{params.script} {input.bvec} {output.bvec} {input.affine_dir}/affine_xfm_ras_*.txt"
