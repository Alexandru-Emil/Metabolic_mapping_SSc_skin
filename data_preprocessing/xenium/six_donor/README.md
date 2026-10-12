# Original six-donor Xenium workflow

Inputs: the QC-filtered six-donor SFE, raw counts and spatial coordinates. Apply the original six-donor QC workflow first. This folder uses unsupervised ResolVI run 17922984, not the seven-donor or semi-supervised alternatives.

1. Run `prepare_resolvi_input.R input_sfe.rds input.h5ad` with Rscript.
2. Create the Python environment from `environment.yml`; run `train_resolvi.py --input input.h5ad --output resolvi.h5ad`.
3. Render `annotate_resolvi.Rmd` with absolute `input_h5ad` and `output_rds` paths, for example: `rmarkdown::render("data_preprocessing/xenium/six_donor/annotate_resolvi.Rmd", params = list(input_h5ad = "/path/resolvi.h5ad", output_rds = "/path/annotation.rds"))`.
4. Run `attach_annotations.R input_sfe.rds annotation.rds annotated_sfe.rds` with Rscript.

The final annotation comes from Yi-Nan Li's 25 February 2026 saved report and includes the 10 August 2026 muscle-label revision. Reference R environment: R 4.5.2, Seurat 5.4.0, SeuratObject 5.3.0, anndataR 1.0.1, dplyr 1.2.0. Changing versions requires explicit opt-in. Existing downstream workflows may also require manuscript-specific annotation columns such as `lv2_anno_new`.

Sources: [ssc-metabolic](https://github.com/ynanli/ssc-metabolic), original trainer commit `c0b9ce6`; final report `XeniumResolviEcFibAnno_20260225_234739.html`.
