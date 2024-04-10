This repository contains the code and processed/supporting data related to the Hu et al. manuscript "A continuum of zinc finger transcription factor retention on native chromatin underlies dynamic genome organization". All raw and processed data generated in this study can be obtained at Genome Sequence Archive (GSA) under accession number: HRA005744.

- CTCF, MAZ, and ZNF143 motif files generated in this study are saved in `motif_file`
- native or non-native CTCF sites are saved in `CTCF_binding_sites`

## ChIP-seq & loMNase-seq

The N-ChIP, X-ChIP and loMNase-seq paired-end sequencing datasets from this study can be preprocessed using the `ChIP-seq/pipeline/Paired-end ChlP-seq_analysis_pipeline.sh` script.

	bash Paired-end_ChIP-seq_analysis_pipeline.sh -i "K562_NChIP_CTCF_75mMNaCl_R1 K562_NChIP_CTCF_75mMNaCl_R2 K562_NChIP_CTCF_75mMNaCl_R3" -p 8 -w ~/data

The public X-ChIP single-end sequencing datasets can be preprocessed using the `ChIP-seq/pipeline/Single-end ChlP-seq_analysis_pipeline.sh` script.

	bash Single-end ChlP-seq_analysis_pipeline.sh -i "K562_XChIP_MAZ_R1 K562_XChIP_MAZ_R2" -p 8 -w ~/data

Mandatory parameters:

	-i|--input chip_prefixes
			prefixes of ChIP-seq samples

	-p|--thread
			number of threads during processing

	-w|--workdir
			workdir

## Micro-C
The `Micro-C/pipeline/Micro-C_hicpro_hg38XX.config` file is the input file of HiC-Pro pipeline. 
The `Micro-C/pipeline/Micro-C_analysis_pipeline.sh` script can preprocess this study's Micro-C sequencing datasets. 

	bash Micro-C_analysis_pipeline.sh -i "K562_MicroC_DMSO_R1 K562_MicroC_DMSO_R2" -p 8 -w ~/data/

Mandatory parameters:

	-i|--input prefixes
			prefixes of Micro-C samples

	-p|--thread
			number of threads during processing

	-w|--workdir
			workdir
