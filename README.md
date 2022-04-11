# Analysis for oral/gut phageome
Scripts for Methods described in [The human oral phageome is highly diverse and rich in jumbo phages](https://www.biorxiv.org/content/10.1101/2020.07.06.186817v1.full.pdf)

### Processing raw reads
The raw reads can be downlaoded from ENA using the accession numbers listed in the `ID` column of [Supplementary Table 1](data/Supplementary_Table1.csv) and also in the `Merged_with` column where processed reads are merged before assembly (e.g. ERR589353 ID is a merge of both ERR589353 and ERR589541).

The raw reads for all samples were trimmed using AlienTrimmer 0.4.0 using default parameters and [Illumina contaminant oligonucleotides](https://gitlab.pasteur.fr/aghozlan/shaman_bioblend/blob/18a17dbb44cece4a8320cce8184adb9966583aaa/alienTrimmerPF8contaminants.fasta). Human contaminant sequences were removed from all samples by discarding reads that mapped against a human reference genome (downloaded from Human Genome Resources at NCBI on 27th February 2017) using Bowtie2 2.2.3 with parameters `-N 1 -k 1 --end-to-end --very-sensitive –phred33 --no-discordant`.

### In silico pipelines
- Generating the phage catalogue: bash_scripts/commands.txt
- Function annotation with HMMs: bash_scripts/protein_function_annotation.sh
- Identifying tRNA using ARAGORN: bash_scripts/extract_trna_from_jumbophages.sh
- Identifying ARGs using CARD database: bash_scripts/annotate_args.sh
- Other data cleaning/wraggling: helper_scripts
  - Step 1: Unfiltered catalog* -> helper_scripts/create_catalogue_dataset.R -> Filtered catalog
  - Step 2: helper_scripts/get_candidate_jumbophages.R
  - Step 3: helper_scripts/collate_functional_annotations.R
  - Step 4: helper_scripts/phage_quantification.R

#### To Do
- Dependency scripts in commands.sh
- Command for vConTACT
- *Recreate the data/virsorter_positive.fasta unfiltered catalog and data/virsorter_positive.ids

### Data analysis and visualisation
Creating phage catalogue and quantification
- Statistical analysis and data visualisation: phage_analysis.R
