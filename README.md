# Analysis for oral/gut phageome
Scripts for Methods described in [The human oral phageome is highly diverse and rich in jumbo phages](https://www.biorxiv.org/content/10.1101/2020.07.06.186817v1.full.pdf)

### In silico pipelines
- Generating the phage catalogue: [bash_scripts/commands.txt](https://github.com/APC-Microbiome-Ireland/phageome_analysis/tree/master/bash_scripts/commands.txt)
- Function annotation with HMMs: [bash_scripts/protein_function_annotation.sh](https://github.com/APC-Microbiome-Ireland/phageome_analysis/tree/master/bash_scripts/protein_function_annotation.sh)
- Identifying tRNA using ARAGORN: [bash_scripts/extract_trna_from_jumbophages.sh](https://github.com/APC-Microbiome-Ireland/phageome_analysis/tree/master/bash_scripts/extract_trna_from_jumbophages.sh)
- Identifying ARGs using CARD database: [bash_scripts/annotate_args.sh](https://github.com/APC-Microbiome-Ireland/phageome_analysis/tree/master/bash_scripts/annotate_args.sh)
- Other data cleaning/wraggling: [helper_scripts](https://github.com/APC-Microbiome-Ireland/phageome_analysis/tree/master/helper_scripts)

#### To Do
- Dependency scripts in commands.sh
- Command for vConTACT
- ../data/virsorter_positive.fasta file

### Data analysis and visualisation
- Statistical analysis and data visualisation: [phage_analysis.R](https://github.com/APC-Microbiome-Ireland/phageome_analysis/blob/master/phage_analysis.R)
