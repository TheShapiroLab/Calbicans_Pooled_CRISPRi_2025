# Calbicans_Pooled_CRISPRi_2025
Scripts for data analysis and visualization for Wensing, Després et al, 2025


Read processing and sgRNA counting was performed on the Digital Research Alliance of Canada computing clusters. The notebooks used to generate the required scripts and analyze counts are provided as is and might require changes depending on your specific programming environment. The basic environment used for these analysis contains the following packages: 
- contourpy-1.2.1
- cycler-0.12.1
- fonttools-4.53.0
- kiwisolver-1.4.5
- matplotlib-3.9.0
- mpmath-1.3.0
- numpy-1.26.4
- packaging-24.1
- pandas-2.2.1
- pillow-10.3.0
- pyparsing-3.1.2
- python_dateutil-2.9.0.post0
- pytz-2024.1
- scipy-1.13.1
- six-1.16.0
- sympy-1.12.1
- tzdata-2024.1 
- seaborn-0.13.2

Downstream analysis and visualization was performed in Jupyter Notebooks, using the following packages: 
- numpy-v1.26.4
- matplotlib-3.9.2
- scipy-1.13.1
- pandas-2.2.2
- seaborn-0.13.2
- biopython-1.85

The repository is divided into five directories:
- 00_sgRNA_annot: Annotates sgRNAs in the library with different information on TSS, GC content, SNPs.
- 01_timecourse_assay: Analyze CRISPRi data from the timecourse assay and generate figure panels
- 02_environmental_screen: Analyze CRISPRi data from the timecourse assay and generate figure panels
- 03_genetic_backgrounds_screen: Analyze CRISPRi data from the clinical isolate assay and generate figure panels
- 04_growth_curves_qPCR: Analyze growth curve data and qPCR data
- 05_WGS: Process reads and call variants for WGS of reference strain and clinical isolates
