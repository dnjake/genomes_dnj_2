## Statistics of Human Genetic Diversity

The project in the genomes_dnj repository characterized series of single nucleotide polymorphisms (SNPs) located in a specific region of chromosome 2.  The work in this repository is intended to extend that study to the whole human genome.  The initial focus is on the development of a few simple statistics that can characterize genetic diversity in the human genome.

## Notebooks

The results of the initial work are presented in a single jupyter notebook chrom_statistics_overview.ipynb located in the notebook_statistics folder.  An html version of that notebook is provided in the root directory.  An online version of the [Statistics of Human Genetic Diversity](https://anaconda.org/dnjake/chrom_statistics_overview/notebook) notebook is available.

## Data
There is a dependency on over 7 gigabytes of hdf5 files of thousand genome phase 3 data.  An autosome_snp_data package that includes the needed data files is available in the [google drive genomes_dnj folder](https://drive.google.com/drive/folders/0B0N6jrX3WLdfTjBialkzRXNzSkE?usp=sharing).  An outline of the process for assembling subpackages into a package executable is provided in a blog post [Assembling Genomes_dnj Packages](https://genomesdnj.blogspot.com/2017/11/assembling-genomesdnj-packages.html).

## Licence

The python code and jupyter notebook are all the work of the owner of this repository.  Reuse of that work is encouraged without any constraint.  The input data all comes from the 1000 genomes project phase 3.  Its reuse is subject to any constraints imposed by that project.
