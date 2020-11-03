# PhyloPhlAn

PhyloPhlAn is an integrated pipeline for large-scale phylogenetic profiling of genomes and metagenomes.
PhyloPhlAn is an accurate, rapid, and easy-to-use method for large-scale microbial genome characterization and phylogenetic analysis at multiple levels of resolution.
PhyloPhlAn can assign both genomes and metagenome-assembled genomes (MAGs) to species-level genome bins (SGBs).
PhyloPhlAn can reconstruct strain-level phylogenies using clade-specific maximally informative phylogenetic markers, and can also scale to very-large phylogenies comprising >17,000 microbial species.

# Installation

## Bioconda

You can install PhyloPhlAn using conda as follows:

~~~Bash
conda install -c bioconda phylophlan
~~~


## Repository

You can clone the PhyloPhlAn repository from GitHub:

~~~Bash
git clone https://github.com/biobakery/phylophlan
cd phylophlan
python setup.py install
~~~

Then remember to check that the [Dependencies and Tools](https://github.com/biobakery/phylophlan/wiki#requirements) are in stalled and available in your system.


# Tutorials

* [PhyloPhlAn User manual](https://github.com/biobakery/phylophlan/wiki)
* [PhyloPhlAn Tutorials](https://github.com/biobakery/biobakery/wiki/PhyloPhlAn3)
* [PhyloPhlAn Example 1: Phylogenetic characterization of isolate genomes of a given species (S. aureus)](https://github.com/biobakery/biobakery/wiki/PhyloPhlAn-3.0:-Example-01:-S.-aureus)
* [PhyloPhlAn Example 2: Prokaryotes Tree of life reconstruction](https://github.com/biobakery/biobakery/wiki/PhyloPhlAn-3.0:-Example-02:-Tree-of-life)
* [PhyloPhlAn Example 3: Metagenomic analysis of the Ethiopian cohort](https://github.com/biobakery/biobakery/wiki/PhyloPhlAn-3.0:-Example-03:-Metagenomic-application)
* [PhyloPhlAn Example 4: High-resolution phylogeny of genomes and MAGs of a known species (E. coli)](https://github.com/biobakery/biobakery/wiki/PhyloPhlAn-3.0:-Example-04:-E.-coli)
* [PhyloPhlAn Example 5: Phylogenetically characterization of an unknown SGB from the Proteobacteria phylum](https://github.com/biobakery/biobakery/wiki/PhyloPhlAn-3.0:-Example-05:-Proteobacteria)


# Support

We provide support through [the bioBakery help forum](https://forum.biobakery.org/) and through the issues tracking system of the [PhyloPhlAn repository](https://github.com/biobakery/phylophlan/issues).


# Citation

If you used PhyloPhlAn please cite the following paper:

**Precise phylogenetic analysis of microbial isolates and genomes from metagenomes using PhyloPhlAn 3.0**  
_Francesco Asnicar_, _Andrew Maltez Thomas_, _Francesco Beghini_, _Claudia Mengoni_, _Serena Manara_, _Paolo Manghi_, _Qiyun Zhu_, _Mattia Bolzan_, _Fabio Cumbo_, _Uyen May_, _Jon G. Sanders_, _Moreno Zolfo_, _Evguenia Kopylova_, _Edoardo Pasolli_, _Rob Knight_, _Siavash Mirarab_, _Curtis Huttenhower_, and _Nicola Segata_  
Nat Commun 11, 2500 (2020)  
DOI: [https://doi.org/10.1038/s41467-020-16366-7](https://doi.org/10.1038/s41467-020-16366-7)
