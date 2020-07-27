# HCA Adult Heart Project


This is the working repository for the processing and analysis of the HCA Heart data in collaboration with:

**Teichmann Lab** - Wellcome Sanger Institute - Cambridge, United Kingdom.
**Hubner Lab** - MDC - Berlin, Germany.
**Noseda Lab** - Imperial College London - London, United Kingdom.
**Seidman Labs** - Harvard Medical School - Boston, MA, United States. 

The repository will host constributions from developers and analysts at Sanger as well as the collaborators in London, Berlin and Boston. 

To have a working copy of the repository, just clone it on your working machine as follow:

```
git clone --recursive https://github.com/cartal/HCA_Heart.git
```

In order to keep visible the contributions from each collaborator, any change, correction or bug fix should be submitted as 
a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) by creating a new branch from the master repo.

For more details please contact Carlos Talavera-López (ct5@sanger.ac.uk) or Monika Litviňuková (monika.litvinukova@mdc-berlin.de). 

## Repo structure

There are three main folders where files for the analysis have been stored:

1. **Archived**: Which includes the scripts to run cellbender, as well as the notebooks for downsteam analyses - **This was not used for the final analysis**.
2. **Analysis**: Contain the notebokks used in the analysis of the paper.
3. **Mapping**: Contains the scripts for the mapping of the reads to reference GRCh38 using `CellRanger-v3.0.2`. 

