# HCA Heart


This is the working repository for the processing and analysis of the HCA Heart data in collaboration with:

1. **Teichmann Lab** - Wellcome Sanger Institute - Cambridge, United Kingdom.
2. **Hubner Lab** - MDC - Berlin, Germany.
3. **Noseda Lab** - Imperial College London - London, United Kingdom.
4. **Seidman Labs** - Harvard Medical School - Boston, MA, United States. 

The repository will host constributions from developers at Sanger as well as the collaborators in London, Berlin and Boston. 

To have a working copy of the repository, just clone it on your working machine as follow:

```
git clone --recursive https://github.com/cartal/HCA_Heart.git
```

In order to keep visible the contributions from each collaborator, any change, correction or bug fix should be submitted as 
a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) by creating a new branch from the master repo.

**Please keep in mind that the code here can't be shared or used outside this consortium without the explicit authorisation from their respective developers.** 

For more details please contact Carlos Talavera-López (ct5@sanger.ac.uk) or Monika Litviňuková (monika.litvinukova@mdc-berlin.de). 

## The analysis

There are two main folders where the notebook for the analyses can be found:

1. **CellBender-based**: Which includes the scripts to run cellbender, as well as the notebooks for downsteam analyses.
2. **CellRanger-based**: Which uses the `filtered_feature_bc_matrix` as main input. No background removal step is performed here.


