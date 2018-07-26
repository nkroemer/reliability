# fmreli toolbox beta v0.2

## What is the fmreli toolbox?

The toolbox was designed to simplify the assessment of reliability and similarity of fMRI data across different sessions and contrasts. It incorporates common measures of global and local reliability. Moreover, it implements the assessment of reliability for cross-sectional designs with only a single run by randomly splitting the trials in half. fmreli offers a graphical user interface (GUI) and incorporates functions provided by [*SPM*](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and and the [*Nifti and ANALYZE toolbox*](https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image). A preprint of the paper detailing the use of the toolbox is available at [bioRxiv.org](https://doi.org/10.1101/215053).

## Dependencies

To use the toolbox, you need to have MATLAB and SPM12 installed (see instructions). 

## Installing the toolbox

1. Make sure that you have [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/download/) installed.
1. Download the latest version of the toolbox on [github](https://github.com/nkroemer/reliability/) and save it into your intended directory.
1. In the MATLAB 'Home' tab, click on the 'Set path' button and select 'Add with subfolders…'. Now, select the fmreli toolbox folder, save, and close. Windows users might have to open MATLAB as an administrator (by right-clicking on the MATLAB icon) to make permanent changes.
Alternatively, you may use the command line to add the path:
```
pathtool
addpath(genpath('fmreli toolbox folder'))
savepath  fmreli toolbox folder/pathdef.m
```

By default, the toolbox comes with the [CONN](https://www.nitrc.org/projects/conn/) atlas (i.e., Harvard-Oxford brain atlas + AAL cerebellum atlas). Other customized atlases can be used, if you provide an 'atlas.nii' and a corresponding 'labels.mat' file containing the atlas labels.

That's it, you are good to go.

## Basic setup

1. After installing, you can call the toolbox by entering *fmreli* in the command line.

1. The toolbox was developed to work primarily with first-level statistics in SPM12.
The required folder structure is (at the moment) as follows:
    1. Project folder containing one subfolder for every participant
    1. Within every participant folder, there should be one folder for every session. The session should be given in the folder name and contain the first-level SPM stats including all outputs.
    For example: C:\user\projects\fmreli\subj_123\new_paradigm_1

1. Next, you can specify the design:
    1. Click on the 'Design' dialog.
    1. Enter the number of subjects.
    1. Click on the 'load subject list' dialog and load a \*.mat file containing the subject IDs as a (N-by-1) vector.
    1. Specify the number of sessions.
    1. Specify the directory where your first-level outputs are stored.
    1. Pick a prefix to be added to all output files and choose a directory where output files from the toolbox will be saved.
    1. Save your study design and close the dialog box. A \*.mat file containing the information on your study design will be saved to your output folder.

1. Now that you have specified the design, you can open it as a template and define the contrasts of interest for the following analysis.
    1. Click on the 'Contrast(s) of interest' dialog.
    1. Define the contrast(s) by the name they were given in the first-level SPM contrast manager. In case you renamed the output files, you can adjust the prefix accordingly in the dialog on the left.
    1. Click on the 'Save contrast definition' dialog. This saves a \*.mat file to your fmreli output folder containing your contrast information.

### Split-half reliability: Assessing reliability within a single contrast and session

If you have only data from one fMRI run/session and are interested in the reliability of one contrast/condition, the toolbox offers the option to split the data.

1. Click on the 'Split-Half' dialog.
1. Load your study design.
1. Indicate the number of the condition you want to split (according to the design information in the SPM.mat).
1. Define a name for the split-half output file.
1. Always check the 'split data' box throughout the use in different modules in analyzing split-half data.
1. Click 'RUN' to start.

### Assessment of global brain reliability

#### Whole-brain/ROI similarity of brain response patterns

The similarity module computes the global similarity of activation maps between participants across runs or sessions.

1. Define the study design and contrasts of interest.
1. OPTIONAL Select 'use ROI' and define the ROI name and directory if you want to restrict the analysis to a subset of voxels.
1. Click 'RUN'.

The output from this module includes heatmaps, density plots, and histograms depicting the similarity between participants across runs or sessions. Moreover, a \*.mat file will be saved containing the correlation coefficients and corresponding p-values.

#### Whole-brain/ROI overlap

In this module, the Dice and Jaccard coefficients of overlap of activated voxels can be computed. 

1. Define the study design and contrasts of interest.
1. Set a p-value used as threshold for significant activation of a voxel. By default, it is applied to the first-level statistics, but you can choose to apply the threshold to the group-level statistics as well.
1. OPTIONAL Select 'use ROI' and define the ROI name and directory if you want to restrict the analysis to a subset of voxels.
1. Click 'RUN'.

In the output folder, you will find a file containing the coefficients for every subject.

### Voxel-wise/ROI reliability

This module is used to calculate the correlation coefficients (intra-class coefficient, Pearson & Spearman correlation coefficients) on the voxel level within subjects and between contrasts.

1. Define the study design and contrasts of interest.
1. Check the boxes to select coefficients you want to analyze.
1. OPTIONAL Select 'use ROI' and define the ROI name and directory if you want to restrict the analysis to a subset of voxels.
1. Click 'RUN'.

This module will create several 3D-nifti output files named after the coefficient you had selected.

### Summary functions

#### Atlas-based reliability

With this module, you can create an output aggregated for ROIs provided by the atlas. The output is a matrix containing the reliability coefficients for every ROI in the atlas for every analysis.

## Bug report, suggestions & questions

Please do not hesitate to email us (<nils.kroemer@uni-tuebingen.de> & <juliane.froehner@tu-dresden.de>) regarding questions, feature requests or bug reports. We are happy to receive your feedback to improve fmreli. In case something does not work or you find a bug or error in the toolbox, please let us know via the [issues page](https://github.com/nkroemer/reliability/issues).
