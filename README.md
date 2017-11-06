# fMRelI Toolbox beta 1.0

## What is the fMRelI toolbox?

The toolbox was designed to simplify the assessment of reliability and similarity of fMRI data across different sessions and contrasts. It incorporates common measures of global and local reliability and also offers the assessment of reliability within a single session and contrast by splitting the data available in half. fMRelI provides a graphical user interface and is based on *SPM* functions and functions from the *Nifti and ANALYZE toolbox*. A preprint of a paper detailing an example of how the toolbox can be employed is available at [bioArXiv.org](paste link here).

## Dependencies

In order to use the toolbox you need to have installed Matlab and SPM12. It currently comes with the Harvard-Oxford brain atlas.
The use of other brain atlases is possible, all you will need is a 'labels.mat' file containing the atlas labels.

## Installing the toolbox

1. Make sure that you have SPM12 installed.
1. Download the latest Version of the toolbox on [github](https://github.com/nkroemer/reliability/tree/fMRelI_beta0.1) and save it into your desired Matlab toolbox directory.
1. If you want to use an additional brain atlas, the files 'atlas.nii', 'atlas.txt' and 'labels.mat' have to be saved to a subfolder called 'atlas'. 
1. In the Matlab 'Home' tab, click on the 'Set path' button and in the dialog box select 'Add with subfolders…'. Now, select the fMRelI toolbox folder, save and close. Windows users might have to open Matlab as an administrator (by right-clicking on the Matlab Symbol) in order for the changes to be permanent.
Alternatively, you may enter this in the command line:
```
pathtool
addpath(genpath('fMRelI toolbox folder'))
savepath  fMRelI toolbox folder/pathdef.m
```

## Basic setup

1. After installing, you can call the toolbox by entering *fmreli* in the command line.

1. In order to use the toolbox you need to have run the first-level statistics in SPM12 on all the contrasts you are interested in.
The required folder structure is as follows:
    1. Master folder containing one subfolder for every participant
    1. Within every participant folder, there should be one folder for every timepoint/session. The timepoint should be defined in the folder name. This folder should contain all of the SPM stats outputs for all contrasts.
    For example: C:\user\projects\fMRelI\subject123\new_paradigm_1

1. Start specifying your design:
    1. Click on the 'Design' dialog.
    1. Enter number of subjects.
    1. Click on the 'load subject list' dialog and load a .mat file containing the subject IDs as a (N-by-1) vector.
    1. Specify number of sessions.
    1. Specify the directory where your first-level outputs are stored.
    1. Pick a prefix to be added to all output files and choose a directory where all the output files from the toolbox will be saved .
    1. Save your study design and close the current dialog box. A .mat file containing the information on your study design will be saved to your output folder.

1. Next, you can open a study design which was created earlier and define the contrasts you are interested in comapring.
    1. Click on the 'Contrast(s) of interest' dialog.
    1. Define the contrast(s) by the name they were given in the first level SPM contrast manager. In case you renamed the output files, you can adjust the prefix accordingly in the dialog on the left.
    1. Click on the 'Save contrast definition' dialog. This saves a mat file to your fMRelI output folder containing your contrast information.

### Split-half reliability: Assessing reliability within a single contrast and session

If you have only data from one fMRI session and are interested in the reliability of one contrast/condition the toolbox offers a function for splitting the data.

1. Click on the 'Split-Half' dialog.
1. Load your study design.
1. Define the data by inserting the onset regressor for your contrast of interest.
1. Define a name for the split data output file.
1. Remember to always check the 'splitted data' box in the different analysis modules when analysing split data.
1. Click 'RUN' to start.

### Assessment of global brain reliability

#### Assessment of similarity

The similarity module allows to compute the global similarity of activation maps between participants over all sessions.

1. Define the study design and contrast definition files.
1. If you are interested in similarity of activation in a single ROI, select 'use ROI' and define the ROI name and directory.
1. Click 'RUN'.

The output from this module includes heatmaps, density plots and histograms depicting the simmilarity between participants and sessions. Moreover there is a .mat file containing the underlying global correlation and p-values.

#### Assessment of overlap

In this module, the Dice and Jaccard coefficients of overlap of activated voxels can be computed. 

1. Define the study design and contrast definition files.
1. Set a p-value which the defines the threshold of significant activation in a voxel. You can choose to apply this threshold to the group-level statistics, by default it is applied to the subject-level statistics.
1. If you are interested in overlap of activation in a single ROI, select 'use ROI' and define the ROI name and directory.
1. Click 'RUN'.

In the output file, you will find the coefficients for every single subject.

### Assessment of voxel-wise reliability

This module is used to calculate the correlation coefficients (Intra-class coefficient, Pearson & Spearman correlation coefficients) on the voxel- or ROI level within subjects and between contrasts.

1. Define the study design and contrast definition files.
1. Check the boxes for the coefficients you are interested in.
1. If you are interested in reliability of activation in a single ROI, select 'use ROI' and define the ROI name and directory.
1. Click 'RUN'.

It will create several 3D-nifti output files named after the coefficient you chose.

### Summary functions

#### Atlas-based reliability

With this module, you can create an output similar to that from the voxel-wise reliability module, but containing information about atlas-defined ROIs instead of the voxel-level. The output is a matrix containing the correlation coefficients for every ROI in the atlas for every comparison, respectively.
