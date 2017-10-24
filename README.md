# fMRelI Toolbox beta

## What is the fMRelI toolbox?

The toolbox was designed to simplify the assessment of reliability and similarity of fMRI data across different sessions and contrasts. It incorporates common measures of global and local reliability and also offers the assessment of reliability within a single session and contrast by splitting the data available in half. fMRelI provides a graphical user interface and is based on *SPM* functions and functions from the *Nifti and ANALYZE toolbox*. A preprint of a paper detailing an example of how the toolbox can be employed is available at [bioArXiv.org](paste link here).

## Dependencies

In order to use the toolbox you need to have installed Matlab and SPM12. It currently comes with the Harvard-Oxford brain atlas.
The use of other brain atlases is possible, all you will need is a 'labels.mat' file containing the atlas labels.

## Installing the toolbox

1. Make sure that you have SPM12 installed.
1. Download the latest Version of the toolbox on [github](https://github.com/nkroemer/reliability/tree/fMRelI_beta0.1) and save it into your desired Matlab toolbox directory.
1. If you want to use an additional brain atlas, the files 'atlas.nii', 'atlas.txt' and 'labels.mat' have to be in a subfolder called 'atlas'. 
1. If you like structure, you can add 'create_xSPM.m', 'fmri_clust_filt.m', 'id.mat', 'load_nii.m', 'save_nii.m' and 'template_3D-4D.mat' to a folder called 'scripts_templates'. 
1. In the Matlab 'Home' tab, click on the 'Set path' button and in the dialog box select 'Add with subfolders…'. Now, select the fMRelI toolbox folder, save and close. Windows users might have to open Matlab as an administrator (by right-clicking on the Matlab Symbol) in order for the changes to be permanent.
Alternatively, you may enter this in the command line:
```
pathtool
addpath(genpath('fMRelI toolbox folder'))
savepath  fMRelI toolbox folder/pathdef.m
```

## Basic setup

1. After installing, you can call the toolbox by entering *fMRelI* in the command line.

In order to use the toolbox you need to have run the first-level statistics in SPM12 on all the contrasts you are interested in

1. Start by specifying your design:
    1. Click on the 'Design' dialog
    1. Enter number of subjects
    1. Click on the 'load subject list' dialog and load a .mat file containing the subject IDs as a (N-by-1) vector
    1. Specify number of sessions
    1. Specify the directory where your first-level outputs are stored
    1. Pick a prefix to be added to all output files and choose a directory where all the output files from the toolbox will be saved 
    1. Save your study design and close the current dialog box.

1. Next, you can open a study design which was created earlier and define the contrasts you are interested in comapring.
    1. Click on the 'Contrast(s) of interest' dialog

### Split-half reliability: Assessing reliability within a single contrast and session

If you have only data from one fMRI session and are interested in the reliability of one contrast/condition the toolbox offers a function for splitting the data.

1. Click on the 'Split-Half' dialog
1.

### Assessment of global brain reliability

### Assessment of voxel-wise reliability



