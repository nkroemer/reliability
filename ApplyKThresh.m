img_path = 'C:\Users\Nils\Downloads\cons\Subject1_mask_day1.nii';
new_img_name = 'Subject1_mask_day1_filtered';
k = 50;

%img_template = spm_vol(img_path);
%[IMG, XYZ] = spm_read_vols(test);

img_template = load_nii(img_path);

[Y,XYZ,C] = fmri_clust_filt(img_template.img,k);

img_target = img_template;

img_target.img = Y;
img_target.fileprefix = new_img_name;

save_nii(img_target,[img_target.fileprefix '.nii'])