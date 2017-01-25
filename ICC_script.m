% calculation of voxel-wise ICC on the basis of 4D images (FourD)
clear all;
dir_results = 'M:\SeSyN\019\Juliane\Toolbox\results_gui';
cd(dir_results);

runs = 3;
nr_subj=126;
for i = 1:runs
    if i == 1
        img_1 = load_nii(sprintf('%d_fourNEU.nii',i));
        img_1 = img_1.img;
        img_1 (~img_1) = nan;
    elseif i == 2
        img_2 = load_nii(sprintf('%d_fourNEU.nii',i));
        img_2 = img_2.img;
        img_2 (~img_2) = nan;
    elseif i == 3
        img_3 = load_nii(sprintf('%d_fourNEU.nii',i));
        img_3 = img_3.img;
        img_3 (~img_3) = nan;
    end;
end;
    
temp_img = 'M:\SeSyN\019\Juliane\Toolbox\results_gui\con_0003.nii';
temp_img=load_nii(temp_img);
dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);

if runs > 1
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                for ind_subj = 1:nr_subj
                    img1_voxel(ind_subj,1) = img_1 (ind_x, ind_y, ind_z, ind_subj);
                    img2_voxel(ind_subj,1) = img_2 (ind_x, ind_y, ind_z, ind_subj);                
                end;              
                %ICC
                if all(~isnan(img1_voxel)) || all(~isnan(img2_voxel))
                    %create data table for ANOVA
                    t = table(img1_voxel, img2_voxel,'VariableNames',{'t1','t2'});
                    %create within-subject variable
                    Time = [1 2];
                    %fit model
                    rm = fitrm(t,'t1-t2~1','WithinDesign',Time);
                    %run ANOVA
                    ranovatbl = ranova(rm);
                    %calculate ICC
                        within = ranovatbl.MeanSq(1);
                        between = ranovatbl.MeanSq(2);
                        ICC = within./(sum(between+within)); 
                else 
                    ICC = 0;
                end;
                
                %matrix with ICC value for each voxel
                ICC_1_2(ind_x, ind_y, ind_z) = ICC;
                                
            end;
        end;
    end;
end;

% save ICC map
target_img = temp_img;
target_img.fileprefix = 'ICC_1_2.nii';
target_img.img = ICC_1_2;
save_nii(target_img,target_img.fileprefix);

if runs >2

for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                   for ind_subj = 1:nr_subj
                    img2_voxel(ind_subj,1) = img_2 (ind_x, ind_y, ind_z, ind_subj);
                    img3_voxel(ind_subj,1) = img_3 (ind_x, ind_y, ind_z, ind_subj);                
                end;              
                %ICC
                if all(~isnan(img2_voxel)) || all(~isnan(img3_voxel))
                    %create data table for ANOVA
                    t = table(img2_voxel, img3_voxel,'VariableNames',{'t1','t2'});
                    %create within-subject variable
                    Time = [1 2];
                    %fit model
                    rm = fitrm(t,'t1-t2~1','WithinDesign',Time);
                    %run ANOVA
                    ranovatbl = ranova(rm);
                    %calculate ICC
                        within = ranovatbl.MeanSq(1);
                        between = ranovatbl.MeanSq(2);
                        ICC = within./(sum(between+within)); 
                else 
                    ICC = 0;
                end;
                
                %matrix with ICC value for each voxel
                ICC_2_3(ind_x, ind_y, ind_z) = ICC;
                          
            end;
        end;
    end;
 % save ICC map
target_img = temp_img;
target_img.fileprefix = 'ICC_2_3.nii';
target_img.img = ICC_2_3;
save_nii(target_img,target_img.fileprefix);   
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                   for ind_subj = 1:nr_subj
                    img1_voxel(ind_subj,1) = img_1 (ind_x, ind_y, ind_z, ind_subj);
                    img3_voxel(ind_subj,1) = img_3 (ind_x, ind_y, ind_z, ind_subj);                
                end;              
                %ICC
                if all(~isnan(img1_voxel)) || all(~isnan(img3_voxel))
                    %create data table for ANOVA
                    t = table(img1_voxel, img3_voxel,'VariableNames',{'t1','t2'});
                    %create within-subject variable
                    Time = [1 2];
                    %fit model
                    rm = fitrm(t,'t1-t2~1','WithinDesign',Time);
                    %run ANOVA
                    ranovatbl = ranova(rm);
                    %calculate ICC
                        within = ranovatbl.MeanSq(1);
                        between = ranovatbl.MeanSq(2);
                        ICC = within./(sum(between+within)); 
                else 
                    ICC = 0;
                end;
                
                %matrix with ICC value for each voxel
                ICC_1_3(ind_x, ind_y, ind_z) = ICC;
                          
                                
            end;
        end;
    end;
% save ICC map
target_img = temp_img;
target_img.fileprefix = 'ICC_1_3.nii';
target_img.img = ICC_1_3;
save_nii(target_img,target_img.fileprefix);
end;