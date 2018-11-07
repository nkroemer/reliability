function out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind)

if sim2mean ==  1
       mean_temp = [];
       mean_temp = mean(TempFourD2,4);
end;

for i = 1:nr_subj
    for j = 1:nr_subj
        temp_nii_1 = [];
        temp_nii_2 = [];
        temp_nii_1 = TempFourD1(:,:,:,i);
        temp_nii_2 = TempFourD2(:,:,:,j);
        if use_roi == 1
            temp_nii_1(~r_roi_ind) = 0;
            temp_nii_2(~r_roi_ind) = 0;
        end;
        temp_1 = temp_nii_1(~isnan(temp_nii_1));
        temp_2 = temp_nii_2(~isnan(temp_nii_2));

        [out.r, out.p] = corrcoef([temp_1,temp_2]);
        out.r_mat(i,j) = out.r(1,2);
        out.p_mat(i,j) = out.p(1,2);   
    end;
    if sim2mean ==  1
       [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
       out.r_mat(i,j+1) = out.r(1,2);
       out.p_mat(i,j+1) = out.p(1,2);      
    end;
end;