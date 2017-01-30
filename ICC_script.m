% calculation of voxel-wise ICC on the basis of 4D images (FourD)
clear all;
dir_results = 'M:\SeSyN\019\Juliane\Toolbox\results_gui';
cd(dir_results);

runs = 3;
nr_subj=126;
for i = 1:runs
    img = load_nii(sprintf('%d_fourNEU.nii',i));
    evalstr = sprintf('img_%d = img;',i);
    eval(evalstr);
    evalstr = sprintf('img_%d = img_%d.img;',i,i);
    eval(evalstr);
end;
    
temp_img = 'M:\SeSyN\019\Juliane\Toolbox\results_gui\con_0003.nii';
temp_img=load_nii(temp_img);
dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);

data=zeros(nr_subj,runs);
ICC_con=zeros(x,y,z);
ICC_abs=zeros(x,y,z);

disp('calculating ICCs')
    for ind_x = 1:x
        for ind_y = 1:y
           for ind_z = 1:z
                for ind_run = 1:runs
                    for ind_subj = 1:nr_subj
                       estr = sprintf('img%d_voxel(ind_subj,1)= img_%d (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run);
                       eval(estr);
                    end;
                    estr = sprintf('data(:,ind_run)=img%d_voxel;',ind_run);
                    eval(estr);
                end;
                    %ICC
                    nsamples=nr_subj*runs;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:runs,
                           grandmean= grandmean + data(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(runs,1);
                    for sess=1:runs
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + data(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:runs
                            subjmean(sub)=subjmean(sub) + data(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./runs;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session
                    
                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:runs
                            WMS = WMS + (data(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (data(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:runs
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= runs.*BMS./(nr_subj-1);
                    WMS= WMS./(runs-1)./nr_subj;
                    JMS= nr_subj.*JMS./(runs-1);
                    EMS= EMS./(runs-1)./(nr_subj-1); 

                    %consistency agreement  
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*EMS); 

                    %absolute agreement 
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                                                   
                    ICC_con(ind_x, ind_y, ind_z) = voxICC_con;
                    ICC_abs(ind_x, ind_y, ind_z) = voxICC_abs;
           end; 
        end;
     end;
    


    
disp('saving ICC images')
 % save ICC maps
target_img = temp_img;
target_img.fileprefix = 'ICC_con.nii';
target_img.img = ICC_con;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = 'ICC_abs.nii';
target_img.img = ICC_abs;
save_nii(target_img,target_img.fileprefix); 
    
disp('finished ICCs')