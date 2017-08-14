function varargout = CorrICC(varargin)
% CORRICC MATLAB code for CorrICC.fig
%      CORRICC, by itself, creates a new CORRICC or raises the existing
%      singleton*.
%
%      H = CORRICC returns the handle to a new CORRICC or the handle to
%      the existing singleton*.
%
%      CORRICC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRICC.M with the given input arguments.
%
%      CORRICC('Property','Value',...) creates a new CORRICC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CorrICC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CorrICC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CorrICC

% Last Modified by GUIDE v2.5 17-Jul-2017 17:13:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrICC_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrICC_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CorrICC is made visible.
function CorrICC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CorrICC (see VARARGIN)

% Choose default command line output for CorrICC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CorrICC wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CorrICC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in design.
function design_Callback(hObject, eventdata, handles)
% hObject    handle to design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
study_design = cellstr(spm_select(1,'mat','load study design'));
load(study_design{1});
assignin('base','study_design',study_design);

% --- Executes on button press in contrast.
function contrast_Callback(hObject, eventdata, handles)
% hObject    handle to contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

% --- Executes on button press in pear.
function pear_Callback(hObject, eventdata, handles)
% hObject    handle to pear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pear


% --- Executes on button press in spea.
function spea_Callback(hObject, eventdata, handles)
% hObject    handle to spea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spea


% --- Executes on button press in cons.
function cons_Callback(hObject, eventdata, handles)
% hObject    handle to cons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cons


% --- Executes on button press in abs.
function abs_Callback(hObject, eventdata, handles)
% hObject    handle to abs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of abs


% --- Executes on button press in roi.
function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi


% --- Executes on button press in roi_dir.
function roi_dir_Callback(hObject, eventdata, handles)
% hObject    handle to roi_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi_dir = spm_select(1,'dir','choose roi directory');
assignin('base','roi_dir',roi_dir);



function roi_name_Callback(hObject, eventdata, handles)
% hObject    handle to roi_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi_name as text
%        str2double(get(hObject,'String')) returns contents of roi_name as a double


% --- Executes during object creation, after setting all properties.
function roi_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in split.
function split_Callback(hObject, eventdata, handles)
% hObject    handle to split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of split

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Starting calculation of correlation measures...');
%% define file seperator 
f = filesep;

%% set parameters
%get information out of study_design
box=pwd;
study_design=evalin('base','study_design');
contrast_def=evalin('base','contrast_def');
runs=str2double(study_design.number_sessions); 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
stats_dir=study_design.stats_directory;
stats_path=study_design.stats_path;
nr_para = study_design.number_parametric;

dir_results = study_design.results_directory;


if runs == 1
    single_run = str2double(study_design.identifier_session);
end;

% get GUI input
roi = get(handles.roi,'value');
if roi == 1
    roi_name=get(handles.roi_name,'String');
    roi_dir=evalin('base','roi_dir');
end;
split = get(handles.split,'value');
cons = get(handles.cons,'value');
abs = get(handles.abs,'value');
pear = get(handles.pear,'value');
spea = get(handles.spea,'value');

%load contrast information
two_cons = contrast_def.two_contrasts;
if two_cons == 0
    con=contrast_def.contrast;
else
    con=contrast_def.contrast1;
    con1=contrast_def.contrast1;
    con2=contrast_def.contrast2;
    con1_count=contrast_def.number_contrast1;
    con2_count=contrast_def.number_contrast2;
end;
cd(dir_results);
%% ROI settings
if roi == 1
    % load ROI
    cd(roi_dir);
    roi_ful = dir(sprintf('%s*',roi_name));
    if isstruct(roi_ful)
        if length(roi_ful)==2
            roi_ful = roi_ful(2).name;
        else
            roi_ful = roi_ful(1).name;
        end;
    end;
    compl = sprintf('%s%s%s%s',roi_dir,f,f, roi_ful);

    % reslice ROI
    disp('...reslicing ROI...');
    stats_filled = sprintf(stats_dir,1);
    temp = sprintf('%s%s%s%s%s%s%s%s%s%s,1',stats_path,f,f,vp{1},f,f,stats_filled,f,f,con);
    matlabbatch{1}.spm.spatial.coreg.write.ref = {temp};
    temp_1 = sprintf('%s,1',compl);
    matlabbatch{1}.spm.spatial.coreg.write.source = {temp_1};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

    % save batch
    save('reslice_roi','matlabbatch');

    % run batch
    spm_jobman('serial',matlabbatch);

    % create index for ROI voxels
    r_roi=dir(sprintf('r%s*',roi_name));
    if length(r_roi)==2
        compl1 = sprintf('%s%s%sr%s.img',roi_dir,f,f,roi_name);
        movefile(compl1,dir_results,'f');
        compl2 = sprintf('%s%s%sr%s.hdr',roi_dir,f,f,roi_name);
        movefile(compl2,dir_results,'f');
        cd(dir_results);    
        r_roi = load_nii(sprintf('r%s.img',roi_name));
        r_roi_ind = r_roi.img==1;
    else
        compl = sprintf('%s%s%sr%s.nii',roi_dir,f,f,roi);
        if ~strcmp(roi_dir,dir_results)
        movefile(compl,dir_results,'f');
        end;
        cd(dir_results);    
        r_roi = load_nii(sprintf('r%s.nii',roi));
        r_roi_ind = r_roi.img==1;        
    end;
    cd(dir_results);

    % load images and set voxels outside ROI = 0
    disp('...loading images and mask with ROI...');
    if split == 0 && two_cons == 0
        for i = 1:runs
            if runs == 1
                i = single_run;
            end;
            img = load_nii(sprintf('4D_%d.nii',i));
            eval(sprintf('img_%d = img;',i));
            eval(sprintf('img_%d = img_%d.img;',i,i));
            for count_subj = 1:nr_subj
                eval(sprintf('img_temp = img_%d(:,:,:,count_subj);',i));
                img_temp(~r_roi_ind) = 0;
                eval(sprintf('img_%d(:,:,:,count_subj) = img_temp;',i));
            end;
        end;
        disp('...loads image dimensions..');
        stats_temp =sprintf(stats_dir,1);
        temp_img = sprintf('%s%s%s%s%s%s%s%s%s%s',stats_path,f,f,vp{1},f,f,stats_temp,f,f,con);
        temp_img=load_nii(temp_img);
        dim = size(temp_img.img);
        x = dim(1);
        y = dim(2);
        z = dim(3);
        if nr_para > 0 
            for ind_para = 1:nr_para
               for i = 1:runs
                    if runs == 1
                        i = single_run;
                    end;
                    img = load_nii(sprintf('4D_par%d_%d.nii',ind_para,i));
                    eval(sprintf('img_%d_par = img;',i));
                    eval(sprintf('img_%d_par = img_%d_par.img;',i,i));
                    for count_subj = 1:nr_subj
                        eval(sprintf('img_temp = img_%d_par(:,:,:,count_subj);',i));
                        img_temp(~r_roi_ind) = 0;
                        eval(sprintf('img_par%d_%d(:,:,:,count_subj) = img_temp;',ind_para,i));
                    end;
                end;
            end;
        end;

    elseif split == 1

        for i = 1:runs 
            if runs == 1
                i = single_run;
            end;

            img1 = load_nii(sprintf('4D_split1_%d.nii',i));
            eval(sprintf('img_%d_split1 = img1;',i));
            eval(sprintf('img_%d_split1 = img_%d_split1.img;',i,i));
            for count_subj = 1:nr_subj
                eval(sprintf('img_temp = img_%d_split1(:,:,:,count_subj);',i));
                img_temp(~r_roi_ind) = 0;
                eval(sprintf('img_%d_split1(:,:,:,count_subj) = img_temp;',i));
            end;

            img2 = load_nii(sprintf('4D_split2_%d.nii',i));
            eval(sprintf('img_%d_split2 = img2;',i));
            eval(sprintf('img_%d_split2 = img_%d_split2.img;',i,i));
            for count_subj = 1:nr_subj
                eval(sprintf('img_temp = img_%d_split2(:,:,:,count_subj);',i));
                img_temp(~r_roi_ind) = 0;
                eval(sprintf('img_%d_split2(:,:,:,count_subj) = img_temp;',i));
            end;
        end;

        if nr_para > 0
            for ind_para = 1:nr_para
                   for i = 1:runs    
                       if runs == 1
                           i = single_run;
                       end;
                        nii = sprintf('4D_split1_par%d_%d.nii',ind_para,i);
                        img1 = load_nii(nii);
                        eval(sprintf('img_%d_par%d_split1 = img1;',i,ind_para));
                        eval(sprintf('img_%d_par%d_split1 = img_%d_par%d_split1.img;',i,ind_para,i,ind_para));
                        for count_subj = 1:nr_subj
                            eval(sprintf('img_temp = img_%d_par%d_split1(:,:,:,count_subj);',i,ind_para));
                            img_temp(~r_roi_ind) = 0;
                            eval(sprintf('img_%d_par%d_split1(:,:,:,count_subj) = img_temp;',i,ind_para));
                        end;                   

                        nii2 = sprintf('4D_split2_par%d_%d.nii',ind_para,i);
                        img2 = load_nii(nii2);
                        eval(sprintf('img_%d_par%d_split2 = img2;',i,ind_para));
                        eval(sprintf('img_%d_par%d_split2 = img_%d_par%d_split2.img;',i,ind_para,i,ind_para));
                        for count_subj = 1:nr_subj
                            eval(sprintf('img_temp = img_%d_par%d_split2(:,:,:,count_subj);',i,ind_para));
                            img_temp(~r_roi_ind) = 0;
                            eval(sprintf('img_%d_par%d_split2(:,:,:,count_subj) = img_temp;',i,ind_para));
                        end;   
                    end; 
            end;
        end;
        disp('...loads image dimensions..');
        stats_temp =sprintf(stats_dir,1);
        temp_img = sprintf('%s%s%s%s%s%s%s%s%s%s',stats_path,f,f,vp{1},f,f,stats_temp,f,f,con);
        temp_img=load_nii(temp_img);
        dim = size(temp_img.img);
        x = dim(1);
        y = dim(2);
        z = dim(3);
    elseif two_cons == 1
        for i = 1:runs    
            if runs == 1
                i = single_run;
            end;
            img1 = load_nii(sprintf('4D_%s_%d.nii',con1,i));
            eval(sprintf('img_%d_con1 = img1;',i));
            eval(sprintf('img_%d_con1 = img_%d_con1.img;',i,i));
            for count_subj = 1:nr_subj
                eval(sprintf('img_temp = img_%d_con1(:,:,:,count_subj);',i));
                img_temp(~r_roi_ind) = 0;
                eval(sprintf('img_%d_con1(:,:,:,count_subj) = img_temp;',i));
            end;

            img2 = load_nii(sprintf('4D_%s_%d.nii',con2,i));
            eval(sprintf('img_%d_con2 = img2;',i));
            eval(sprintf('img_%d_con2 = img_%d_con2.img;',i,i));
            for count_subj = 1:nr_subj
                eval(sprintf('img_temp = img_%d_con2(:,:,:,count_subj);',i));
                img_temp(~r_roi_ind) = 0;
                eval(sprintf('img_%d_con2(:,:,:,count_subj) = img_temp;',i));
            end;
        end;
        disp('...loads image dimensions..');
        stats_temp =sprintf(stats,1);
        temp_img = sprintf('%s%s%s%s%s%s%s%s%s%s',path,f,f,vp{1},f,f,stats_temp,f,f,con1);
        temp_img=load_nii(temp_img);
        dim = size(temp_img.img);
        x = dim(1);
        y = dim(2);
        z = dim(3);
    end;    
else
    % load 4D images whole-brain, without ROI
    disp('...loads image dimensions..');
    stats_temp =sprintf(stats_dir,1);
    temp_img = sprintf('%s%s%s%s%s%s%s%s%s%s',stats_path,f,f,vp{1},f,f,stats_temp,f,f,con);
    temp_img=load_nii(temp_img);
    dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);
    if split == 0 && two_cons == 0
        for i = 1:runs
            img = load_nii(sprintf('4D_%d.nii',i));
            eval(sprintf('img_%d = img.img;',i));
            if nr_para > 0
                for ind_para = 1:nr_para
                    img = load_nii(sprintf('4D_par%d_%d.nii',ind_para,i));
                    eval(sprintf('img_par%d_%d = img.img;',ind_para,i));
                end;
            end;
        end;
    elseif split == 1
        for i = 1:runs 
            if runs == 1
                i = single_run;
            end;
            img1 = load_nii(sprintf('4D_split1_%d.nii',i));    
            img2 = load_nii(sprintf('4D_split2_%d.nii',i));
            eval(sprintf('img_%d_split1 = img1.img;',i));
            eval(sprintf('img_%d_split2 = img2.img;',i));
            if nr_para > 0
                for ind_para = 1:nr_para
                    img1 = load_nii(sprintf('4D_split1_par%d_%d.nii',ind_para,i));
                    img1 = load_nii(sprintf('4D_split2_par%d_%d.nii',ind_para,i));
                    eval(sprintf('img_%d_par%d_split1 = img1.img;',i,ind_para));
                    eval(sprintf('img_%d_par%d_split2 = img2.img;',i,ind_para)); 
                end;
            end;
        end;                                    
    elseif two_cons == 1
        for i = 1:runs    
            if runs == 1
                i = single_run;
            end;
            img1 = load_nii(sprintf('4D_%s_%d.nii',con1,i));
            img2 = load_nii(sprintf('4D_%s_%d.nii',con2,i));
            eval(sprintf('img_%d_con1 = img1.img;',i));
            eval(sprintf('img_%d_con2 = img2.img;',i));            
        end;            
    end;
end;

%% create correlation maps 
cd (dir_results); 
if pear == 1 || spea == 1
% initiate summary
cols = {};
summary = [];
if split == 0 && two_cons == 0 && runs > 1    
    count_comp = 0;
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                count_comp = count_comp +1;
                fprintf('...creates correlation maps for session %d and session %d...\n',i_run,i_run+i_sec);
                %load 4D images
                eval(sprintf('one=img_%d;',i_run));
                one (~one) = nan;
                eval(sprintf('second=img_%d;',i_run+i_sec));
                second (~second) = nan;

                % create correlation vectors 
                r_vec_pear_1_2 = zeros(x,y,z);
                r_vec_spea_1_2 = zeros(x,y,z);
                z_r_vec_pear_1_2 = zeros(x,y,z);
                z_r_vec_spea_1_2 = zeros(x,y,z);

                for ind_x = 1:x
                     fprintf('...voxel x = %d,...\n',ind_x) 
                   for ind_y = 1:y
                        for ind_z = 1:z
                            first_voxel = one (ind_x, ind_y, ind_z, :);
                            second_voxel = second (ind_x, ind_y, ind_z, :);
                            % Pearson
                            if pear==1
                                r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                                if isnan (r(1,2))
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                else
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
                                end;
                            end;
                            %Spearman
                            if spea==1
                                first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                                second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                                r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                                r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                                z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);                    
                            end;
                        end;
                    end;
                end;

            % save correlation maps
            if pear == 1
            target_img = temp_img;
            file = sprintf('CorrMaps_pear_%d_%d.nii',i_run,i_run+i_sec);
            target_img.fileprefix = file;
            target_img.img = r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_CorrMaps_pear_%d_%d.nii',i_run,i_run+i_sec);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);      
            end;
            
            if spea == 1
            target_img = temp_img;
            file = sprintf('CorrMaps_spea_%d_%d.nii',i_run,i_run+i_sec);
            target_img.fileprefix = file;
            target_img.img = r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps_spea_%d_%d.nii',i_run,i_run+i_sec);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);                      
            end;
            
            %compute z mean and inverse Fisher's transformation
            if pear == 1
            mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_pear_%d_%d',i_run,i_run+i_sec);
            end;
            if spea == 1
            mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_spea_%d_%d',i_run,i_run+i_sec);    
            end;
            end; 
        end;
    end;
    % calculate average correlation for all comparisons
    avg_corr_4D_pear = zeros(x,y,z,count_comp);
    avg_corr_4D_spea = zeros(x,y,z,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                ind_comp = ind_comp+1;
                if pear == 1
                map = load_nii(sprintf('z_CorrMaps_pear_%d_%d.nii',i_run,i_run+i_sec));
                avg_corr_4D_pear(:,:,:,ind_comp) = map.img;
                end;
                if spea == 1
                map = load_nii(sprintf('z_CorrMaps_spea_%d_%d.nii',i_run,i_run+i_sec));
                avg_corr_4D_spea(:,:,:,ind_comp) = map.img;
                end;                
            end; 
        end;
    end;
    avg_pear = zeros(x,y,z);
    avg_spea = zeros(x,y,z);
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                if pear == 1
                    avg_temp = mean(avg_corr_4D_pear(ind_x,ind_y,ind_z,:));
                    avg_pear(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                end;
                if spea == 1
                    avg_temp = mean(avg_corr_4D_spea(ind_x,ind_y,ind_z,:));
                    avg_spea(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                end;                
            end;
        end,
    end;
    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = 'CorrMaps_pear.nii';
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = 'CorrMaps_spea.nii';
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;

if nr_para > 0    
    count_comp = 0;
    for i_par = 1:nr_para
        for i_run = 1:runs
           for i_sec = 1:runs-i_run
                if i_run+1 <= runs
                count_comp = count_comp +1;
                fprintf('...creates correlation maps for parametric modulator %d in session %d and session %d...\n',i_par,i_run,i_run+i_sec);
                %load 4D images
                eval(sprintf('one=img_par%d_%d;',i_par,i_run));
                one (~one) = nan;
                eval(sprintf('second=img_par%d_%d;',i_par,i_run+i_sec));
                second (~second) = nan;

                % create correlation vectors 
                r_vec_pear_1_2 = zeros(x,y,z);
                r_vec_spea_1_2 = zeros(x,y,z);
                z_r_vec_pear_1_2 = zeros(x,y,z);
                z_r_vec_spea_1_2 = zeros(x,y,z);

                for ind_x = 1:x
                    fprintf('...voxel x = %d,...\n',ind_x) 
                    for ind_y = 1:y
                        for ind_z = 1:z
                            first_voxel = one (ind_x, ind_y, ind_z, :);
                            second_voxel = second (ind_x, ind_y, ind_z, :);
                            % Pearson
                            if pear==1
                                r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                                if isnan (r(1,2))
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                else
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
                                end;
                            end;
                            %Spearman
                            if spea==1
                                first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                                second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                                r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                                r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                                z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);                    
                            end;
                        end;
                    end;
                end;

                % save correlation maps
                target_img = temp_img;
                file = sprintf('CorrMaps_pear_par%d_%d_%d.nii',i_par,i_run,i_run+i_sec);
                target_img.fileprefix = file;
                target_img.img = r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);  

                target_img = temp_img;
                file = sprintf('z_CorrMaps_pear_par%d_%d_%d.nii',i_par,i_run,i_run+i_sec);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);                     

                target_img = temp_img;
                file = sprintf('CorrMaps_spea_par%d_%d_%d.nii',i_par,i_run,i_run+i_sec);
                target_img.fileprefix = file;
                target_img.img = r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                target_img = temp_img;
                file = sprintf('z_CorrMaps_spea_par%d_%d_%d.nii',i_par,i_run,i_run+i_sec);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);                      

                %compute z mean and inverse Fisher's transformation
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear_par%d_%d_%d',i_par,i_run,i_run+i_sec);

                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea_par%d_%d_%d',i_par,i_run,i_run+i_sec);    

                end; 
            end;
        end;
    end;
    
    % calculate average correlation for all comparisons
    for i_par = 1:nr_para
        avg_corr_4D_pear = zeros(x,y,z,count_comp);
        avg_corr_4D_spea = zeros(x,y,z,count_comp);
        ind_comp = 0;
        for i_run = 1:runs
            for i_sec = 1:runs-i_run
                if i_run+i_sec <= runs
                    ind_comp = ind_comp+1;
                    if pear == 1
                        map = load_nii(sprintf('z_CorrMaps_pear_par%d_%d_%d.nii',i_par,i_run,i_run+i_sec));
                        avg_corr_4D_pear(:,:,:,ind_comp) = map.img;
                    end;
                    if spea == 1
                        map = load_nii(sprintf('z_CorrMaps_spea_par%d_%d_%d.nii',i_par,i_run,i_run+i_sec));
                        avg_corr_4D_spea(:,:,:,ind_comp) = map.img;
                    end;                
                end; 
            end;
        end;
        avg_pear = zeros(x,y,z);
        avg_spea = zeros(x,y,z);
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z
                    if pear == 1
                        avg_temp = mean(avg_corr_4D_pear(ind_x,ind_y,ind_z,:));
                        avg_pear(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;
                    if spea == 1
                        avg_temp = mean(avg_corr_4D_spea(ind_x,ind_y,ind_z,:));
                        avg_spea(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;                
                end;
            end,
        end;
        if pear == 1
            target_img = temp_img;
            target_img.fileprefix = sprintf('CorrMaps_pear_par%d.nii',i_par);
            target_img.img = avg_pear;
            save_nii(target_img,target_img.fileprefix); 
        end;
        if spea == 1
            target_img = temp_img;
            target_img.fileprefix = sprintf('CorrMaps_spea_par%d.nii',i_par);
            target_img.img = avg_spea;
            save_nii(target_img,target_img.fileprefix);        
        end;   
    end;
end;
%% based on split half
elseif split == 1
    for i = 1:runs
        if runs == 1
            i = single_run;
        end;
        fprintf('...creates correlation maps for splitted session %d...\n',i);
        
        %load 4D images
        eval(sprintf('one=img_%d_split1;',i));
        one (~one) = nan;
        eval(sprintf('second=img_%d_split2;',i));
        second (~second) = nan;

        % create correlation vectors 
        r_vec_pear_1_2 = zeros(x,y,z);
        r_vec_spea_1_2 = zeros(x,y,z);
        z_r_vec_pear_1_2 = zeros(x,y,z);
        z_r_vec_spea_1_2 = zeros(x,y,z);
    
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z
                    first_voxel = one (ind_x, ind_y, ind_z, :);
                    second_voxel = second (ind_x, ind_y, ind_z, :);
                    % Pearson
                    if pear==1
                        r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                        if isnan (r(1,2))
                            r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                            z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                        else
                      % correction for understimation of reliability via
                      % split-half
                            r = (2.*r(1,2))./(1+r(1,2));
                            r_vec_pear_1_2(ind_x, ind_y, ind_z) = r;
                            z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r);
                        end;
                    end;
                    %Spearman
                    if spea==1
                        first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                        second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                        r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                        % correction for understimation of reliability via
                        % split-half
                        r = (2.*r)./(1+r);
                        r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                        z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);  
                    end;
                end;
            end;
        end;
        
        % save correlation map
        target_img = temp_img;
        file = sprintf('CorrMaps_pear_%d_split.nii',i);
        target_img.fileprefix = file;
        target_img.img = r_vec_pear_1_2;
        save_nii(target_img,target_img.fileprefix);   

        target_img = temp_img;
        file = sprintf('CorrMaps_spea_%d_split.nii',i);
        target_img.fileprefix = file;
        target_img.img = r_vec_spea_1_2;
        save_nii(target_img,target_img.fileprefix);

        target_img = temp_img;
        file = sprintf('z_CorrMaps_pear_%d_split.nii',i);
        target_img.fileprefix = file;
        target_img.img = z_r_vec_pear_1_2;
        save_nii(target_img,target_img.fileprefix);   

        target_img = temp_img;
        file = sprintf('z_CorrMaps_spea_%d_split.nii',i);
        target_img.fileprefix = file;
        target_img.img = z_r_vec_spea_1_2;
        save_nii(target_img,target_img.fileprefix);
                    
        %compute z mean and inverse Fisher's transformation
        mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_pear_split_%d',i);
        mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_spea_split_%d',i);  
    end; 
    
    % calculate average correlation for all comparisons
        avg_corr_4D_pear = zeros(x,y,z,runs);
        avg_corr_4D_spea = zeros(x,y,z,runs);
        for i_run = 1:runs
            if pear == 1
                map = load_nii(sprintf('z_CorrMaps_pear_%d_split.nii',i_run));
                avg_corr_4D_pear(:,:,:,i_run) = map.img;
            end;
            if spea == 1
                map = load_nii(sprintf('z_CorrMaps_spea_%d_split.nii',i_run));
                avg_corr_4D_spea(:,:,:,i_run) = map.img;
            end;                
        end;
        avg_pear = zeros(x,y,z);
        avg_spea = zeros(x,y,z);
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z
                    if pear == 1
                        avg_temp = mean(avg_corr_4D_pear(ind_x,ind_y,ind_z,:));
                        avg_pear(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;
                    if spea == 1
                        avg_temp = mean(avg_corr_4D_spea(ind_x,ind_y,ind_z,:));
                        avg_spea(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;                
                end;
            end;
        end;
        if pear == 1
            target_img = temp_img;
            target_img.fileprefix = 'CorrMaps_pear_split.nii';
            target_img.img = avg_pear;
            save_nii(target_img,target_img.fileprefix); 
        end;
        if spea == 1
            target_img = temp_img;
            target_img.fileprefix = 'CorrMaps_spea_split.nii';
            target_img.img = avg_spea;
            save_nii(target_img,target_img.fileprefix);        
        end;   
% comparison of splitted parametric regressor    
    if nr_para > 0
       for ind_para = 1:nr_para 
           for i = 1:runs
                if runs == 1
                    i = single_run;
                end;
                fprintf('...creates correlation maps for parametric contrast in splitted session %d...\n',i);

                %load 4D images
                eval(sprintf('one=img_%d_par%d_split1;',i,ind_para));
                one (~one) = nan;
                eval(sprintf('second=img_%d_par%d_split2;',i,ind_para));
                second (~second) = nan;

                % create correlation vectors 
                r_vec_pear_1_2 = zeros(x,y,z);
                r_vec_spea_1_2 = zeros(x,y,z);
                z_r_vec_pear_1_2 = zeros(x,y,z);
                z_r_vec_spea_1_2 = zeros(x,y,z);
    
                for ind_x = 1:x
                    for ind_y = 1:y
                        for ind_z = 1:z
                            first_voxel = one (ind_x, ind_y, ind_z, :);
                            second_voxel = second (ind_x, ind_y, ind_z, :);
                            % Pearson
                            if pear==1
                                r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                                if isnan (r(1,2))
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                else
                              % correction for understimation of reliability via
                              % split-half
                                    r = (2.*r(1,2))./(1+r(1,2));
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = r;
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r);
                                end;
                            end;
                            %Spearman
                            if spea==1
                                first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                                second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                                r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                                % correction for understimation of reliability via
                                % split-half
                                r = (2.*r)./(1+r);
                                r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                                z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);  
                            end;
                        end;
                    end;
                end;
        
                % save correlation map
                target_img = temp_img;
                file = sprintf('CorrMaps_pear_par%d_%d_split.nii',ind_para,i);
                target_img.fileprefix = file;
                target_img.img = r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);   

                target_img = temp_img;
                file = sprintf('CorrMaps_spea_par%d_%d_split.nii',ind_para,i);
                target_img.fileprefix = file;
                target_img.img = r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                target_img = temp_img;
                file = sprintf('z_CorrMaps_pear_par%d_%d_split.nii',ind_para,i);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);   

                target_img = temp_img;
                file = sprintf('z_CorrMaps_spea_par%d_%d_split.nii',ind_para,i);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                %compute z mean and inverse Fisher's transformation
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear_par%d_split_%d',ind_para,i);
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea_par%d_split_%d',ind_para,i);  
            end; 
    
    % calculate average correlation for all comparisons
        avg_corr_4D_pear = zeros(x,y,z,runs);
        avg_corr_4D_spea = zeros(x,y,z,runs);
        for i_run = 1:runs
            if pear == 1
                map = load_nii(sprintf('z_CorrMaps_pear_par%d_%d_split.nii',ind_para,i_run));
                avg_corr_4D_pear(:,:,:,i_run) = map.img;
            end;
            if spea == 1
                map = load_nii(sprintf('z_CorrMaps_spea_par%d_%d_split.nii',ind_para,i_run));
                avg_corr_4D_spea(:,:,:,i_run) = map.img;
            end;                
        end;
        avg_pear = zeros(x,y,z);
        avg_spea = zeros(x,y,z);
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z
                    if pear == 1
                        avg_temp = mean(avg_corr_4D_pear(ind_x,ind_y,ind_z,:));
                        avg_pear(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;
                    if spea == 1
                        avg_temp = mean(avg_corr_4D_spea(ind_x,ind_y,ind_z,:));
                        avg_spea(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;                
                end;
            end;
        end;
        if pear == 1
            target_img = temp_img;
            target_img.fileprefix = sprintf('CorrMaps_pear_par%d_split.nii',ind_para);
            target_img.img = avg_pear;
            save_nii(target_img,target_img.fileprefix); 
        end;
        if spea == 1
            target_img = temp_img;
            target_img.fileprefix = sprintf('CorrMaps_spea_par%d_split.nii',ind_para);
            target_img.img = avg_spea;
            save_nii(target_img,target_img.fileprefix);        
        end;   
        end;
    end;
%% two contrasts out of one statistic    
elseif two_cons == 1
    % compare contrasts within sessions
    for i_run = 1:runs
        fprintf('...\n compare contrasts of session %d\n...',i_run)
    
        if runs == 1
            i_run = single_run;
        end;
        eval(sprintf('one = img_%d_con1',i_run));
        one = one.img;
        one (~one) = nan;

        eval(sprintf('second = img_%d_con2',i_run));
        second = second.img;
        second (~second) = nan;
        
        % create correlation vectors 
        r_vec_pear_1_2 = zeros(x,y,z);
        r_vec_spea_1_2 = zeros(x,y,z);
        z_r_vec_pear_1_2 = zeros(x,y,z);
        z_r_vec_spea_1_2 = zeros(x,y,z);

        for ind_x = 1:x
            fprintf('..correlations voxels x = %d\n',ind_x)
            for ind_y = 1:y
                for ind_z = 1:z
                    first_voxel = one (ind_x, ind_y, ind_z, :);
                    second_voxel = second (ind_x, ind_y, ind_z, :);
                    % Pearson
                    if pear==1
                        r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                        if isnan (r(1,2))
                            r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                            z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                        else
                            r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
                            z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
                        end;
                    end;
                    %Spearman
                    if spea==1
                        first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                        second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                        r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                        r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                        z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);                    
                    end;
                end;
            end;
        end;
        
            % save correlation map
            fprintf('...saving maps...\n')
            target_img = temp_img;
            file = sprintf('CorrMaps_pear_%d_con1_con2.nii',i_run);
            target_img.fileprefix = file;
            target_img.img = r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_CorrMaps_pear_%d_con1_con2.nii',i_run);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);                     

            target_img = temp_img;
            file = sprintf('CorrMaps_spea_%d_con1_con2.nii',i_run);
            target_img.fileprefix = file;
            target_img.img = r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps_spea_%d_con1_con2.nii',i_run);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);  

            %compute z mean and inverse Fisher's transformation
            mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_pear_%d_con1_con2',i_run);
            mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_spea_%d_con1_con2',i_run);

    end;
    
    % calculate average correlation for all comparisons
        avg_corr_4D_pear = zeros(x,y,z,runs);
        avg_corr_4D_spea = zeros(x,y,z,runs);
        for i_run = 1:runs
            if pear == 1
                map = load_nii(sprintf('z_CorrMaps_pear_%d_con1_con2.nii',i_run));
                avg_corr_4D_pear(:,:,:,i_run) = map.img;
            end;
            if spea == 1
                map = load_nii(sprintf('z_CorrMaps_spea_%d_con1_con2.nii',i_run));
                avg_corr_4D_spea(:,:,:,i_run) = map.img;
            end;                
        end;
        avg_pear = zeros(x,y,z);
        avg_spea = zeros(x,y,z);
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z
                    if pear == 1
                        avg_temp = mean(avg_corr_4D_pear(ind_x,ind_y,ind_z,:));
                        avg_pear(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;
                    if spea == 1
                        avg_temp = mean(avg_corr_4D_spea(ind_x,ind_y,ind_z,:));
                        avg_spea(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;                
                end;
            end;
        end;
        if pear == 1
            target_img = temp_img;
            target_img.fileprefix = 'CorrMaps_pear_con1_con2.nii';
            target_img.img = avg_pear;
            save_nii(target_img,target_img.fileprefix); 
        end;
        if spea == 1
            target_img = temp_img;
            target_img.fileprefix = 'CorrMaps_spea_con1_con2.nii';
            target_img.img = avg_spea;
            save_nii(target_img,target_img.fileprefix);        
        end;   
    
    % compare contrasts between sessions
    for i_con = 1:2
        eval(sprintf('con=con%d;',i_con));
        eval(sprintf('con_count = con%d_count;',i_con));
        count_comp = 0;
    for i_run = 1:runs
        for count = 1:runs-1
            if i_run+count <= runs
            if runs > 1
                count_comp = count_comp+1;
                fprintf('...\n compare contrast %s between sessions\n...',con)
                eval(sprintf('one = img_%d_con1',i_run));
                one = one.img;
                one (~one) = nan;

                eval(sprintf('second = img_%d_con1',i_run+count));
                second = second.img;
                second (~second) = nan;
                % create correlation vectors 
                r_vec_pear_1_2 = zeros(x,y,z);
                r_vec_spea_1_2 = zeros(x,y,z);
                z_r_vec_pear_1_2 = zeros(x,y,z);
                z_r_vec_spea_1_2 = zeros(x,y,z);

                for ind_x = 1:x
                    fprintf('..correlations voxels x = %d\n',ind_x)
                    for ind_y = 1:y
                        for ind_z = 1:z
                            first_voxel = one (ind_x, ind_y, ind_z, :);
                            second_voxel = second (ind_x, ind_y, ind_z, :);
                            % Pearson
                            if pear==1
                                r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                                if isnan (r(1,2))
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                                else
                                    r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
                                    z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
                                end;
                            end;
                            %Spearman
                            if spea==1
                                first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                                second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                                r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                                r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                                z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);                    
                            end;
                        end;
                    end;
                end;

                % save correlation map
                fprintf('...saving maps...\n')
                target_img = temp_img;
                file = sprintf('CorrMaps_pear_%d_%d_con%d.nii',i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);  

                target_img = temp_img;
                file = sprintf('z_CorrMaps_pear_%d_%d_con%d.nii',i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);                     

                target_img = temp_img;
                file = sprintf('CorrMaps_spea_%d_%d_con%d.nii',i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                target_img = temp_img;
                file = sprintf('z_CorrMaps_spea_%d_%d_con%d.nii',i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);  

                %compute z mean and inverse Fisher's transformation
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear_%d_%d_con%d',i_run,i_run+count,con_count);
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea_%d_%d_con%d',i_run,i_run+count,con_count);
            end;
            end;
        end;
    end;
    end;
 
    % calculate average correlation for all comparisons
    for i_con = 1:2
        eval(sprintf('con=con%d;',i_con));
        eval(sprintf('con_count = con%d_count;',i_con));

        avg_corr_4D_pear = zeros(x,y,z,count_comp);
        avg_corr_4D_spea = zeros(x,y,z,count_comp);
        ind_comp = 0;
        for i_run = 1:runs
            for i_sec = 1:runs-i_run
                if i_run+i_sec <= runs
                    ind_comp = ind_comp+1;
                    if pear == 1
                        map = load_nii(sprintf('z_CorrMaps_pear_%d_%d_con%d.nii',i_run,i_run+i_sec,con_count));
                        avg_corr_4D_pear(:,:,:,ind_comp) = map.img;
                    end;
                    if spea == 1
                        map = load_nii(sprintf('z_CorrMaps_spea_%d_%d_con%d.nii',i_run,i_run+i_sec,con_count));
                        avg_corr_4D_spea(:,:,:,ind_comp) = map.img;
                    end;                
                end; 
            end;
        end;
        avg_pear = zeros(x,y,z);
        avg_spea = zeros(x,y,z);
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z
                    if pear == 1
                        avg_temp = mean(avg_corr_4D_pear(ind_x,ind_y,ind_z,:));
                        avg_pear(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;
                    if spea == 1
                        avg_temp = mean(avg_corr_4D_spea(ind_x,ind_y,ind_z,:));
                        avg_spea(ind_x,ind_y,ind_z) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
                    end;                
                end;
            end,
        end;
        if pear == 1
            target_img = temp_img;
            target_img.fileprefix = sprintf('CorrMaps_pear_con%d.nii',con_count);
            target_img.img = avg_pear;
            save_nii(target_img,target_img.fileprefix); 
        end;
        if spea == 1
            target_img = temp_img;
            target_img.fileprefix = sprintf('CorrMaps_spea_con%d.nii',con_count);
            target_img.img = avg_spea;
            save_nii(target_img,target_img.fileprefix);        
        end;    
    end;
    
    
    % parametric modulators in two contrasts not yet implemented
%     for i_run = 1:runs
%         fprintf('...\n compare parametric modulators for contrasts of session %d\n...',i_run)
%                 eval(sprintf('one = img_%d_con1',i_run));
%                 one = one.img;
%                 one (~one) = nan;
% 
%                 eval(sprintf('second = img_%d_con1',i_run+count));
%                 second = second.img;
%                 second (~second) = nan;
%         % create correlation vectors 
%         r_vec_pear_1_2 = zeros(x,y,z);
%         r_vec_spea_1_2 = zeros(x,y,z);
%         z_r_vec_pear_1_2 = zeros(x,y,z);
%         z_r_vec_spea_1_2 = zeros(x,y,z);
% 
%         for ind_x = 1:x
%             fprintf('..correlations voxels x = %d\n',ind_x)
%             for ind_y = 1:y
%                 for ind_z = 1:z
%                     first_voxel = one (ind_x, ind_y, ind_z, :);
%                     second_voxel = second (ind_x, ind_y, ind_z, :);
%                     % Pearson
%                     if pear==1
%                         r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
%                         if isnan (r(1,2))
%                             r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
%                             z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
%                         else
%                             r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
%                             z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
%                         end;
%                     end;
%                     %Spearman
%                     if spea==1
%                         first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
%                         second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
%                         r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
%                         r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
%                         z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);                    
%                     end;
%                 end;
%             end;
%         end;
%         
%         if runs == 1
%             % save correlation map
%             fprintf('...saving maps...\n')
%             target_img = temp_img;
%             file = sprintf('%s_pear_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = r_vec_pear_1_2;
%             save_nii(target_img,target_img.fileprefix);  
% 
%             target_img = temp_img;
%             file = sprintf('z_%s_pear_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = z_r_vec_pear_1_2;
%             save_nii(target_img,target_img.fileprefix);                     
% 
%             target_img = temp_img;
%             file = sprintf('%s_spea_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = r_vec_spea_1_2;
%             save_nii(target_img,target_img.fileprefix);
% 
%             target_img = temp_img;
%             file = sprintf('z_%s_spea_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = z_r_vec_spea_1_2;
%             save_nii(target_img,target_img.fileprefix);  
% 
%             %compute z mean and inverse Fisher's transformation
%             mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
%             mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%             summary(1,end+1)=mean_r;
%             cols{1,end+1} = sprintf('mean_pear_%d_con%d_con%d',single_run,con1_count,con2_count);
%             mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
%             mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%             summary(1,end+1)=mean_r;
%             cols{1,end+1} = sprintf('mean_spea_%d_con%d_con%d',single_run,con1_count,con2_count);
%           
%         else
%             % save correlation map
%             fprintf('...saving maps...\n')
%             target_img = temp_img;
%             file = sprintf('%s_pear_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = r_vec_pear_1_2;
%             save_nii(target_img,target_img.fileprefix);  
% 
%             target_img = temp_img;
%             file = sprintf('z_%s_pear_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = z_r_vec_pear_1_2;
%             save_nii(target_img,target_img.fileprefix);                     
% 
%             target_img = temp_img;
%             file = sprintf('%s_spea_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = r_vec_spea_1_2;
%             save_nii(target_img,target_img.fileprefix);
% 
%             target_img = temp_img;
%             file = sprintf('z_%s_spea_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
%             target_img.fileprefix = file;
%             target_img.img = z_r_vec_spea_1_2;
%             save_nii(target_img,target_img.fileprefix);  
% 
%             %compute z mean and inverse Fisher's transformation
%             mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
%             mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%             summary(1,end+1)=mean_r;
%             cols{1,end+1} = sprintf('mean_pear_%d_con%d_con%d',i_run,con1_count,con2_count);
%             mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
%             mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%             summary(1,end+1)=mean_r;
%             cols{1,end+1} = sprintf('mean_spea_%d_con%d_con%d',i_run,con1_count,con2_count);
% 
%         end;
%     end;
%     % compare contrasts between sessions
%     for i_con = 1:2
%         eval(sprintf('cons=cons%d;',i_con));
%         eval(sprintf('con_count = cons%d_count;',i_con));
%     for i_run = 1:runs
%         for count = 1:runs-1
%             if i_run+count <= runs
%             if runs > 1
%                 fprintf('...\n compare contrast %s between sessions\n...',cons)
%                 eval(sprintf('one = img_%d_con1',i_run));
%                 one = one.img;
%                 one (~one) = nan;
% 
%                 eval(sprintf('second = img_%d_con1',i_run+count));
%                 second = second.img;
%                 second (~second) = nan;
%                 % create correlation vectors 
%                 r_vec_pear_1_2 = zeros(x,y,z);
%                 r_vec_spea_1_2 = zeros(x,y,z);
%                 z_r_vec_pear_1_2 = zeros(x,y,z);
%                 z_r_vec_spea_1_2 = zeros(x,y,z);
% 
%                 for ind_x = 1:x
%                     fprintf('..correlations voxels x = %d\n',ind_x)
%                     for ind_y = 1:y
%                         for ind_z = 1:z
%                             first_voxel = one (ind_x, ind_y, ind_z, :);
%                             second_voxel = second (ind_x, ind_y, ind_z, :);
%                             % Pearson
%                             if pear==1
%                                 r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
%                                 if isnan (r(1,2))
%                                     r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
%                                     z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
%                                 else
%                                     r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
%                                     z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
%                                 end;
%                             end;
%                             %Spearman
%                             if spea==1
%                                 first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
%                                 second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
%                                 r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
%                                 r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
%                                 z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);                    
%                             end;
%                         end;
%                     end;
%                 end;
% 
%                 % save correlation map
%                 fprintf('...saving maps...\n')
%                 target_img = temp_img;
%                 file = sprintf('%s_pear_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
%                 target_img.fileprefix = file;
%                 target_img.img = r_vec_pear_1_2;
%                 save_nii(target_img,target_img.fileprefix);  
% 
%                 target_img = temp_img;
%                 file = sprintf('z_%s_pear_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
%                 target_img.fileprefix = file;
%                 target_img.img = z_r_vec_pear_1_2;
%                 save_nii(target_img,target_img.fileprefix);                     
% 
%                 target_img = temp_img;
%                 file = sprintf('%s_spea_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
%                 target_img.fileprefix = file;
%                 target_img.img = r_vec_spea_1_2;
%                 save_nii(target_img,target_img.fileprefix);
% 
%                 target_img = temp_img;
%                 file = sprintf('z_%s_spea_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
%                 target_img.fileprefix = file;
%                 target_img.img = z_r_vec_spea_1_2;
%                 save_nii(target_img,target_img.fileprefix);  
% 
%                 %compute z mean and inverse Fisher's transformation
%                 mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
%                 mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%                 summary(1,end+1)=mean_r;
%                 cols{1,end+1} = sprintf('mean_pear_%d_%d_con%d',i_run,i_run+count,con_count);
%                 mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
%                 mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%                 summary(1,end+1)=mean_r;
%                 cols{1,end+1} = sprintf('mean_spea_%d_%d_con%d',i_run,i_run+count,con_count);
%             end;
%             end;
%         end;
%     end;
%     end;
end;


results_corr=dataset({summary(1,:),cols{:}});
assignin('base','results_corr',results_corr);
    
cd(dir_results);             
save results_corr.mat results_corr            
disp('...finished correlation maps')
end;
%% ICCs
if cons == 1 || abs == 1
data=zeros(nr_subj,runs);
ICC_con_ROI=zeros(x,y,z);
ICC_abs_ROI=zeros(x,y,z);
z_ICC_con_ROI=zeros(x,y,z);
z_ICC_abs_ROI=zeros(x,y,z);

% initiating summary
summary = [];
cols = {};

disp('...calculating ICCs...')
if split == 0 && two_cons == 0
    for ind_x = 1:x
        fprintf('..scanning voxel x = %d ...\n',ind_x)
        for ind_y = 1:y
           for ind_z = 1:z
                for ind_run = 1:runs
                    if runs == 1
                        ind_run = single_run;
                    end;
                    for ind_subj = 1:nr_subj
                       eval(sprintf('img%d_voxel(ind_subj,1)= img_%d (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run));
                    end;
                    eval(sprintf('data(:,ind_run)=img%d_voxel;',ind_run));
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
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                    z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;
                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                    z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;
           end; 
        end;
     end;
disp('saving ICC images')
 % save ICC maps
target_img = temp_img;
target_img.fileprefix = 'ICC_con.nii';
target_img.img = ICC_con_ROI;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = 'ICC_abs.nii';
target_img.img = ICC_abs_ROI;
save_nii(target_img,target_img.fileprefix); 

target_img = temp_img;
target_img.fileprefix = 'z_ICC_con.nii';
target_img.img = z_ICC_con_ROI;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = 'z_ICC_abs.nii';
target_img.img = z_ICC_abs_ROI;
save_nii(target_img,target_img.fileprefix); 

% computing average ICC
disp('...computing ROI total ICC...')
ROI_subj = zeros(nr_subj,runs);
for ind_subj = 1:nr_subj
    for ind_run = 1:runs
        if runs == 1
            ind_run = single_run;
        end;
        eval(sprintf('img%d_%d = img_%d (:,:,:,ind_subj);',ind_run,ind_subj,ind_run));
        eval(sprintf('ROI_subj(ind_subj,ind_run) = mean(img%d_%d(:));',ind_run,ind_subj));
    end;
end;                
                    %ICC
                    nsamples=nr_subj*runs;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:runs,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(runs,1);
                    for sess=1:runs
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:runs
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
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
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= 'ICC_con';
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}='z_ICC_con';
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= 'ICC_abs';
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}='z_ICC_abs';
                    summary(1,end+1)=z_ICC_abs_ROI_total;                    
                    end;

% computing ROI ICC
if roi == 1
disp('...computing ROI total ICC...')
ROI_subj = zeros(nr_subj,runs);
for ind_subj = 1:nr_subj
    for ind_run = 1:runs
        if runs == 1
            ind_run = single_run;
        end;
        eval(sprintf('img%d_%d = img_%d (:,:,:,ind_subj);',ind_run,ind_subj,ind_run));
        eval(sprintf('ROI_subj(ind_subj,ind_run) = mean(img%d_%d(r_roi_ind));',ind_run,ind_subj));
    end;
end;                
                    %ICC
                    nsamples=nr_subj*runs;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:runs,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(runs,1);
                    for sess=1:runs
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:runs
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
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
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= 'ICC_con';
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}='z_ICC_con';
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= 'ICC_abs';
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}='z_ICC_abs';
                    summary(1,end+1)=z_ICC_abs_ROI_total;                    
                    end;
end;
                    assignin('base','summary',summary);
                    assignin('base','cols',cols);

                    results_ICC=dataset({summary(1,:),cols{:}});
                    assignin('base','results_ICC',results_ICC);

                    cd(dir_results);             
                    save results_ICC.mat results_ICC;

if nr_para > 0
    for ind_par = 1:nr_para
        for ind_x = 1:x
            fprintf('..scanning voxel x = %d ...\n',ind_x)
            for ind_y = 1:y
               for ind_z = 1:z
                    for ind_run = 1:runs
                        if runs == 1
                            ind_run = single_run;
                        end;
                        for ind_subj = 1:nr_subj
                           estr = sprintf('img%d_par%d_voxel(ind_subj,1)= img_par%d_%d (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_par,ind_par,ind_run);
                           eval(estr);
                        end;
                        estr = sprintf('data(:,ind_run)=img%d_par%d_voxel;',ind_run,ind_par);
                        eval(estr);
                    end;
                        %ICC over all sessions
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
                        if cons==1
                        voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                        ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                        end;

                        %absolute agreement 
                        if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                           runs.* (JMS-EMS)./nr_subj);
                        ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                        end;
               end; 
            end;
        end;
        disp('saving ICC images')
         % save ICC maps
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_con_par%d.nii',ind_par);
        target_img.img = ICC_con_ROI;
        save_nii(target_img,target_img.fileprefix); 

        clear target_img;
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_abs_par%d.nii',ind_par);
        target_img.img = ICC_abs_ROI;
        save_nii(target_img,target_img.fileprefix); 

        target_img = temp_img;
        target_img.fileprefix = sprintf('z_ICC_con_par%d.nii',ind_par);
        target_img.img = z_ICC_con_ROI;
        save_nii(target_img,target_img.fileprefix); 

        clear target_img;
        target_img = temp_img;
        target_img.fileprefix = sprintf('z_ICC_abs_par%d.nii',ind_par);
        target_img.img = z_ICC_abs_ROI;
        save_nii(target_img,target_img.fileprefix); 

        % computing ROI ICC
        if roi == 1
        disp('...computing ROI total ICC...')
        ROI_subj = zeros(nr_subj,runs);
        for ind_subj = 1:nr_subj
            for ind_run = 1:runs
                if runs == 1
                    ind_run = single_run;
                end;
                eval(sprintf('img%d_%d = img_par%d_%d (:,:,:,ind_subj);',ind_run,ind_subj,ind_par,ind_run));
                eval(sprintf('ROI_subj(ind_subj,ind_run) = mean(img%d_%d(r_roi_ind));',ind_run,ind_subj));
            end;
        end;                
                    %ICC
                    nsamples=nr_subj*runs;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:runs,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(runs,1);
                    for sess=1:runs
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:runs
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
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
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_ROI_par%d',ind_par);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_ROI_par%d',ind_par);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_ROI_par%d',ind_par);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_ROI_par%d',ind_par);
                    summary(1,end+1)=z_ICC_abs_ROI_total;                    
                    end;   

        end;
        % computing average ICC
        disp('...computing total ICC...')
        ROI_subj = zeros(nr_subj,runs);
        for ind_subj = 1:nr_subj
            for ind_run = 1:runs
                if runs == 1
                    ind_run = single_run;
                end;
                eval(sprintf('img%d_%d = img_par%d_%d (:,:,:,ind_subj);',ind_run,ind_subj,ind_par,ind_run));
                eval(sprintf('ROI_subj(ind_subj,ind_run) = mean(img%d_%d(:));',ind_run,ind_subj));
            end;
        end;                
                    %ICC
                    nsamples=nr_subj*runs;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:runs,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(runs,1);
                    for sess=1:runs
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:runs
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
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
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_par%d',ind_par);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_par%d',ind_par);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_par%d',ind_par);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_par%d',ind_par);
                    summary(1,end+1)=z_ICC_abs_ROI_total;                    
                    end;   
end;
                    assignin('base','summary',summary);
                    assignin('base','cols',cols);

                    results_ICC=dataset({summary(1,:),cols{:}});
                    assignin('base','results_ICC',results_ICC);

                    cd(dir_results);             
                    save results_ICC_par.mat results_ICC;
end;
%calculation ICC for each comparison
if runs > 1
    for ind_run = 1:runs
        for ind_sec = 1:runs-1
            if ind_run + ind_sec <= runs
            for ind_x = 1:x
                fprintf('...x = %d...',ind_x)
                for ind_y = 1:y
                    for ind_z = 1:z
                        data = zeros(nr_subj,2);
                    for ind_subj = 1:nr_subj
                       eval(sprintf('img%d_voxel(ind_subj,1)= img_%d(ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run));
                       eval(sprintf('img%d_voxel(ind_subj,1)= img_%d(ind_x, ind_y, ind_z, ind_subj);',ind_run+ind_sec,ind_run+ind_sec));                       
                    end;
                    eval(sprintf('data(:,1)=img%d_voxel;',ind_run,ind_run));
                    eval(sprintf('data(:,2)=img%d_voxel;',ind_run+ind_sec,ind_run+ind_sec));                    
                    %ICC
    
                        nsamples=nr_subj*runs;
                        
                        grandmean=0;
                        for sub=1:nr_subj,     
                            for sess=1:2,
                               grandmean= grandmean + data(sub,sess);
                            end
                        end;
                        grandmean=grandmean./nsamples;

                        sessionmean=zeros(2,1);
                        for sess=1:2
                            for sub=1:nr_subj,  
                                sessionmean(sess) = sessionmean(sess) + data(sub,sess);
                            end
                            sessionmean(sess)=sessionmean(sess)./nr_subj;
                        end

                        subjmean=zeros(nr_subj,1);
                        for sub=1:nr_subj
                            for sess=1:2
                                subjmean(sub)=subjmean(sub) + data(sub,sess);
                            end
                              subjmean(sub)=subjmean(sub)./2;
                        end

                        % mean squares
                        BMS=0; % between subject
                        WMS=0; % within subject 
                        EMS=0; % error
                        JMS=0; % session

                        for sub=1:nr_subj,    
                            BMS = BMS + (subjmean(sub)-grandmean).^2;
                            for sess=1:2
                                WMS = WMS + (data(sub,sess)-subjmean(sub)).^2;
                                EMS = EMS + (data(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                            end
                        end;

                        for sess=1:2
                            JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                        end;

                        %define the true value of the mean square.
                        BMS= 2.*BMS./(nr_subj-1);
                        WMS= WMS./(2-1)./nr_subj;
                        JMS= nr_subj.*JMS./(2-1);
                        EMS= EMS./(2-1)./(nr_subj-1); 
                        
                        fprintf('...calculate voxel ICC x = %d, y = %d, z = %d...\n',ind_x,ind_y,ind_z)
                        %consistency agreement  
                        if cons==1
                        voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                        ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                        end;

                        %absolute agreement 
                        if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                           2.* (JMS-EMS)./nr_subj);
                        ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                        end;

                    end;
                 end; 
              end;
                             disp('saving ICC images')
                             % save ICC maps
                            target_img = temp_img;
                            file1 = sprintf('ICC_con_%d_%d.nii',ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file1;
                            target_img.img = ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file2 = sprintf('ICC_abs_%d_%d.nii',ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file2;
                            target_img.img = ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            target_img = temp_img;
                            file3 = sprintf('z_ICC_con_%d_%d.nii',ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file3;
                            target_img.img = z_ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file4 = sprintf('z_ICC_abs_%d_%d.nii',ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file4;
                            target_img.img = z_ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix);   
                            
     % computing ROI ICC for each comparison
     if roi == 1 
        disp('...computing ROI total ICC...')
            ROI_subj = zeros(nr_subj,2);
            for ind_subj = 1:nr_subj
                eval(sprintf('temp = img_%d(:,:,:,ind_subj);',ind_run));
                ROI_subj(ind_subj,1) = mean(temp(r_roi_ind));
                eval(sprintf('temp1 = img_%d(:,:,:,ind_subj);',ind_run+ind_sec));
                ROI_subj(ind_subj,2) = mean(temp1(r_roi_ind));
            end;                
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./2;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session
                    
                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:2
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_%d_%d',ind_run,ind_run+ind_sec);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_%d_%d',ind_run,ind_run+ind_sec);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    end;
                    

                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_%d_%d',ind_run,ind_run+ind_sec);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_%d_%d',ind_run,ind_run+ind_sec);
                    summary(1,end+1)=z_ICC_abs_ROI_total;
                    end;
    end;
            end;
                    assignin('base','summary',summary);
                    assignin('base','cols',cols);

                    results_ICC=dataset({summary(1,:),cols{:}});
                    assignin('base','results_ICC',results_ICC);

                    cd(dir_results);             
                    save results_ICC.mat results_ICC;
     end;                            
  end;
end;

if nr_para > 0
    fprintf('...ICC parametric modulator...\n')
    for ind_para = 1:nr_para
%calculation ICC for each comparison
if runs > 1
    for ind_run = 1:runs
        for ind_sec = 1:runs-1
            if ind_run + ind_sec <= runs

            for ind_x = 1:x
                fprintf('...x = %d...',ind_x)
                for ind_y = 1:y
                    for ind_z = 1:z
                        data = zeros(nr_subj,2);
                    for ind_subj = 1:nr_subj
                       eval(sprintf('img%d_voxel(ind_subj,1)= img_par%d_%d(ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_para,ind_run));
                       eval(sprintf('img%d_voxel(ind_subj,1)= img_par%d_%d(ind_x, ind_y, ind_z, ind_subj);',ind_run+ind_sec,ind_para,ind_run+ind_sec));                       
                    end;
                    eval(sprintf('data(:,1)=img%d_voxel;',ind_run,ind_run));
                    eval(sprintf('data(:,2)=img%d_voxel;',ind_run+ind_sec,ind_run+ind_sec));                    
                    %ICC
    
                        nsamples=nr_subj*runs;
                        
                        grandmean=0;
                        for sub=1:nr_subj,     
                            for sess=1:2,
                               grandmean= grandmean + data(sub,sess);
                            end
                        end;
                        grandmean=grandmean./nsamples;

                        sessionmean=zeros(2,1);
                        for sess=1:2
                            for sub=1:nr_subj,  
                                sessionmean(sess) = sessionmean(sess) + data(sub,sess);
                            end
                            sessionmean(sess)=sessionmean(sess)./nr_subj;
                        end

                        subjmean=zeros(nr_subj,1);
                        for sub=1:nr_subj
                            for sess=1:2
                                subjmean(sub)=subjmean(sub) + data(sub,sess);
                            end
                              subjmean(sub)=subjmean(sub)./2;
                        end

                        % mean squares
                        BMS=0; % between subject
                        WMS=0; % within subject 
                        EMS=0; % error
                        JMS=0; % session

                        for sub=1:nr_subj,    
                            BMS = BMS + (subjmean(sub)-grandmean).^2;
                            for sess=1:2
                                WMS = WMS + (data(sub,sess)-subjmean(sub)).^2;
                                EMS = EMS + (data(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                            end
                        end;

                        for sess=1:2
                            JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                        end;

                        %define the true value of the mean square.
                        BMS= 2.*BMS./(nr_subj-1);
                        WMS= WMS./(2-1)./nr_subj;
                        JMS= nr_subj.*JMS./(2-1);
                        EMS= EMS./(2-1)./(nr_subj-1); 
                        
                        fprintf('...calculate voxel ICC x = %d, y = %d, z = %d...\n',ind_x,ind_y,ind_z)
                        %consistency agreement  
                        if cons==1
                        voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                        ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                        end;

                        %absolute agreement 
                        if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                           2.* (JMS-EMS)./nr_subj);
                        ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                        end;

                    end;
                 end; 
              end;
                             disp('saving ICC images')
                             % save ICC maps
                            target_img = temp_img;
                            file1 = sprintf('ICC_con_par%d_%d_%d.nii',ind_para,ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file1;
                            target_img.img = ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file2 = sprintf('ICC_abs_par%d_%d_%d.nii',ind_para,ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file2;
                            target_img.img = ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            target_img = temp_img;
                            file3 = sprintf('z_ICC_con_par%d_%d_%d.nii',ind_para,ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file3;
                            target_img.img = z_ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file4 = sprintf('z_ICC_abs_par%d_%d_%d.nii',ind_para,ind_run,ind_run+ind_sec);
                            target_img.fileprefix = file4;
                            target_img.img = z_ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix);   
            end;
                            
     % computing ROI ICC for each comparison
     if roi == 1 
        disp('...computing ROI total ICC...')
            ROI_subj = zeros(nr_subj,2);
            for ind_subj = 1:nr_subj
                eval(sprintf('temp = img_par%d_%d(:,:,:,ind_subj);',ind_para,ind_run));
                ROI_subj(ind_subj,1) = mean(temp(r_roi_ind));
                eval(sprintf('temp1 = img_par%d_%d(:,:,:,ind_subj);',ind_para,ind_run+ind_sec));
                ROI_subj(ind_subj,2) = mean(temp1(r_roi_ind));
            end;                
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./2;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session
                    
                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:2
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_par%d_%d_%d',ind_para,ind_run,ind_run+ind_sec);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_par%d_%d_%d',ind_para,ind_run,ind_run+ind_sec);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    end;
                    

                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_par%d_%d_%d',ind_para,ind_run,ind_run+ind_sec);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_par%d_%d_%d',ind_para,ind_run,ind_run+ind_sec);
                    summary(1,end+1)=z_ICC_abs_ROI_total;
                    end;
    end;
        end;

     end;                            
  end;
end;
end;
                    assignin('base','summary',summary);
                    assignin('base','cols',cols);

                    results_ICC=dataset({summary(1,:),cols{:}});
                    assignin('base','results_ICC',results_ICC);

                    cd(dir_results);             
                    save results_ICC.mat results_ICC;
%% based on split-half
elseif split == 1
for ind = 1:runs
     for ind_x = 1:x
          fprintf('...x = %d...',ind_x)
        for ind_y = 1:y
           for ind_z = 1:z
                    if runs == 1
                        ind_run = single_run;
                    else
                        ind_run = ind;
                    end;
                    eval(sprintf('data_%d = zeros(nr_subj,2);',ind_run));
                    for ind_subj = 1:nr_subj
                       eval(sprintf('img%d_voxel_split1(ind_subj,1)= img_%d_split1 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run));
                       eval(sprintf('img%d_voxel_split2(ind_subj,1)= img_%d_split2 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run));                       
                    end;
                    eval(sprintf('data_%d(:,1)=img%d_voxel_split1;',ind_run,ind_run));
                    eval(sprintf('data_%d(:,2)=img%d_voxel_split2;',ind_run,ind_run));                    
                    %ICC
    
                        nsamples=nr_subj*runs;
                        
                        grandmean=0;
                        for sub=1:nr_subj,     
                            for sess=1:2,
                               eva=sprintf('grandmean= grandmean + data_%d(sub,sess);',ind);
                               eval(eva);
                            end
                        end;
                        grandmean=grandmean./nsamples;

                        sessionmean=zeros(2,1);
                        for sess=1:2
                            for sub=1:nr_subj,  
                                e = sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess);',ind);
                                eval(e);
                            end
                            sessionmean(sess)=sessionmean(sess)./nr_subj;
                        end

                        subjmean=zeros(nr_subj,1);
                        for sub=1:nr_subj
                            for sess=1:2
                                ev = sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess);',ind);
                                eval(ev);
                            end
                              subjmean(sub)=subjmean(sub)./2;
                        end

                        % mean squares
                        BMS=0; % between subject
                        WMS=0; % within subject 
                        EMS=0; % error
                        JMS=0; % session

                        for sub=1:nr_subj,    
                            BMS = BMS + (subjmean(sub)-grandmean).^2;
                            for sess=1:2
                                evs = sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2;',ind);
                                eval(evs);
                                evst = sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;',ind);
                                eval(evst);
                            end
                        end;

                        for sess=1:2
                            JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                        end;

                        %define the true value of the mean square.
                        BMS= 2.*BMS./(nr_subj-1);
                        WMS= WMS./(2-1)./nr_subj;
                        JMS= nr_subj.*JMS./(2-1);
                        EMS= EMS./(2-1)./(nr_subj-1); 
                        
                        fprintf('...calculate voxel ICC x = %d, y = %d, z = %d...\n',ind_x,ind_y,ind_z)
                        %consistency agreement  
                        if cons==1
                        voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                        ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                        end;

                        %absolute agreement 
                        if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                           2.* (JMS-EMS)./nr_subj);
                        ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                        end;
                                                
            end;
         end; 
      end;
                            disp('saving ICC images')
                             % save ICC maps
                            target_img = temp_img;
                            file1 = sprintf('ICC_con_%d_split.nii',ind);
                            target_img.fileprefix = file1;
                            target_img.img = ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file2 = sprintf('ICC_abs_%d_split.nii',ind);
                            target_img.fileprefix = file2;
                            target_img.img = ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            target_img = temp_img;
                            file3 = sprintf('z_ICC_con_%d_split.nii',ind);
                            target_img.fileprefix = file3;
                            target_img.img = z_ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file4 = sprintf('z_ICC_abs_%d_split.nii',ind);
                            target_img.fileprefix = file4;
                            target_img.img = z_ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix);              
end;
   
     

     % computing ROI ICC
if roi == 1
disp('...computing ROI total ICC...')

for ind_runs = 1:runs
    if runs == 1
        ind_runs = single_run;
    end;
    ROI_subj = zeros(nr_subj,2);
    for ind_subj = 1:nr_subj
        for ind_split = 1:2
            eval(sprintf('img%d_%d = img_%d_split%d (:,:,:,ind_subj);',ind_split,ind_subj,ind_run,ind_split));
            eval(sprintf('ROI_subj(ind_subj,ind_split) = mean(img%d_%d(r_roi_ind));',ind_split,ind_subj));
        end;
    end;                
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./2;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session
                    
                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:2
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_ROI_split%d',ind_runs);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_ROI_split%d',ind_runs);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    end;
                    

                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_ROI_split%d',ind_runs);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_ROI_split%d',ind_runs);
                    summary(1,end+1)=z_ICC_abs_ROI_total;
                    end;


end;
end;
                    assignin('base','summary',summary);
                    assignin('base','cols',cols);

                    results_ICC=dataset({summary(1,:),cols{:}});
                    assignin('base','results_ICC',results_ICC);

                    cd(dir_results);             
                    save results_ICC_split.mat results_ICC;
if nr_para > 0
    fprintf('...ICC parametric modulator for splitted sessions...\n')
        for ind_para = 1:nr_para
             for ind_run = 1:runs
             if runs == 1
                 ind_run = single_run;
             end;
             for ind_x = 1:x
                for ind_y = 1:y
                   for ind_z = 1:z
                 fprintf('... x = %d, y = %d, z = %d ...\n',ind_x,ind_y,ind_z)

                            eval(sprintf('data_%d = zeros(nr_subj,2);',ind_run));
                            for ind_subj = 1:nr_subj
                               eval(sprintf('img%d_par%d_voxel_split1(ind_subj,1)= img_%d_par%d_split1 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_para,ind_run,ind_para));
                               eval(sprintf('img%d_par%d_voxel_split2(ind_subj,1)= img_%d_par%d_split2 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_para,ind_run,ind_para));
                            end;
                            eval(sprintf('data_%d(:,1)=img%d_par%d_voxel_split1;',ind_run,ind_run,ind_para));
                            eval(sprintf('data_%d(:,2)=img%d_par%d_voxel_split2;',ind_run,ind_run,ind_para));

                  %ICC
                                nsamples=nr_subj*2;

                                grandmean=0;
                                for sub=1:nr_subj,     
                                    for sess=1:2,
                                       eva=sprintf('grandmean= grandmean + data_%d(sub,sess);',ind_run);
                                       eval(eva);
                                    end
                                end;
                                grandmean=grandmean./nsamples;

                                sessionmean=zeros(2,1);
                                for sess=1:2
                                    for sub=1:nr_subj,  
                                        e = sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess);',ind_run);
                                        eval(e);
                                    end
                                    sessionmean(sess)=sessionmean(sess)./nr_subj;
                                end

                                subjmean=zeros(nr_subj,1);
                                for sub=1:nr_subj
                                    for sess=1:2
                                        ev = sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess);',ind_run);
                                        eval(ev);
                                    end
                                      subjmean(sub)=subjmean(sub)./2;
                                end

                                % mean squares
                                BMS=0; % between subject
                                WMS=0; % within subject 
                                EMS=0; % error
                                JMS=0; % session

                                for sub=1:nr_subj,    
                                    BMS = BMS + (subjmean(sub)-grandmean).^2;
                                    for sess=1:2
                                        evs = sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2;',ind_run);
                                        eval(evs);
                                        evst = sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;',ind_run);
                                        eval(evst);
                                    end
                                end;

                                for sess=1:2
                                    JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                                end;

                                %define the true value of the mean square.
                                BMS= 2.*BMS./(nr_subj-1);
                                WMS= WMS./(2-1)./nr_subj;
                                JMS= nr_subj.*JMS./(2-1);
                                EMS= EMS./(2-1)./(nr_subj-1); 

                                %consistency agreement  
                                if cons==1
                                voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                                ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                                z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                                end;

                                %absolute agreement 
                                if abs==1
                                voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                                   2.* (JMS-EMS)./nr_subj);
                                ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                                z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                                end;
                    end;
                end; 
             end;

                                             disp('saving ICC images')
                                     % save ICC maps
                                    target_img = temp_img;
                                    file1 = sprintf('ICC_con_%d_split_par%d.nii',ind_run,ind_para);
                                    target_img.fileprefix = file1;
                                    target_img.img = ICC_con_ROI;
                                    save_nii(target_img,target_img.fileprefix); 

                                    clear target_img;
                                    target_img = temp_img;
                                    file2 = sprintf('ICC_abs_%d_split_par%d.nii',ind_run,ind_para);
                                    target_img.fileprefix = file2;
                                    target_img.img = ICC_abs_ROI;
                                    save_nii(target_img,target_img.fileprefix); 

                                    target_img = temp_img;
                                    file3 = sprintf('z_ICC_con_%d_split_par%d.nii',ind_run,ind_para);
                                    target_img.fileprefix = file3;
                                    target_img.img = z_ICC_con_ROI;
                                    save_nii(target_img,target_img.fileprefix); 

                                    clear target_img;
                                    target_img = temp_img;
                                    file4 = sprintf('z_ICC_abs_%d_split_par%d.nii',ind_run,ind_para);
                                    target_img.fileprefix = file4;
                                    target_img.img = z_ICC_abs_ROI;
                                    save_nii(target_img,target_img.fileprefix); 
        end;
        end;   
        
if roi == 1        
for ind_para = 1:nr_para

for ind_runs = 1:runs
    if runs == 1
        ind_runs = single_run;
    end;
    ROI_subj = zeros(nr_subj,2);
    for ind_subj = 1:nr_subj
        for ind_split = 1:2
            eval(sprintf('img%d_%d = img_%d_par%d_split%d (:,:,:,ind_subj);',ind_split,ind_subj,ind_runs,ind_para,ind_split));
            eval(sprintf('ROI_subj(ind_subj,ind_split) = mean(img%d_%d(r_roi_ind));',ind_split,ind_subj));
        end;
    end;                
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./2;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session
                    
                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:2
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    
                    cols{1,end+1}= sprintf('ICC_con_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=z_ICC_con_ROI_total;                    
                    end;
                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    
                    cols{1,end+1}= sprintf('ICC_abs_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=z_ICC_abs_ROI_total;                     
                    end;
                    
                  
                    
end;                  
                    assignin('base','summary',summary);
                    assignin('base','cols',cols);

                    results_ICC=dataset({summary(1,:),cols{:}});
                    assignin('base','results_ICC',results_ICC);

                    cd(dir_results);             
                    sprintf('save results_ICC_split_par%d.mat results_ICC',ind_para);
end;
end;
end;

%% two contrasts out of one statistic
elseif two_cons == 1
for ind = 1:runs
    
    if runs == 1
        ind = single_run;
    end;
    disp('...comparison of two contrasts in session %d...', ind);
 
     eval(sprintf('data_%d = zeros(nr_subj,2);',ind));
    for ind_x = 1:x
          fprintf('...x = %d...\n',ind_x)
        for ind_y = 1:y
            for ind_z = 1:z
                    for ind_subj = 1:nr_subj
                       eval(sprintf('img%d_voxel_con1(ind_subj,1)= img_%d_con1 (ind_x, ind_y, ind_z, ind_subj);',ind,ind));
                       eval(sprintf('img%d_voxel_con2(ind_subj,1)= img_%d_con2 (ind_x, ind_y, ind_z, ind_subj);',ind,ind));
                    end;
                    eval(sprintf('data_%d(:,1)=img%d_voxel_con1;',ind,ind));
                    eval(sprintf('data_%d(:,2)=img%d_voxel_con2;',ind,ind));
                    %ICC
    
                        nsamples=nr_subj*runs;
                        
                        grandmean=0;
                        for sub=1:nr_subj,     
                            for sess=1:2,
                               eva=sprintf('grandmean= grandmean + data_%d(sub,sess);',ind);
                               eval(eva);
                            end
                        end;
                        grandmean=grandmean./nsamples;

                        sessionmean=zeros(2,1);
                        for sess=1:2
                            for sub=1:nr_subj,  
                                e = sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess);',ind);
                                eval(e);
                            end
                            sessionmean(sess)=sessionmean(sess)./nr_subj;
                        end

                        subjmean=zeros(nr_subj,1);
                        for sub=1:nr_subj
                            for sess=1:2
                                ev = sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess);',ind);
                                eval(ev);
                            end
                              subjmean(sub)=subjmean(sub)./2;
                        end

                        % mean squares
                        BMS=0; % between subject
                        WMS=0; % within subject 
                        EMS=0; % error
                        JMS=0; % session

                        for sub=1:nr_subj,    
                            BMS = BMS + (subjmean(sub)-grandmean).^2;
                            for sess=1:2
                                evs = sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2;',ind);
                                eval(evs);
                                evst = sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;',ind);
                                eval(evst);
                            end
                        end;

                        for sess=1:2
                            JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                        end;

                        %define the true value of the mean square.
                        BMS= 2.*BMS./(nr_subj-1);
                        WMS= WMS./(2-1)./nr_subj;
                        JMS= nr_subj.*JMS./(2-1);
                        EMS= EMS./(2-1)./(nr_subj-1); 
                        
                        %consistency agreement  
                        if cons==1
                        voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                        ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                        end;

                        %absolute agreement 
                        if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                           2.* (JMS-EMS)./nr_subj);
                        ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                        end;
                                                
            end;
        end; 
    end;
                            disp('saving ICC images')
                             % save ICC maps
                            target_img = temp_img;
                            file1 = sprintf('ICC_con_%d_2cons.nii',ind);
                            target_img.fileprefix = file1;
                            target_img.img = ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file2 = sprintf('ICC_abs_%d_2cons.nii',ind);
                            target_img.fileprefix = file2;
                            target_img.img = ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            target_img = temp_img;
                            file3 = sprintf('z_ICC_con_%d_2cons.nii',ind);
                            target_img.fileprefix = file3;
                            target_img.img = z_ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file4 = sprintf('z_ICC_abs_%d_2cons.nii',ind);
                            target_img.fileprefix = file4;
                            target_img.img = z_ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix);              
end;
   
     

     % computing ROI ICC
disp('...computing ROI total ICC...')

if roi == 1
for ind_runs = 1:runs
    if runs == 1
        ind_runs = single_run;
    end;
    ROI_subj = zeros(nr_subj,2);
  for ind_con = 1:2
    for ind_subj = 1:nr_subj
            eval(sprintf('img%d_%d = img_%d_con%d (:,:,:,ind_subj);',ind_con,ind_subj,ind_runs,ind_con));
            eval(sprintf('ROI_subj(ind_subj,ind_con) = mean(img%d_%d(r_roi_ind));',ind_con,ind_subj));
    end;
  end;
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./2;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session
                    
                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:2
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_ROI_%d',ind_runs);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_ROI_%d',ind_runs);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    end;
                    

                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_ROI_%d',ind_runs);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_ROI_%d',ind_runs);
                    summary(1,end+1)=z_ICC_abs_ROI_total;
                    end;
end;
end;
                    assignin('base','summary',summary);
                    assignin('base','cols',cols);

                    results_ICC=dataset({summary(1,:),cols{:}});
                    assignin('base','results_ICC',results_ICC);

                    cd(dir_results);             
                    save results_ICC_2cons.mat results_ICC
   

disp('finished ICC calculation')
end;
end;
cd(box);    
disp('DONE');






% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
axes(hObject);
imshow('logo.png');


