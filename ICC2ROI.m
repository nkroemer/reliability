function varargout = ICC2ROI(varargin)
% ICC2ROI MATLAB code for ICC2ROI.fig
%      ICC2ROI, by itself, creates a new ICC2ROI or raises the existing
%      singleton*.
%
%      H = ICC2ROI returns the handle to a new ICC2ROI or the handle to
%      the existing singleton*.
%
%      ICC2ROI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ICC2ROI.M with the given input arguments.
%
%      ICC2ROI('Property','Value',...) creates a new ICC2ROI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ICC2ROI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ICC2ROI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ICC2ROI

% Last Modified by GUIDE v2.5 02-Jun-2017 14:40:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ICC2ROI_OpeningFcn, ...
                   'gui_OutputFcn',  @ICC2ROI_OutputFcn, ...
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


% --- Executes just before ICC2ROI is made visible.
function ICC2ROI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ICC2ROI (see VARARGIN)

% Choose default command line output for ICC2ROI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ICC2ROI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ICC2ROI_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on button press in def.
function def_Callback(hObject, eventdata, handles)
% hObject    handle to def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

function corr_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to corr_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corr_threshold as text
%        str2double(get(hObject,'String')) returns contents of corr_threshold as a double

% --- Executes during object creation, after setting all properties.
function corr_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corr_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k_Callback(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k as text
%        str2double(get(hObject,'String')) returns contents of k as a double


% --- Executes during object creation, after setting all properties.
function k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in abs.
function abs_Callback(hObject, eventdata, handles)
% hObject    handle to abs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of abs


% --- Executes on button press in cons.
function cons_Callback(hObject, eventdata, handles)
% hObject    handle to cons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cons

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Starting creation of ROI out of ICC maps...');
%% define file seperator 
f = filesep;

study_design = evalin('base','study_design');
contrast_def = evalin('base','contrast_def');
box_path=pwd;

%% get GUI input
abs = get(handles.abs,'Value');
cons = get(handles.cons,'Value');

%% get study design info
results_dir = study_design.results_directory;
stats_path = study_design.stats_path;
load(study_design.subject_list); % loads vp list
runs = str2double(study_design.number_sessions);
stats_filled = sprintf(study_design.stats_directory,1);
two_cons = contrast_def.two_contrasts;
if two_cons == 1
    con = contrast_def.contrast1;
    con1_count = contrast_def.number_contrast1;
    con2_count = contrast_def.number_contrast2;
    con1 = contrast_def.contrast1;
    con2 = contrast_def.contrast1;
else
    con = contrast_def.contrast;
end;

%% reslice atlas and move to results directory
disp('...reslicing atlas...');
cd('atlas')
atlas_dir = pwd;
atlas_name = 'atlas.nii';
atlas_compl = sprintf('%s%s%s%s',atlas_dir,f,f,atlas_name);
    matlabbatch{1}.spm.spatial.coreg.write.ref = {sprintf('%s%s%s%s%s%s%s%s%s%s,1',stats_path,f,f,vp{1},f,f,stats_filled,f,f,con)};
    matlabbatch{1}.spm.spatial.coreg.write.source = {sprintf('%s,1',atlas_compl)};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

    % save batch
    save('batch_reslice_atlas','matlabbatch');

    % run batch
    spm_jobman('serial',matlabbatch);

% load resliced atlas
atlas_r = sprintf('%s%s%sr%s',atlas_dir,f,f,atlas_name);
movefile (atlas_r,results_dir,'f');
cd(results_dir);
atlas = load_nii(sprintf('r%s',atlas_name));
atlas = atlas.img;

% load atlas labels
cd(atlas_dir);
load('labels.mat');

% get GUI input
k = str2double(get(handles.k,'String'));
r_thres = str2double(get(handles.corr_threshold,'String'));

% loads template image to generate image dimensions
temp_img = sprintf('%s%s%s%s%s%s%s%s%s%s',stats_path,f,f,vp{1},f,f,stats_filled,f,f,con); % Beispielbild Kontrast
temp_img=load_nii(temp_img);
dim = size(temp_img.img);
x = dim(1);
y = dim(2);
z = dim(3);
    
%% create ROI
%loads map
cd(results_dir);
if two_cons == 0
if abs == 1
    ICC_map = load_nii('ICC_abs.nii');
    % create correlation matrix   
    ICC_roi = zeros(x,y,z);

    for ind_x = 1:x
        fprintf('...voxel x = %d...\n',ind_x);
        for ind_y = 1:y
            for ind_z = 1:z
                if ICC_map.img(ind_x,ind_y,ind_z) >= r_thres      
                    ICC_roi(ind_x,ind_y,ind_z) = 1;
                else
                    ICC_roi(ind_x,ind_y,ind_z) = 0;
                end;
            end;
        end;
    end;

    % apply cluster threshold
    fprintf('..apply cluster threshold k = %d...\n',k);
    cd(box_path);
    [Y,XYZ,C] = fmri_clust_filt(ICC_roi,k);

    % save ICC_ROIs
    cd(study_design.results_directory) 
    cluster_name = sprintf('cluster_%d_%0.2f_ICCabs.mat',k,r_thres);
    save(cluster_name,'C');
    target_img = temp_img;
    target_img.fileprefix = sprintf('ICCabs_ROI_%d_%0.2f.nii',k,r_thres);
    target_img.img = Y;
    save_nii(target_img,target_img.fileprefix);  

    % load reliable clusters 
    nr_clusters = length(C);

    for i = 1:nr_clusters
        fprintf('...label cluster %d...\n',i);
        vec_voxels = C(i).XYZ;
        labels_voxels = cell(length(vec_voxels),1);
        for j = 1:length(vec_voxels)
            x = vec_voxels(1,j);
            y = vec_voxels(2,j);
            z = vec_voxels(3,j);
            lab = atlas(x,y,z);
            if lab ~= 0
                labels_voxels{j,1} = Labels{lab};
            else
                labels_voxels{j,1} = 'no label';
            end;
        end;
        C(i).labels = labels_voxels;
                % find main region in cluster via most frequent label
        lab_neu = {};
        for i_lab = 1:length(labels_voxels)
            if i_lab == 1
                lab_neu{1,1} = labels_voxels{1,1};
            end;
            for j = 1:length(lab_neu)
                tmp = 0;
                if strcmp(lab_neu{j,1},labels_voxels{i_lab,1}) == 1        
                   break
                else
                   tmp = i_lab;
                end;
            end;
            if tmp ~= 0
                lab_neu{end+1,1} = labels_voxels{tmp};
            end;
        end;

        for k_lab = 1:length(lab_neu)
            count = 0;
            for l = 1:length(labels_voxels)
                if strcmp(lab_neu{k_lab,1},labels_voxels{l,1}) == 1
                    count = count+1;
                end;
            end;
            lab_neu{k_lab,2} = count;
        end;

        numbers = [];
        for m = 1:length(lab_neu(:,2))
            numbers(m,1) = lab_neu{m,2};
        end;
        [numbers_sort,sort_idx] = sort(numbers);
        labels_sorted = lab_neu(sort_idx,1);
        if strcmp(labels_sorted(end,1),'no label') && length(labels_sorted)>1
            main_label = labels_sorted(end-1,1);
        else
            main_label = labels_sorted(end,1);
        end;
        
        C(i).main_label = main_label;
        clearvars vec_voxels labels_voxels
    end;

% conversion to MNI coordinates
ex4d = dir('4D*.nii');
nifti_header = load_untouch_header_only(ex4d(1).name);

transformation_matrix = [nifti_header.hist.srow_x;nifti_header.hist.srow_y;nifti_header.hist.srow_z];

for ind_clust = 1:length(C)
    for coordinates = 1:length(C(ind_clust).XYZ)
        mni_coordinates(coordinates,1:3) = transformation_matrix *[C(ind_clust).XYZ(1,coordinates) C(ind_clust).XYZ(2,coordinates) C(ind_clust).XYZ(3,coordinates) 1]';
    end;
        mni_coordinates = mni_coordinates';
        C(ind_clust).mni = mni_coordinates;
        clearvars mni_coordinates
end;


    cluster_name = sprintf('cluster_%d_%0.2f_ICCabs.mat',k,r_thres);
    save(cluster_name,'C');
% calculate activation estimate for each cluster and each subject
cd(results_dir);
nr_cluster = length(C);

%collect cluster coordinates
for count_cl = 1:nr_cluster
    name = strtok(C(count_cl).main_label{1});
    name(regexp(name,'-'))=[];
    eval(sprintf('%s_xyz = C(count_cl).XYZ;',name));
end;

%load 4D images
for count_runs = 1:runs
    file = sprintf('4D_%d.nii',count_runs);
    eval(sprintf('FourD_%d = load_nii(file);',count_runs));
end;

%extract and average activation in cluster
for count_runs = 1:runs
    for count_cl = 1:nr_cluster
        fprintf('...extracting data for cluster %d...\n',count_cl)
        name = strtok(C(count_cl).main_label{1});
        name(regexp(name,'-'))=[];
        eval(sprintf('vox = size(%s_xyz,2);',name));
        sum_vox = [];
        for ind_subj = 1:length(vp)
        for ix = 1:vox
            eval(sprintf('sum_vox = [sum_vox;FourD_%d.img(%s_xyz(1,ix),%s_xyz(2,ix),%s_xyz(3,ix),ind_subj)];',count_runs,name,name,name));
        end;
        BOLD = mean(sum_vox);
        eval(sprintf('BOLD_%s{ind_subj,1} = vp{ind_subj};',name));
        eval(sprintf('BOLD_%s{ind_subj,1+count_runs} = BOLD;',name));        
        end;
        eval(sprintf('save BOLD_%s_ICCabsROI_%d_%0.2f.mat BOLD_%s;',name,k,r_thres,name));
    end;
end;        
end;
if cons == 1
    ICC_map2 = load_nii('ICC_con.nii');
    % create correlation matrix   
    ICC_roi = zeros(x,y,z);   
    for ind_x = 1:x
        fprintf('...voxel x = %d...\n',ind_x)
        for ind_y = 1:y
            for ind_z = 1:z
                if ICC_map2.img(ind_x,ind_y,ind_z) >= r_thres      
                    ICC_roi(ind_x,ind_y,ind_z) = 1;
                else
                    ICC_roi(ind_x,ind_y,ind_z) = 0;
                end;
            end;
        end;
    end;

    % apply cluster threshold
    fprintf('..apply cluster threshold k = %d...\n',k);
    cd(box_path);
    [Y,XYZ,C] = fmri_clust_filt(ICC_roi,k);


    % save Corr_ROIs
    cd(study_design.results_directory) 
    cluster_name = sprintf('cluster_%d_%0.2f_ICCcon.mat',k,r_thres);
    save(cluster_name,'C');
    target_img = temp_img;
    target_img.fileprefix = sprintf('ICCcon_ROI_%d_%0.2f.nii',k,r_thres);
    target_img.img = Y;
    save_nii(target_img,target_img.fileprefix);  



    % load reliable clusters 
    nr_clusters = length(C);

    for i = 1:nr_clusters
        fprintf('...label cluster %d...',i);
        vec_voxels = C(i).XYZ;
        labels_voxels = cell(length(vec_voxels),1);
        for j = 1:length(vec_voxels)
            x = vec_voxels(1,j);
            y = vec_voxels(2,j);
            z = vec_voxels(3,j);
            lab = atlas(x,y,z);
            if lab ~= 0
                labels_voxels{j,1} = Labels{lab};
            else
                labels_voxels{j,1} = 'no label';
            end;
        end;
        C(i).labels = labels_voxels;
        % find main region in cluster via most frequent label
        lab_neu = {};
        for i_lab = 1:length(labels_voxels)
            if i_lab == 1
                lab_neu{1,1} = labels_voxels{1,1};
            end;
            for j = 1:length(lab_neu)
                tmp = 0;
                if strcmp(lab_neu{j,1},labels_voxels{i_lab,1}) == 1        
                   break
                else
                   tmp = i_lab;
                end;
            end;
            if tmp ~= 0
                lab_neu{end+1,1} = labels_voxels{tmp};
            end;
        end;

        for kc = 1:length(lab_neu)
            count = 0;
            for l = 1:length(labels_voxels)
                if strcmp(lab_neu{kc,1},labels_voxels{l,1}) == 1
                    count = count+1;
                end;
            end;
            lab_neu{kc,2} = count;
        end;

        numbers = [];
        for m = 1:length(lab_neu(:,2))
            numbers(m,1) = lab_neu{m,2};
        end;
        [numbers_sort,sort_idx] = sort(numbers);
        labels_sorted = lab_neu(sort_idx,1);
        if strcmp(labels_sorted(end,1),'no label') && length(labels_sorted)>1
            main_label = labels_sorted(end-1,1);
        else
            main_label = labels_sorted(end,1);
        end;
        
        C(i).main_label = main_label;
        clearvars vec_voxels labels_voxels
    end;

% conversion to MNI coordinates
ex4d = dir('4D*.nii');
nifti_header = load_untouch_header_only(ex4d(1).name);

transformation_matrix = [nifti_header.hist.srow_x;nifti_header.hist.srow_y;nifti_header.hist.srow_z];

for ind_clust = 1:length(C)
    for coordinates = 1:length(C(ind_clust).XYZ)
        mni_coordinates(coordinates,1:3) = transformation_matrix *[C(ind_clust).XYZ(1,coordinates) C(ind_clust).XYZ(2,coordinates) C(ind_clust).XYZ(3,coordinates) 1]';
    end;
        mni_coordinates = mni_coordinates';
        C(ind_clust).mni = mni_coordinates;
        clearvars mni_coordinates
end;


    cluster_name = sprintf('cluster_%d_%0.2f_ICCcon.mat',k,r_thres);
    save(cluster_name,'C');    
    
% calculate activation estimate for each cluster and each subject
cd(results_dir);
nr_cluster = length(C);

%collect cluster coordinates
for count_cl = 1:nr_cluster
    name = strtok(C(count_cl).main_label{1});
    name(regexp(name,'-'))=[];
    eval(sprintf('%s_xyz = C(count_cl).XYZ;',name));
end;

%load 4D images
for count_runs = 1:runs
    file = sprintf('4D_%d.nii',count_runs);
    eval(sprintf('FourD_%d = load_nii(file);',count_runs));
end;

%extract and average activation in cluster
for count_runs = 1:runs
    for count_cl = 1:nr_cluster
        fprintf('...extracting data for cluster %d...\n',count_cl)
        name = strtok(C(count_cl).main_label{1});
        name(regexp(name,'-'))=[];
        eval(sprintf('vox = size(%s_xyz,2);',name));
        sum_vox = [];
        for ind_subj = 1:length(vp)
        for x = 1:vox
            eval(sprintf('sum_vox = [sum_vox;FourD_%d.img(%s_xyz(1,x),%s_xyz(2,x),%s_xyz(3,x),ind_subj)];',count_runs,name,name,name));
        end;
        BOLD = mean(sum_vox);
        eval(sprintf('BOLD_%s{ind_subj,1} = vp{ind_subj};',name));
        eval(sprintf('BOLD_%s{ind_subj,1+count_runs} = BOLD;',name));        
        end;
        eval(sprintf('save BOLD_%s_ICCconROI_%d_%0.2f.mat BOLD_%s;',name,k,r_thres,name));
    end;
end;    
end;
else
for i_con = 1:2
if abs == 1
    eval(sprintf('con_number = con%d_count;',i_con)); 
    contrast = sprintf('con%d',con_number);
    ICC_map = load_nii(sprintf('ICC_abs_%s.nii',contrast));
    % create correlation matrix   
    ICC_roi = zeros(x,y,z);

    for ind_x = 1:x
        fprintf('...voxel x = %d...\n',ind_x);
        for ind_y = 1:y
            for ind_z = 1:z
                if ICC_map.img(ind_x,ind_y,ind_z) >= r_thres      
                    ICC_roi(ind_x,ind_y,ind_z) = 1;
                else
                    ICC_roi(ind_x,ind_y,ind_z) = 0;
                end;
            end;
        end;
    end;

    % apply cluster threshold
    fprintf('..apply cluster threshold k = %d...\n',k);
    cd(box_path);
    [Y,XYZ,C] = fmri_clust_filt(ICC_roi,k);


    % save Corr_ROIs
    cd(study_design.results_directory) 
    cluster_name = sprintf('cluster_%d_%0.2f_ICCabs_%s.mat',k,r_thres,contrast);
    save(cluster_name,'C');
    target_img = temp_img;
    target_img.fileprefix = sprintf('ICCabs_%s_ROI_%d_%0.2f.nii',contrast,k,r_thres);
    target_img.img = Y;
    save_nii(target_img,target_img.fileprefix);  



    % load reliable clusters 
    nr_clusters = length(C);

    for i = 1:nr_clusters
        fprintf('...label cluster %d...\n',i);
        vec_voxels = C(i).XYZ;
        labels_voxels = cell(length(vec_voxels),1);
        for j = 1:length(vec_voxels)
            x = vec_voxels(1,j);
            y = vec_voxels(2,j);
            z = vec_voxels(3,j);
            lab = atlas(x,y,z);
            if lab ~= 0
                labels_voxels{j,1} = Labels{lab};
            else
                labels_voxels{j,1} = 'no label';
            end;
        end;
        C(i).labels = labels_voxels;
                % find main region in cluster via most frequent label
        lab_neu = {};
        for i_lab = 1:length(labels_voxels)
            if i_lab == 1
                lab_neu{1,1} = labels_voxels{1,1};
            end;
            for j = 1:length(lab_neu)
                tmp = 0;
                if strcmp(lab_neu{j,1},labels_voxels{i_lab,1}) == 1        
                   break
                else
                   tmp = i_lab;
                end;
            end;
            if tmp ~= 0
                lab_neu{end+1,1} = labels_voxels{tmp};
            end;
        end;

        for k_lab = 1:length(lab_neu)
            count = 0;
            for l = 1:length(labels_voxels)
                if strcmp(lab_neu{k_lab,1},labels_voxels{l,1}) == 1
                    count = count+1;
                end;
            end;
            lab_neu{k_lab,2} = count;
        end;

        numbers = [];
        for m = 1:length(lab_neu(:,2))
            numbers(m,1) = lab_neu{m,2};
        end;
        [numbers_sort,sort_idx] = sort(numbers);
        labels_sorted = lab_neu(sort_idx,1);
        if strcmp(labels_sorted(end,1),'no label') && length(labels_sorted)>1
            main_label = labels_sorted(end-1,1);
        else
            main_label = labels_sorted(end,1);
        end;
        
        C(i).main_label = main_label;
        clearvars vec_voxels labels_voxels
    end;

% conversion to MNI coordinates
ex4d = dir('4D*.nii');
nifti_header = load_untouch_header_only(ex4d(1).name);

transformation_matrix = [nifti_header.hist.srow_x;nifti_header.hist.srow_y;nifti_header.hist.srow_z];

for ind_clust = 1:length(C)
    for coordinates = 1:length(C(ind_clust).XYZ)
        mni_coordinates(coordinates,1:3) = transformation_matrix *[C(ind_clust).XYZ(1,coordinates) C(ind_clust).XYZ(2,coordinates) C(ind_clust).XYZ(3,coordinates) 1]';
    end;
        mni_coordinates = mni_coordinates';
        C(ind_clust).mni = mni_coordinates;
        clearvars mni_coordinates
end;


    cluster_name = sprintf('cluster_%d_%0.2f_ICCabs_%s.mat',k,r_thres,contrast);
    save(cluster_name,'C');
% calculate activation estimate for each cluster and each subject
cd(results_dir);
nr_cluster = length(C);

%collect cluster coordinates
for count_cl = 1:nr_cluster
    name = strtok(C(count_cl).main_label{1});
    name(regexp(name,'-'))=[];
    eval(sprintf('%s_xyz = C(count_cl).XYZ;',name));
end;

%load 4D images
eval(sprintf('cont = con%d;',i_con));
for count_runs = 1:runs
    file = sprintf('4D_%s_%d.nii',cont,count_runs);
    eval(sprintf('FourD_%d = load_nii(file);',count_runs));
end;

%extract and average activation in cluster
for count_runs = 1:runs
    for count_cl = 1:nr_cluster
        fprintf('...extracting data for cluster %d...\n',count_cl)
        name = strtok(C(count_cl).main_label{1});
        name(regexp(name,'-'))=[];
        eval(sprintf('vox = size(%s_xyz,2);',name));
        sum_vox = [];
        for ind_subj = 1:length(vp)
        for x = 1:vox
            eval(sprintf('sum_vox = [sum_vox;FourD_%d.img(%s_xyz(1,x),%s_xyz(2,x),%s_xyz(3,x),ind_subj)];',count_runs,name,name,name));
        end;
        BOLD = mean(sum_vox);
        eval(sprintf('BOLD_%s{ind_subj,1} = vp{ind_subj};',name));
        eval(sprintf('BOLD_%s{ind_subj,1+count_runs} = BOLD;',name));        
        end;
        eval(sprintf('save BOLD_%s_ICCabsROI_%s_%d_%0.2f.mat BOLD_%s;',name,contrast,k,r_thres,name));
    end;
end;        
end;
if cons == 1
    eval(sprintf('con_number = con%d_count;',i_con)); 
    contrast = sprintf('con%d',con_number);
    ICC_map2 = load_nii(sprintf('ICC_con_%s.nii',contrast));
    % create correlation matrix   
    ICC_roi = zeros(x,y,z);   
    x = size(ICC_map2.img,1);
    y = size(ICC_map2.img,2);
    z = size(ICC_map2.img,3);
    
    for ind_x = 1:x
        fprintf('...voxel x = %d...\n',ind_x)
        for ind_y = 1:y
            for ind_z = 1:z
                if ICC_map2.img(ind_x,ind_y,ind_z) >= r_thres      
                    ICC_roi(ind_x,ind_y,ind_z) = 1;
                else
                    ICC_roi(ind_x,ind_y,ind_z) = 0;
                end;
            end;
        end;
    end;

    % apply cluster threshold
    fprintf('..apply cluster threshold k = %d...\n',k);
    cd(box_path);
    [Y,XYZ,C] = fmri_clust_filt(ICC_roi,k);


    % save Corr_ROIs
    cd(study_design.results_directory) 
    cluster_name = sprintf('cluster_%d_%0.2f_ICCcon_%s.mat',k,r_thres,contrast);
    save(cluster_name,'C');
    target_img = temp_img;
    target_img.fileprefix = sprintf('ICCcon_%s_ROI_%d_%0.2f.nii',contrast,k,r_thres);
    target_img.img = Y;
    save_nii(target_img,target_img.fileprefix);  



    % load reliable clusters 
    nr_clusters = length(C);

    for i = 1:nr_clusters
        fprintf('...label cluster %d...',i);
        vec_voxels = C(i).XYZ;
        labels_voxels = cell(length(vec_voxels),1);
        for j = 1:length(vec_voxels)
            x = vec_voxels(1,j);
            y = vec_voxels(2,j);
            z = vec_voxels(3,j);
            lab = atlas(x,y,z);
            if lab ~= 0
                labels_voxels{j,1} = Labels{lab};
            else
                labels_voxels{j,1} = 'no label';
            end;
        end;
        C(i).labels = labels_voxels;
        % find main region in cluster via most frequent label
        lab_neu = {};
        for i_lab = 1:length(labels_voxels)
            if i_lab == 1
                lab_neu{1,1} = labels_voxels{1,1};
            end;
            for j = 1:length(lab_neu)
                tmp = 0;
                if strcmp(lab_neu{j,1},labels_voxels{i_lab,1}) == 1        
                   break
                else
                   tmp = i_lab;
                end;
            end;
            if tmp ~= 0
                lab_neu{end+1,1} = labels_voxels{tmp};
            end;
        end;

        for kc = 1:length(lab_neu)
            count = 0;
            for l = 1:length(labels_voxels)
                if strcmp(lab_neu{kc,1},labels_voxels{l,1}) == 1
                    count = count+1;
                end;
            end;
            lab_neu{kc,2} = count;
        end;

        numbers = [];
        for m = 1:length(lab_neu(:,2))
            numbers(m,1) = lab_neu{m,2};
        end;
        [numbers_sort,sort_idx] = sort(numbers);
        labels_sorted = lab_neu(sort_idx,1);
        if strcmp(labels_sorted(end,1),'no label') && length(labels_sorted)>1
            main_label = labels_sorted(end-1,1);
        else
            main_label = labels_sorted(end,1);
        end;
        
        C(i).main_label = main_label;
        clearvars vec_voxels labels_voxels
    end;

% conversion to MNI coordinates
ex4d = dir('4D*.nii');
nifti_header = load_untouch_header_only(ex4d(1).name);

transformation_matrix = [nifti_header.hist.srow_x;nifti_header.hist.srow_y;nifti_header.hist.srow_z];

for ind_clust = 1:length(C)
    for coordinates = 1:length(C(ind_clust).XYZ)
        mni_coordinates(coordinates,1:3) = transformation_matrix *[C(ind_clust).XYZ(1,coordinates) C(ind_clust).XYZ(2,coordinates) C(ind_clust).XYZ(3,coordinates) 1]';
    end;
        mni_coordinates = mni_coordinates';
        C(ind_clust).mni = mni_coordinates;
        clearvars mni_coordinates
end;


    cluster_name = sprintf('cluster_%d_%0.2f_ICCcon_%s.mat',k,r_thres,contrast);
    save(cluster_name,'C');    
    
% calculate activation estimate for each cluster and each subject
cd(results_dir);
nr_cluster = length(C);

%collect cluster coordinates

for count_cl = 1:nr_cluster
    name = strtok(C(count_cl).main_label{1});
    name(regexp(name,'-'))=[];
    eval(sprintf('%s_xyz = C(count_cl).XYZ;',name));
end;

eval(sprintf('cont = con%d;',i_con));

%load 4D images
for count_runs = 1:runs
    file = sprintf('4D_%s_%d.nii',cont,count_runs);
    eval(sprintf('FourD_%d = load_nii(file);',count_runs));
end;

%extract and average activation in cluster
for count_runs = 1:runs
    for count_cl = 1:nr_cluster
        fprintf('...extracting data for cluster %d...\n',count_cl)
        name = strtok(C(count_cl).main_label{1});
        name(regexp(name,'-'))=[];
        eval(sprintf('vox = size(%s_xyz,2);',name));
        sum_vox = [];
        for ind_subj = 1:length(vp)
        for x = 1:vox
            eval(sprintf('sum_vox = [sum_vox;FourD_%d.img(%s_xyz(1,x),%s_xyz(2,x),%s_xyz(3,x),ind_subj)];',count_runs,name,name,name));
        end;
        BOLD = mean(sum_vox);
        eval(sprintf('BOLD_%s{ind_subj,1} = vp{ind_subj};',name));
        eval(sprintf('BOLD_%s{ind_subj,1+count_runs} = BOLD;',name));        
        end;
        eval(sprintf('save BOLD_%s_ICCconROI_%s_%d_%0.2f.mat BOLD_%s;',name,contrast,k,r_thres,name));
    end;
end;    
end;    
end;
end;



cd(box_path);
disp('DONE');


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
axes(hObject);
imshow('logo.png');