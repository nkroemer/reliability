function varargout = CorrMaps(varargin)
% CORRMAPS MATLAB code for CorrMaps.fig
%      CORRMAPS, by itself, creates a new CORRMAPS or raises the existing
%      singleton*.
%
%      H = CORRMAPS returns the handle to a new CORRMAPS or the handle to
%      the existing singleton*.
%
%      CORRMAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRMAPS.M with the given input arguments.
%
%      CORRMAPS('Property','Value',...) creates a new CORRMAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CorrMaps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CorrMaps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CorrMaps

% Last Modified by GUIDE v2.5 02-Jun-2017 14:35:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrMaps_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrMaps_OutputFcn, ...
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


% --- Executes just before CorrMaps is made visible.
function CorrMaps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CorrMaps (see VARARGIN)

% Choose default command line output for CorrMaps
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CorrMaps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CorrMaps_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on button press in split.
function split_Callback(hObject, eventdata, handles)
% hObject    handle to split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of split


function prefix_Callback(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix as text
%        str2double(get(hObject,'String')) returns contents of prefix as a double


% --- Executes during object creation, after setting all properties.
function prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pearson.
function pearson_Callback(hObject, eventdata, handles)
% hObject    handle to pearson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pearson


% --- Executes on button press in spearman.
function spearman_Callback(hObject, eventdata, handles)
% hObject    handle to spearman (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spearman

% --- Executes on button press in con_def.
function con_def_Callback(hObject, eventdata, handles)
% hObject    handle to con_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Starting calculation of correlation maps...');
%% define file seperator 
f = filesep;
%% set parameters
%get study design information
study_design=evalin('base','study_design');
contrast_def=evalin('base','contrast_def');

runs=str2double(study_design.number_sessions); 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
stats=study_design.stats_directory;
path=study_design.stats_path;
dir_results=study_design.results_directory;
if runs == 1
    single_run = str2double(study_design.identifier_session);
end;
par = study_design.number_parametric;

%get GUI input
pear = get(handles.pearson,'value');
spea = get(handles.spearman,'value');
split = get(handles.split,'value');
name = get(handles.prefix,'String');

%get contrast information
two_cons = contrast_def.two_contrasts;
if two_cons == 0
    con=contrast_def.contrast;
    disp('...loads image dimensions..');
    stats_temp =sprintf(stats,1);
    temp_img = sprintf('%s%s%s%s%s%s%s%s%s%s',path,f,f,vp{1},f,f,stats_temp,f,f,con);
    temp_img=load_nii(temp_img);
    dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);    
else
    con1=contrast_def.contrast1;
    con1_count=contrast_def.number_contrast1;
    con2=contrast_def.contrast2;
    con2_count=contrast_def.number_contrast2;
    
    disp('...loads image dimensions..');
    if runs == 1
        stats_temp = sprintf(stats,single_run);
    else
        stats_temp =sprintf(stats,1);
    end;
    temp_img = sprintf('%s%s%s%s%s%s%s%s%s%s',path,f,f,vp{1},f,f,stats_temp,f,f,con1);
    temp_img=load_nii(temp_img);
    dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);
end;

boxpath=pwd;
cd (dir_results); 

% initiate summary
cols = {};
summary = [];

%% create correlation maps
%% 'normal' design
if split == 0 && two_cons == 0 && runs > 1
    
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
        if i_run+i_sec <= runs
            fprintf('...creates correlation maps for session %d and session %d...\n',i_run,i_run+i_sec);
            %load 4D images
            file1 = sprintf('4D_%d.nii',i_run);
            one = load_nii(file1);
            one = one.img;
            one (~one) = nan;

            file2 = sprintf('4D_%d.nii',i_run+i_sec);
            second = load_nii(file2);
            second = second.img;
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
    file = sprintf('%s_pear_%d_%d.nii',name,i_run,i_run+i_sec);
    target_img.fileprefix = file;
    target_img.img = r_vec_pear_1_2;
    save_nii(target_img,target_img.fileprefix);  
         
    target_img = temp_img;
    file = sprintf('z_%s_pear_%d_%d.nii',name,i_run,i_run+i_sec);
    target_img.fileprefix = file;
    target_img.img = z_r_vec_pear_1_2;
    save_nii(target_img,target_img.fileprefix);                     

    target_img = temp_img;
    file = sprintf('%s_spea_%d_%d.nii',name,i_run,i_run+i_sec);
    target_img.fileprefix = file;
    target_img.img = r_vec_spea_1_2;
    save_nii(target_img,target_img.fileprefix);
                    
    target_img = temp_img;
    file = sprintf('z_%s_spea_%d_%d.nii',name,i_run,i_run+i_sec);
    target_img.fileprefix = file;
    target_img.img = z_r_vec_spea_1_2;
    save_nii(target_img,target_img.fileprefix);                      
    
    %compute z mean and inverse Fisher's transformation
    mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
    mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
    summary(1,end+1)=mean_r;
    cols{1,end+1} = sprintf('mean_pear_%d_%d',i_run,i_run+i_sec);
    
    mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
    mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
    summary(1,end+1)=mean_r;
    cols{1,end+1} = sprintf('mean_spea_%d_%d',i_run,i_run+i_sec);    

        end; 
        end;
    end;
 
if par > 0    
    for i_par = 1:par
    for i_run = 1:runs
       for i_sec = 1:runs-i_run

        if i_run+1 <= runs
            fprintf('...creates correlation maps for parametric modulator %d in session %d and session %d...\n',i_par,i_run,i_run+i_sec);
            %load 4D images
            file1 = sprintf('4D_par%d_%d.nii',i_par,i_run);
            one = load_nii(file1);
            one = one.img;
            one (~one) = nan;

            file2 = sprintf('4D_par%d_%d.nii',i_par,i_run+i_sec);
            second = load_nii(file2);
            second = second.img;
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
    file = sprintf('%s_pear_par%d_%d_%d.nii',name,i_par,i_run,i_run+i_sec);
    target_img.fileprefix = file;
    target_img.img = r_vec_pear_1_2;
    save_nii(target_img,target_img.fileprefix);  
         
    target_img = temp_img;
    file = sprintf('z_%s_pear_par%d_%d_%d.nii',name,i_par,i_run,i_run+i_sec);
    target_img.fileprefix = file;
    target_img.img = z_r_vec_pear_1_2;
    save_nii(target_img,target_img.fileprefix);                     

    target_img = temp_img;
    file = sprintf('%s_spea_par%d_%d_%d.nii',name,i_par,i_run,i_run+i_sec);
    target_img.fileprefix = file;
    target_img.img = r_vec_spea_1_2;
    save_nii(target_img,target_img.fileprefix);
                    
    target_img = temp_img;
    file = sprintf('z_%s_spea_par%d_%d_%d.nii',name,i_par,i_run,i_run+i_sec);
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
end;
    
    
%% based on split half
elseif split == 1
    for i = 1:runs
        fprintf('...creates correlation maps for splitted session %d...\n',i);
        if runs == 1 
            img1=sprintf('4D_split1_%d.nii',single_run);
            one = load_nii(img1);
            one = one.img;
            one (~one) = nan;

            img2=sprintf('4D_split2_%d.nii',single_run);
            second = load_nii(img2);
            second = second.img;
            second (~second) = nan;
        else
            img1=sprintf('4D_split1_%d.nii',i);
            one = load_nii(img1);
            one = one.img;
            one (~one) = nan;

            img2=sprintf('4D_split2_%d.nii',i);
            second = load_nii(img2);
            second = second.img;
            second (~second) = nan;
        end;

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
        
        if runs == 1
        % save correlation map
        target_img = temp_img;
        file = sprintf('%s_pear_%d_split.nii',name,single_run);
        target_img.fileprefix = file;
        target_img.img = r_vec_pear_1_2;
        save_nii(target_img,target_img.fileprefix);   

        target_img = temp_img;
        file = sprintf('%s_spea_%d_split.nii',name,single_run);
        target_img.fileprefix = file;
        target_img.img = r_vec_spea_1_2;
        save_nii(target_img,target_img.fileprefix);

        target_img = temp_img;
        file = sprintf('z_%s_pear_%d_split.nii',name,single_run);
        target_img.fileprefix = file;
        target_img.img = z_r_vec_pear_1_2;
        save_nii(target_img,target_img.fileprefix);   

        target_img = temp_img;
        file = sprintf('z_%s_spea_%d_split.nii',name,single_run);
        target_img.fileprefix = file;
        target_img.img = z_r_vec_spea_1_2;
        save_nii(target_img,target_img.fileprefix);
                    
        %compute z mean and inverse Fisher's transformation
        mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_pear_split_%d',single_run);
        mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_spea_split_%d',single_run);

        else
        % save correlation map
        target_img = temp_img;
        file = sprintf('%s_pear_%d_split.nii',name,i);
        target_img.fileprefix = file;
        target_img.img = r_vec_pear_1_2;
        save_nii(target_img,target_img.fileprefix);   

        target_img = temp_img;
        file = sprintf('%s_spea_%d_split.nii',name,i);
        target_img.fileprefix = file;
        target_img.img = r_vec_spea_1_2;
        save_nii(target_img,target_img.fileprefix);

        target_img = temp_img;
        file = sprintf('z_%s_pear_%d_split.nii',name,i);
        target_img.fileprefix = file;
        target_img.img = z_r_vec_pear_1_2;
        save_nii(target_img,target_img.fileprefix);   

        target_img = temp_img;
        file = sprintf('z_%s_spea_%d_split.nii',name,i);
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
                  
    
    end;  
% comparison of splitted parametric regressor    
    if par > 0
        for ind_par = 1:par
            for i = 1:runs
                fprintf('...creates correlation maps for parametric modulators in session %d...\n',i);

                if runs == 1
                    img1=sprintf('4D_split1_par%d_%d.nii',ind_par,single_run);
                    one = load_nii(img1);
                    one = one.img;
                    one (~one) = nan;

                    img2=sprintf('4D_split2_par%d_%d.nii',ind_par,single_run);
                    second = load_nii(img2);
                    second = second.img;
                    second (~second) = nan;
                else
                    img1=sprintf('4D_split1_par%d_%d.nii',ind_par,i);
                    one = load_nii(img1);
                    one = one.img;
                    one (~one) = nan;

                    img2=sprintf('4D_split2_par%d_%d.nii',ind_par,i);
                    second = load_nii(img2);
                    second = second.img;
                    second (~second) = nan;
                end;
                
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
                
                if runs == 1
                % save correlation map
                target_img = temp_img;
                file = sprintf('%s_pear_par%d_%d_split.nii',name,ind_par,single_run);
                target_img.fileprefix = file;
                target_img.img = r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);   

                target_img = temp_img;
                file = sprintf('%s_spea_par%d_%d_split.nii',name,ind_par,single_run);
                target_img.fileprefix = file;
                target_img.img = r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                target_img = temp_img;
                file = sprintf('z_%s_pear_par%d_%d_split.nii',name,ind_par,single_run);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);   

                target_img = temp_img;
                file = sprintf('z_%s_spea_par%d_%d_split.nii',name,ind_par,single_run);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);
                        
                %compute z mean and inverse Fisher's transformation
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear_split_%d_par%d',single_run,ind_par);
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea_split_%d_par%d',single_run, ind_par);

                else
                % save correlation map
                target_img = temp_img;
                file = sprintf('%s_pear_par%d_%d_split.nii',name,ind_par,i);
                target_img.fileprefix = file;
                target_img.img = r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);   

                target_img = temp_img;
                file = sprintf('%s_spea_par%d_%d_split.nii',name,ind_par,i);
                target_img.fileprefix = file;
                target_img.img = r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                target_img = temp_img;
                file = sprintf('z_%s_pear_par%d_%d_split.nii',name,ind_par,i);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);   

                target_img = temp_img;
                file = sprintf('z_%s_spea_par%d_%d_split.nii',name,ind_par,i);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);
                        
                %compute z mean and inverse Fisher's transformation
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear_split_%d_par%d',i,ind_par);
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea_split_%d_par%d',i, ind_par);

                end;
            end;  
        end;
    end;

%% two contrasts out of one statistic    
elseif two_cons == 1
    % compare contrasts within sessions
    for i_run = 1:runs
        fprintf('...\n compare contrasts of session %d\n...',i_run)
    
        if runs == 1
            one = load_nii(sprintf('4D_%s_%d.nii',con1,single_run));
            one = one.img;
            one (~one) = nan;

            second = load_nii(sprintf('4D_%s_%d.nii',con2,single_run));
            second = second.img;
            second (~second) = nan;
        else
            one = load_nii(sprintf('4D_%s_%d.nii',con1,i_run));
            one = one.img;
            one (~one) = nan;

            second = load_nii(sprintf('4D_%s_%d.nii',con2,i_run));
            second = second.img;
            second (~second) = nan;
        end;
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
        
        if runs == 1
            % save correlation map
            fprintf('...saving maps...\n')
            target_img = temp_img;
            file = sprintf('%s_pear_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_%s_pear_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);                     

            target_img = temp_img;
            file = sprintf('%s_spea_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_%s_spea_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);  

            %compute z mean and inverse Fisher's transformation
            mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_pear_%d_con%d_con%d',single_run,con1_count,con2_count);
            mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_spea_%d_con%d_con%d',single_run,con1_count,con2_count);
          
        else
            % save correlation map
            fprintf('...saving maps...\n')
            target_img = temp_img;
            file = sprintf('%s_pear_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_%s_pear_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);                     

            target_img = temp_img;
            file = sprintf('%s_spea_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_%s_spea_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);  

            %compute z mean and inverse Fisher's transformation
            mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_pear_%d_con%d_con%d',i_run,con1_count,con2_count);
            mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_spea_%d_con%d_con%d',i_run,con1_count,con2_count);

        end;
    end;
    % compare contrasts between sessions
    for i_con = 1:2
        eval(sprintf('con=con%d;',i_con));
        eval(sprintf('con_count = con%d_count;',i_con));
    for i_run = 1:runs
        for count = 1:runs-1
            if i_run+count <= runs
            if runs > 1
                fprintf('...\n compare contrast %s between sessions\n...',con)
                one = load_nii(sprintf('4D_%s_%d.nii',con,i_run));
                one = one.img;
                one (~one) = nan;

                second = load_nii(sprintf('4D_%s_%d.nii',con,i_run+count));
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
                file = sprintf('%s_pear_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);  

                target_img = temp_img;
                file = sprintf('z_%s_pear_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);                     

                target_img = temp_img;
                file = sprintf('%s_spea_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                target_img = temp_img;
                file = sprintf('z_%s_spea_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
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
    for i_run = 1:runs
        fprintf('...\n compare parametric modulators for contrasts of session %d\n...',i_run)
    
        if runs == 1
            one = load_nii(sprintf('4D_%s_%d.nii',con1,single_run));
            one = one.img;
            one (~one) = nan;

            second = load_nii(sprintf('4D_%s_%d.nii',con2,single_run));
            second = second.img;
            second (~second) = nan;
        else
            one = load_nii(sprintf('4D_%s_%d.nii',con1,i_run));
            one = one.img;
            one (~one) = nan;

            second = load_nii(sprintf('4D_%s_%d.nii',con2,i_run));
            second = second.img;
            second (~second) = nan;
        end;
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
        
        if runs == 1
            % save correlation map
            fprintf('...saving maps...\n')
            target_img = temp_img;
            file = sprintf('%s_pear_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_%s_pear_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);                     

            target_img = temp_img;
            file = sprintf('%s_spea_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_%s_spea_%d_con%d_con%d.nii',name,single_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);  

            %compute z mean and inverse Fisher's transformation
            mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_pear_%d_con%d_con%d',single_run,con1_count,con2_count);
            mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_spea_%d_con%d_con%d',single_run,con1_count,con2_count);
          
        else
            % save correlation map
            fprintf('...saving maps...\n')
            target_img = temp_img;
            file = sprintf('%s_pear_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_%s_pear_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_pear_1_2;
            save_nii(target_img,target_img.fileprefix);                     

            target_img = temp_img;
            file = sprintf('%s_spea_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_%s_spea_%d_con%d_con%d.nii',name,i_run,con1_count,con2_count);
            target_img.fileprefix = file;
            target_img.img = z_r_vec_spea_1_2;
            save_nii(target_img,target_img.fileprefix);  

            %compute z mean and inverse Fisher's transformation
            mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_pear_%d_con%d_con%d',i_run,con1_count,con2_count);
            mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_spea_%d_con%d_con%d',i_run,con1_count,con2_count);

        end;
    end;
    % compare contrasts between sessions
    for i_con = 1:2
        eval(sprintf('con=con%d;',i_con));
        eval(sprintf('con_count = con%d_count;',i_con));
    for i_run = 1:runs
        for count = 1:runs-1
            if i_run+count <= runs
            if runs > 1
                fprintf('...\n compare contrast %s between sessions\n...',con)
                one = load_nii(sprintf('4D_%s_%d.nii',con,i_run));
                one = one.img;
                one (~one) = nan;

                second = load_nii(sprintf('4D_%s_%d.nii',con,i_run+count));
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
                file = sprintf('%s_pear_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);  

                target_img = temp_img;
                file = sprintf('z_%s_pear_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = z_r_vec_pear_1_2;
                save_nii(target_img,target_img.fileprefix);                     

                target_img = temp_img;
                file = sprintf('%s_spea_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
                target_img.fileprefix = file;
                target_img.img = r_vec_spea_1_2;
                save_nii(target_img,target_img.fileprefix);

                target_img = temp_img;
                file = sprintf('z_%s_spea_%d_%d_con%d.nii',name,i_run,i_run+count,con_count);
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
end;


results_corr=dataset({summary(1,:),cols{:}});
assignin('base','results_corr',results_corr);
    
cd(dir_results);             
save results_corr.mat results_corr            
disp('...finished correlation maps')
cd(boxpath);


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
axes(hObject);
imshow('logo.png');
