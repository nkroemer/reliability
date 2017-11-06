function varargout = Similarity(varargin)
% SIMILARITY MATLAB code for Similarity.fig
%      SIMILARITY, by itself, creates a new SIMILARITY or raises the existing
%      singleton*.
%
%      H = SIMILARITY returns the handle to a new SIMILARITY or the handle to
%      the existing singleton*.
%
%      SIMILARITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMILARITY.M with the given input arguments.
%
%      SIMILARITY('Property','Value',...) creates a new SIMILARITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Similarity_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Similarity_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Similarity

% Last Modified by GUIDE v2.5 26-Jul-2017 09:47:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Similarity_OpeningFcn, ...
                   'gui_OutputFcn',  @Similarity_OutputFcn, ...
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


% --- Executes just before Similarity is made visible.
function Similarity_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Similarity (see VARARGIN)

% Choose default command line output for Similarity
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Similarity wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Similarity_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
set(handles.name_roi,'TooltipString','Type in name of ROI file WITHOUT suffix');


% --- Executes on button press in design.
function design_Callback(hObject, eventdata, handles)
% hObject    handle to design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
study_design = cellstr(spm_select(1,'mat','load study design'));
load(study_design{1});
assignin('base','study_design',study_design);

% --- Executes on button press in contrast_def.
function contrast_def_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

% --- Executes on button press in use_roi.
function use_roi_Callback(hObject, eventdata, handles)
% hObject    handle to use_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_roi



function name_roi_Callback(hObject, eventdata, handles)
% hObject    handle to name_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name_roi as text
%        str2double(get(hObject,'String')) returns contents of name_roi as a double


% --- Executes during object creation, after setting all properties.
function name_roi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dir_roi.
function dir_roi_Callback(hObject, eventdata, handles)
% hObject    handle to dir_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi_dir = spm_select(1,'dir','choose roi directory');
assignin('base','roi_dir',roi_dir);


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
%% define file seperator 
f = filesep;

disp('Starting calculations of similarity...');
study_design=evalin('base','study_design');
contrast_def = evalin('base','contrast_def');
box_path=evalin('base','box_path');
cd(box_path);

%% get study design info
results_dir = study_design.results_directory;
runs = str2double(study_design.number_sessions);
nr_subj = str2double(study_design.number_subjects);
load(study_design.subject_list);
stats=study_design.stats_directory;
path=study_design.stats_path;

%% get GUI input
split = get(handles.split,'Value');
use_roi = get(handles.use_roi,'Value');
if use_roi == 1
    str = 'ROI';
else
    str = '';
end;

%% load contrast information
two_cons = contrast_def.two_contrasts;
if two_cons == 0
    con=contrast_def.contrast;
    nr_para = study_design.number_parametric;

else
    con=[contrast_def.contrast1 contrast_def.contrast_format];
    con1=contrast_def.contrast1;
    con2=contrast_def.contrast2;
    con1_count=contrast_def.contrast1_number;
    con2_count=contrast_def.contrast2_number;
    %nr_para1 = study_design.number_parametric1;
    %nr_para2 = study_design.number_parametric2;
    
end;

%% load and reslice ROI
if use_roi == 1
    
    disp('...load and reslice ROI...')
    roi_pur=get(handles.name_roi,'String');
    roi_dir=evalin('base','roi_dir');
    cd(roi_dir);
    roi_name = dir(sprintf('%s*',roi_pur));
    if isstruct(roi_name)
        if length(roi_name)==2
            roi_name = roi_name(2).name;
        else
            roi_name = roi_name(1).name;
        end;
    end;
    roi_compl = [roi_dir f roi_name];
    disp('...reslicing ROI...');
    stats_filled = sprintf(stats,1);
    temp = [path f id{1} f stats_filled f con ',1'];
    matlabbatch{1}.spm.spatial.coreg.write.ref = {temp};
    matlabbatch{1}.spm.spatial.coreg.write.source = {sprintf('%s,1',roi_compl)};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    % save batch
    save('batch_reslice_roi','matlabbatch');

    % run batch
    spm_jobman('serial',matlabbatch);
    r_roi=dir(sprintf('r%s*',roi_pur));
    if length(r_roi)==2
        if ~strcmp(roi_dir,results_dir)
        compl1 = [roi_dir f 'r' roi_pur '.img'];
        movefile(compl1,results_dir,'f');
        compl2 = [roi_dir f 'r' roi_pur '.hdr'];
        movefile(compl2,results_dir,'f');
        end;
        cd(results_dir);    
        r_roi = load_nii(sprintf('r%s.img',roi_pur));
        r_roi_ind = r_roi.img==1;
    else
        if ~strcmp(roi_dir,results_dir)
        compl = [roi_dir f  'r' roi_pur '.nii'];
        movefile(compl,results_dir,'f');
        end;
        cd(results_dir);    
        r_roi = load_nii(sprintf('r%s.nii',roi));
        r_roi_ind = r_roi.img==1;        
    end;
end;

% load 4D images
if two_cons == 0 && split == 0
    for ind_run = 1:runs
        fprintf('...load 4D image for run %d...\n',ind_run);
        file = [results_dir f '4D_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file;',ind_run));
        eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d = FourD%d.img;',ind_run,ind_run));
        
        %calculation of mean activation 
        eval(sprintf('dims = size(FourD%d(:,:,:,1));',ind_run));
        x = dims(1);
        y = dims(2);
        z = dims(3);
        eval(sprintf('mean_FourD%d = zeros(x,y,z);',ind_run));
        fprintf('...generate mean activation per voxel for run %d...\n',ind_run); 
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z 
                    eval(sprintf('mean_FourD%d(ind_x,ind_y,ind_z) = mean(FourD%d(ind_x,ind_y,ind_z,:));',ind_run,ind_run));
                end;
            end;
        end;
        eval(sprintf('mean_FourD%d = mean_FourD%d(~isnan(mean_FourD%d));',ind_run,ind_run,ind_run));
    end;
    % load parametric modulator
    if nr_para > 0
        for ind_para = 1:nr_para
            for ind_run = 1:runs
                fprintf('...load 4D parametric image for run %d...\n',ind_run);
                file = [results_dir f '4D_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d_par = file;',ind_run));
                eval(sprintf('FourD%d_par = load_nii(FourD%d_par);',ind_run,ind_run));
                eval(sprintf('FourD%d_par%d = FourD%d_par.img;',ind_run,ind_para,ind_run));

                %calculation of mean activation 
                eval(sprintf('dims = size(FourD%d_par%d(:,:,:,1));',ind_run,ind_para));
                x = dims(1);
                y = dims(2);
                z = dims(3);
                eval(sprintf('mean_FourD%d_par = zeros(x,y,z);',ind_run));
                fprintf('...generate mean activation per voxel for run %d...\n',ind_run); 
                for ind_x = 1:x
                    for ind_y = 1:y
                        for ind_z = 1:z 
                            eval(sprintf('mean_FourD%d_par(ind_x,ind_y,ind_z) = mean(FourD%d_par%d(ind_x,ind_y,ind_z,:));',ind_run,ind_run,ind_para));
                        end;
                    end;
                end;
                eval(sprintf('mean_FourD%d_par%d = mean_FourD%d_par(~isnan(mean_FourD%d_par));',ind_run,ind_para,ind_run,ind_run));
            end;
        end;
    end;
elseif two_cons == 1
     for ind_run = 1:runs
        fprintf('...load 4D image for run %d and contrast %s...\n',ind_run,con1);
        file1 = [results_dir f '4D_' con1 '_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file1;',ind_run));
        eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d_%d = FourD%d.img;',ind_run,con1_count,ind_run));
        eval(sprintf('dims = size(FourD%d_%d(:,:,:,1));',ind_run,con1_count));
        x = dims(1);
        y = dims(2);
        z = dims(3);
        eval(sprintf('mean_FourD%d_%d = zeros(x,y,z);',ind_run,con1_count));
        fprintf('...generate mean activation per voxel for run %d and contrast %s...\n',ind_run,con1); 
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z 
                    eval(sprintf('mean_FourD%d_%d(ind_x,ind_y,ind_z) = mean(FourD%d_%d(ind_x,ind_y,ind_z,:));',ind_run,con1_count,ind_run,con1_count));
                end;
            end;
        end;
        eval(sprintf('mean_FourD%d_%d = mean_FourD%d_%d(~isnan(mean_FourD%d_%d));',ind_run,con1_count,ind_run,con1_count,ind_run,con1_count));

        fprintf('...load 4D image for run %d and contrast %s...\n',ind_run,con2);
        file2 = [results_dir f '4D_' con2 '_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file2;',ind_run));
        eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d_%d = FourD%d.img;',ind_run,con2_count,ind_run));
        eval(sprintf('dims = size(FourD%d_%d(:,:,:,1));',ind_run,con2_count));
        x = dims(1);
        y = dims(2);
        z = dims(3);
        fprintf('...generate mean activation per voxel for run %d and contrast %s...\n',ind_run,con2); 
        eval(sprintf('mean_FourD%d_%d = zeros(x,y,z);',ind_run,con2_count));
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z 
                    eval(sprintf('mean_FourD%d_%d(ind_x,ind_y,ind_z) = mean(FourD%d_%d(ind_x,ind_y,ind_z,:));',ind_run,con2_count,ind_run,con2_count));
                end;
            end;
        end;
        eval(sprintf('mean_FourD%d_%d = mean_FourD%d_%d(~isnan(mean_FourD%d_%d));',ind_run,con2_count,ind_run,con2_count,ind_run,con2_count));
        
    end;
%     if nr_para1 > 0
%         for ind_para = 1:nr_para1
%             for ind_run = 1:runs
%                 fprintf('...load 4D image parametric for run %d and contrast %s...\n',ind_run,con1);
%                 file1 = [results_dir f '4D_' con1 '_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
%                 eval(sprintf('FourD%d = file1;',ind_run));
%                 eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
%                 eval(sprintf('FourD%d_%d_par%d = FourD%d.img;',ind_run,con1_count,ind_para,ind_run));
%                 eval(sprintf('dims = size(FourD%d_%d_par%d(:,:,:,1));',ind_run,con1_count,ind_para));
%                 x = dims(1);
%                 y = dims(2);
%                 z = dims(3);
%                 eval(sprintf('mean_FourD%d_%d = zeros(x,y,z);',ind_run,con1_count));
%                 fprintf('...generate mean activation per voxel for run %d and contrast %s...\n',ind_run,con1); 
%                 for ind_x = 1:x
%                     for ind_y = 1:y
%                         for ind_z = 1:z 
%                             eval(sprintf('mean_FourD%d_%d(ind_x,ind_y,ind_z) = mean(FourD%d_%d_par%d(ind_x,ind_y,ind_z,:));',ind_run,con1_count,ind_run,con1_count,ind_para));
%                         end;
%                     end;
%                 end;
%                 eval(sprintf('mean_FourD%d_%d_par%d = mean_FourD%d_%d(~isnan(mean_FourD%d_%d));',ind_run,con1_count,ind_para,ind_run,con1_count,ind_run,con1_count));
%             end;
%         end;
%     end;
%     if nr_para2 > 0
%         for ind_para = 1:nr_para2
%             for ind_run = 1:runs    
%                 fprintf('...load 4D image parametric for run %d and contrast %s...',ind_run,con2);
%                 file2 = [results_dir f '4D_' con2 '_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
%                 eval(sprintf('FourD%d = file2;',ind_run));
%                 eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
%                 eval(sprintf('FourD%d_%d_par%d = FourD%d.img;',ind_run,con2_count,ind_para,ind_run));
%                 eval(sprintf('dims = size(FourD%d_%d_par%d(:,:,:,1));',ind_run,con2_count,ind_para));
%                 x = dims(1);
%                 y = dims(2);
%                 z = dims(3);
%                 fprintf('...generate mean activation per voxel for run %d and contrast %s...\n',ind_run,con2); 
%                 eval(sprintf('mean_FourD%d_%d = zeros(x,y,z);',ind_run,con2_count));
%                 for ind_x = 1:x
%                     for ind_y = 1:y
%                         for ind_z = 1:z 
%                             eval(sprintf('mean_FourD%d_%d(ind_x,ind_y,ind_z) = mean(FourD%d_%d_par%d(ind_x,ind_y,ind_z,:));',ind_run,con2_count,ind_run,con2_count,ind_para));
%                         end;
%                     end;
%                 end;
%                 eval(sprintf('mean_FourD%d_%d_par%d = mean_FourD%d_%d(~isnan(mean_FourD%d_%d));',ind_run,con2_count,ind_para,ind_run,con2_count,ind_run,con2_count));
% 
%             end;
%         end;
%     end;
elseif split == 1
for ind_run = 1:runs
        fprintf('...load 4D image for run %d and split 1...\n',ind_run);
        file1 = [results_dir f '4D_split1_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file1;',ind_run));
        eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d_split1 = FourD%d.img;',ind_run,ind_run));
        eval(sprintf('dims = size(FourD%d_split1(:,:,:,1));',ind_run));
        x = dims(1);
        y = dims(2);
        z = dims(3);
        eval(sprintf('mean_FourD%d_split1 = zeros(x,y,z);',ind_run));
        fprintf('...generate mean activation per voxel for run %d and split1...\n',ind_run); 
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z 
                    eval(sprintf('mean_FourD%d_split1(ind_x,ind_y,ind_z) = mean(FourD%d_split1(ind_x,ind_y,ind_z,:));',ind_run,ind_run));
                end;
            end;
        end;
        eval(sprintf('mean_FourD%d_split1 = mean_FourD%d_split1(~isnan(mean_FourD%d_split1));',ind_run,ind_run,ind_run));

        fprintf('...load 4D image for run %d and split 2...\n',ind_run);
        file2 = [results_dir f '4D_split2_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file2;',ind_run));
        eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d_split2 = FourD%d.img;',ind_run,ind_run));
        eval(sprintf('dims = size(FourD%d_split2(:,:,:,1));',ind_run));
        x = dims(1);
        y = dims(2);
        z = dims(3);
        eval(sprintf('mean_FourD%d_split2 = zeros(x,y,z);',ind_run));
        fprintf('...generate mean activation per voxel for run %d and split2...\n',ind_run); 
        for ind_x = 1:x
            for ind_y = 1:y
                for ind_z = 1:z 
                    eval(sprintf('mean_FourD%d_split2(ind_x,ind_y,ind_z) = mean(FourD%d_split2(ind_x,ind_y,ind_z,:));',ind_run,ind_run));
                end;
            end;
        end;
        eval(sprintf('mean_FourD%d_split2 = mean_FourD%d_split2(~isnan(mean_FourD%d_split2));',ind_run,ind_run,ind_run));
end;  
if nr_para > 0 
    for ind_para = 1:nr_para
        for ind_run = 1:runs
                fprintf('...load 4D image parametric for run %d and split 1...\n',ind_run);
                file1 = [results_dir f '4D_split1_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];               
                eval(sprintf('FourD%d = file1;',ind_run));
                eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d_split1_par%d = FourD%d.img;',ind_run,ind_para,ind_run));
                eval(sprintf('dims = size(FourD%d_split1_par%d(:,:,:,1));',ind_run,ind_para));
                x = dims(1);
                y = dims(2);
                z = dims(3);
                eval(sprintf('mean_FourD%d_split1 = zeros(x,y,z);',ind_run));
                fprintf('...generate mean activation per voxel for run %d and split1...\n',ind_run); 
                for ind_x = 1:x
                    for ind_y = 1:y
                        for ind_z = 1:z 
                            eval(sprintf('mean_FourD%d_split1(ind_x,ind_y,ind_z) = mean(FourD%d_split1_par%d(ind_x,ind_y,ind_z,:));',ind_run,ind_run,ind_para));
                        end;
                    end;
                end;
                eval(sprintf('mean_FourD%d_split1_par%d = mean_FourD%d_split1(~isnan(mean_FourD%d_split1));',ind_run,ind_para,ind_run,ind_run));

                fprintf('...load 4D image parametric for run %d and split 2...\n',ind_run);
                file2 = [results_dir f '4D_split2_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file2;',ind_run));
                eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d_split2_par%d = FourD%d.img;',ind_run,ind_para,ind_run));
                eval(sprintf('dims = size(FourD%d_split2_par%d(:,:,:,1));',ind_run,ind_para));
                x = dims(1);
                y = dims(2);
                z = dims(3);
                eval(sprintf('mean_FourD%d_split2 = zeros(x,y,z);',ind_run));
                fprintf('...generate mean activation per voxel for run %d and split2...\n',ind_run); 
                for ind_x = 1:x
                    for ind_y = 1:y
                        for ind_z = 1:z 
                            eval(sprintf('mean_FourD%d_split2(ind_x,ind_y,ind_z) = mean(FourD%d_split2_par%d(ind_x,ind_y,ind_z,:));',ind_run,ind_run,ind_para));
                        end;
                    end;
                end;
                eval(sprintf('mean_FourD%d_split2_par%d = mean_FourD%d_split2(~isnan(mean_FourD%d_split2));',ind_run,ind_para,ind_run,ind_run));
        end;
    end;
end;
end;

%% calculation of similarity
if two_cons == 0 && split == 0
for ind_run = 1:runs
    for count = 0:runs-1
        if ind_run+count <= runs
                fprintf('...compare session %d to %d...\n',ind_run,ind_run+count);
            for i = 1:nr_subj
                for j = 1:nr_subj
                    temp_nii_1 = [];
                    temp_nii_2 = [];
                    eval(sprintf('temp_nii_1 = FourD%d(:,:,:,i);',ind_run));
                    eval(sprintf('temp_nii_2 = FourD%d(:,:,:,j);',ind_run+count));
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
               eval(sprintf('mean_temp = mean_FourD%d;',ind_run+count)); 
               [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
               out.r_mat(i,j+1) = out.r(1,2);
               out.p_mat(i,j+1) = out.p(1,2);      
            end;
            
        
        cd(results_dir);
        roi_name=get(handles.name_roi,'String');

        if use_roi == 1
            eval(sprintf('save similarity_%s_%s-%d-%d.mat out',roi_name,strtok(strtok(strtok(contrast_def.contrast,'.'),'.'),'.'),ind_run, ind_run+count));           
        else
            eval(sprintf('save similarity_%s-%d-%d.mat out',strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count));
        end;
        
        % color matrix    
        f1 = figure;
        set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
        subplot(2,2,1:2);
        colormap('jet');
        imagesc(out.r_mat);
        caxis([-1,1]);
        colorbar;
        xlabel(sprintf('subjects session %d (last column: similarity to mean image)',ind_run+count));
        ylabel(sprintf('subjects session %d',ind_run));
        if use_roi == 1
            name1=sprintf('similarity-%s-%s-%d-%d',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
        else
            name1=sprintf('similarity-%s-%d-%d',strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run,ind_run+count);
        end;
        title(name1);
        % histograms
        subplot(2,2,3);
        triu_r_mat = triu(out.r_mat,1);
        vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
        h1 = histogram(diag(out.r_mat));
        hold on;
        h2 = histogram(vec_off_diag_r_mat);
        h1.Normalization = 'probability';
        h1.BinWidth = 0.1;
        h2.Normalization = 'probability';
        h2.BinWidth = 0.1;
        xlabel('similarity');
        ylabel('frequency in percentage');
        if use_roi == 1
            name2=sprintf('histogram-%s-%s-%d-%d',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
        else
            name2=sprintf('histogram-%s-%d-%d',strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
        end;
        title(name2);
                
        % ecdf - densitiy plots
        subplot(2,2,4);
        [f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
        plot(x,f,'LineWidth',2,'Color','blue')
        hold on
        plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on;
        clear x f flo fup
        [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
        plot(x,f,'LineWidth',2,'Color','red')
        hold on
        plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
        hold off;
        if use_roi == 1
            name3=sprintf('cumulative-density-%s-%s-%d-%d',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
        else
            name3=sprintf('cumulative-density-%s-%d-%d',strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
        end;    
        title(name3);
        
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fig,name1,'-dpng','-r0')

        close(f1);
        clearvars out f1 ;
        end;
    end;
end;
if nr_para > 0 
    for ind_para = 1:nr_para
        for ind_run = 1:runs
            for count = 0:runs-1
                if ind_run+count <= runs
                        fprintf('...compare parametric session %d to %d...\n',ind_run,ind_run+count);
                    for i = 1:nr_subj
                        for j = 1:nr_subj
                            temp_nii_1 = [];
                            temp_nii_2 = [];
                            eval(sprintf('temp_nii_1 = FourD%d_par%d(:,:,:,i);',ind_run,ind_para));
                            eval(sprintf('temp_nii_2 = FourD%d_par%d(:,:,:,j);',ind_run+count,ind_para));
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
                       eval(sprintf('mean_temp = mean_FourD%d_par%d;',ind_run+count,ind_para)); 
                       [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
                       out.r_mat(i,j+1) = out.r(1,2);
                       out.p_mat(i,j+1) = out.p(1,2);      
                    end;


                cd(results_dir);
                roi_name=get(handles.name_roi,'String');

                if use_roi == 1
                    eval(sprintf('save similarity-par%d-%s-%s-%d-%d.mat out',ind_para,roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count));           
                else
                    eval(sprintf('save similarity-par%d-%s-%d-%d.mat out',ind_para,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count));
                end;

                % color matrix    
                f1 = figure;
                set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
                subplot(2,2,1:2);
                colormap('jet');
                imagesc(out.r_mat);
                caxis([-1,1]);
                colorbar;
                xlabel(sprintf('subjects session %d (last column: similarity to mean image)',ind_run+count));
                ylabel(sprintf('subjects session %d',ind_run));
                if use_roi == 1
                    name1=sprintf('similarity-par%d-%s-%s-%d-%d',ind_para,strtok(strtok(contrast_def.contrast,'.'),'.'),roi_name,ind_run, ind_run+count);
                else
                    name1=sprintf('similarity-par%d-%s-%d-%d',ind_para,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run,ind_run+count);
                end;
                title(name1);
                % histograms
                subplot(2,2,3);
                triu_r_mat = triu(out.r_mat,1);
                vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
                h1 = histogram(diag(out.r_mat));
                hold on;
                h2 = histogram(vec_off_diag_r_mat);
                h1.Normalization = 'probability';
                h1.BinWidth = 0.1;
                h2.Normalization = 'probability';
                h2.BinWidth = 0.1;
                xlabel('similarity');
                ylabel('frequency in percentage');
                if use_roi == 1
                    name2=sprintf('histogram-%s-%s-%d-%d',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
                else
                    name2=sprintf('histogram-%s%d-%d',strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
                end;
                title(name2);

                % ecdf - densitiy plots
                subplot(2,2,4);
[f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
        plot(x,f,'LineWidth',2,'Color','blue')
        hold on
        plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on;
        clear x f flo fup
        [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
        plot(x,f,'LineWidth',2,'Color','red')
        hold on
        plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
        hold off;
                if use_roi == 1
                    name3=sprintf('cumulative-density-%s-%s-%d-%d',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
                else
                    name3=sprintf('cumulative-density-%s-%d-%d',strtok(strtok(contrast_def.contrast,'.'),'.'),ind_run, ind_run+count);
                end;    
                title(name3);

                fig = gcf;
                fig.PaperPositionMode = 'auto';
                print(fig,name1,'-dpng','-r0')

                close(f1);
                clearvars out f1 ;
                end;
            end;
        end;
    end;
end;
elseif two_cons == 1
    for ind_run = 1:runs
        for count = 0:runs-1
            if ind_run+count <= runs
                %con1
                    fprintf('...compare session %d to %d in contrast %d...\n',ind_run,ind_run+count,con1_count);
                for i = 1:nr_subj
                    for j = 1:nr_subj
                        temp_nii_1 = [];
                        temp_nii_2 = [];
                        eval(sprintf('temp_nii_1 = FourD%d_%d(:,:,:,i);',ind_run,con1_count));
                        eval(sprintf('temp_nii_2 = FourD%d_%d(:,:,:,j);',ind_run+count,con1_count));
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
               eval(sprintf('mean_temp = mean_FourD%d_%d;',ind_run+count,con1_count)); 
               [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
               out.r_mat(i,j+1) = out.r(1,2);
               out.p_mat(i,j+1) = out.p(1,2);  
                end;

            cd(results_dir);
            roi_name=get(handles.name_roi,'String');

            if use_roi == 1
                eval(sprintf('save similarity_%s_%d-%d_con%d.mat out',roi_name,ind_run, ind_run+count,con1_count));           
            else
                eval(sprintf('save similarity_%d-%d_con%d.mat out',ind_run, ind_run+count,con1_count));
            end;

            % color matrix    
        f1 = figure;
        set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
        subplot(2,2,1:2);
        colormap('jet');
        imagesc(out.r_mat);
                caxis([-1,1]);
                colorbar;
                xlabel(sprintf('subject session %d (last column: similarity to mean image)',ind_run+count));
        ylabel(sprintf('subject session %d',ind_run));
            if use_roi == 1
                name1=sprintf('similarity-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con1_count);
            else
                name1=sprintf('similarity-%d-%d-con%d',ind_run, ind_run+count,con1_count);
            end;
        title(name1);

            % histograms
        subplot(2,2,3);
        triu_r_mat = triu(out.r_mat,1);
        vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
        h1 = histogram(diag(out.r_mat));
        hold on;
        h2 = histogram(vec_off_diag_r_mat);
        h1.Normalization = 'probability';
        h1.BinWidth = 0.1;
        h2.Normalization = 'probability';
        h2.BinWidth = 0.1;
        xlabel('similarity');
        ylabel('frequency in percentage');
            if use_roi == 1
                name2=sprintf('histograms-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con1_count);
            else
                name2=sprintf('histograms-%d-%d-con%d',ind_run, ind_run+count,con1_count);
            end;
        title(name2);

            % ecdf - densitiy plots
            subplot(2,2,4);
[f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
        plot(x,f,'LineWidth',2,'Color','blue')
        hold on
        plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on;
        clear x f flo fup
        [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
        plot(x,f,'LineWidth',2,'Color','red')
        hold on
        plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
        hold off;
            if use_roi == 1
                name3=sprintf('density-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con1_count);
            else
                name3=sprintf('density-%d-%d-con%d',ind_run, ind_run+count,con1_count);
            end;     
            title(name3);
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fig,name1,'-dpng','-r0')
        close(f1);
        clearvars out f1;
            %con2
                   fprintf('...compare session %d to %d in contrast %d...\n',ind_run,ind_run+count,con2_count);
                for i = 1:nr_subj
                    for j = 1:nr_subj
                        temp_nii_1 = [];
                        temp_nii_2 = [];
                        eval(sprintf('temp_nii_1 = FourD%d_%d(:,:,:,i);',ind_run,con2_count));
                        eval(sprintf('temp_nii_2 = FourD%d_%d(:,:,:,j);',ind_run+count,con2_count));
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
               eval(sprintf('mean_temp = mean_FourD%d_%d;',ind_run+count,con2_count)); 
               [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
               out.r_mat(i,j+1) = out.r(1,2);
               out.p_mat(i,j+1) = out.p(1,2);                      
                end;

            cd(results_dir);
            roi_name=get(handles.name_roi,'String');

            if use_roi == 1
                eval(sprintf('save similarity_%s_%d-%d-con%d.mat out',roi_name,ind_run, ind_run+count,con2_count));           
            else
                eval(sprintf('save similarity_%d-%d-con%d.mat out',ind_run, ind_run+count,con2_count));
            end;

            % color matrix    
        f1 = figure;
        set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
        subplot(2,2,1:2);
        colormap('jet');
        imagesc(out.r_mat);
                caxis([-1,1]);
                colorbar;
                xlabel(sprintf('subject session %d (last column: similarity to mean image)',ind_run+count));
        ylabel(sprintf('subject session %d',ind_run));

            if use_roi == 1
                name1=sprintf('similarity-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con2_count);
            else
                name1=sprintf('similarity-%d-%d-con%d',ind_run, ind_run+count,con2_count);
            end;
        title(name1);

            % histograms
        subplot(2,2,3);
        triu_r_mat = triu(out.r_mat,1);
        vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
        h1 = histogram(diag(out.r_mat));
        hold on;
        h2 = histogram(vec_off_diag_r_mat);
        h1.Normalization = 'probability';
        h1.BinWidth = 0.1;
        h2.Normalization = 'probability';
        h2.BinWidth = 0.1;
        xlabel('similarity');
        ylabel('frequency in percentage');
            if use_roi == 1
                name2=sprintf('histograms-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con2_count);
            else
                name2=sprintf('histograms-%d-%d-con%d',ind_run, ind_run+count,con2_count);
            end;
        title(name2);

            % ecdf - densitiy plots
            subplot(2,2,4);
[f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
        plot(x,f,'LineWidth',2,'Color','blue')
        hold on
        plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on;
        clear x f flo fup
        [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
        plot(x,f,'LineWidth',2,'Color','red')
        hold on
        plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
        hold off;
            if use_roi == 1
                name3=sprintf('density-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con2_count);
            else
                name3=sprintf('density-%d-%d-con%d',ind_run, ind_run+count,con2_count);
            end;     
            title(name3);
        
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fig,name1,'-dpng','-r0')
        close(fig);  
        clearvars out f1;
            end;
            
        end;
    end;
    % compares two contrasts within one session
    for  i_run = 1:runs
        fprintf('...compare session %d in contrast %d and %d...\n',i_run,con1_count,con2_count);
        for i = 1:nr_subj
            for j = 1:nr_subj
                temp_nii_1 = [];
                temp_nii_2 = [];
                eval(sprintf('temp_nii_1 = FourD%d_%d(:,:,:,i);',i_run,con1_count));
                eval(sprintf('temp_nii_2 = FourD%d_%d(:,:,:,j);',i_run,con2_count));
                if use_roi == 1
                    temp_nii_1(~r_roi_ind) = 0;
                    temp_nii_2(~r_roi_ind) = 0;
                end;
                temp_1 = temp_nii_1(~isnan(temp_nii_1));
                temp_2 = temp_nii_2(~isnan(temp_nii_2));

                [out.r, out.p] = corrcoef([temp_1,temp_2]);
                out.r_mat(i,j) = out.r(1,2);
                out.p_mat(i,j) = out.p(1,2);
                eval(sprintf('mean_temp = mean_FourD%d_%d;',i_run,con1_count)); 
               [out.r, out.p] = corrcoef([temp_2,mean_temp]); 
               out.r_mat(j,nr_subj+2) = out.r(1,2);
               out.p_mat(j,nr_subj+2) = out.p(1,2); 
            end;
               eval(sprintf('mean_temp = mean_FourD%d_%d;',i_run,con2_count)); 
               [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
               out.r_mat(i,j+1) = out.r(1,2);
               out.p_mat(i,j+1) = out.p(1,2);   
              
        end;

            cd(results_dir);
            roi_name=get(handles.name_roi,'String');

            if use_roi == 1
                eval(sprintf('save similarity_%s_%d-con%d-con%d.mat out',roi_name,i_run,con1_count,con2_count));           
            else
                eval(sprintf('save similarity_%d-con%d-con%d.mat out',i_run,con1_count,con2_count));
            end;

            % color matrix    
        f1 = figure;
        set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
        subplot(2,2,1:2);
        colormap('jet');
        imagesc(out.r_mat);
                caxis([-1,1]);
                colorbar;
                xlabel(sprintf('subject session %d (last column: similarity to mean image)',i_run));
        ylabel(sprintf('subject session %d',i_run));

            if use_roi == 1
                name1=sprintf('similarity-%s-%d-con%d-con%d',roi_name,i_run,con1_count,con2_count);
            else
                name1=sprintf('similarity-%d-con%d-con%d',i_run,con1_count,con2_count);
            end;
           title(name1);
        % histograms
        subplot(2,2,3);
        triu_r_mat = triu(out.r_mat,1);
        vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
        h1 = histogram(diag(out.r_mat));
        hold on;
        h2 = histogram(vec_off_diag_r_mat);
        h1.Normalization = 'probability';
        h1.BinWidth = 0.1;
        h2.Normalization = 'probability';
        h2.BinWidth = 0.1;
        xlabel('similarity');
        ylabel('frequency in percentage');
            if use_roi == 1
                name2=sprintf('histograms-%s-%d-con%d-con%d',roi_name,i_run,con1_count,con2_count);
            else
                name2=sprintf('histograms-%d-con%d-con%d',i_run,con1_count,con2_count);
            end;
            title(name2);
                
        % ecdf - densitiy plots
        subplot(2,2,4);
[f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
        plot(x,f,'LineWidth',2,'Color','blue')
        hold on
        plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on;
        clear x f flo fup
        [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
        plot(x,f,'LineWidth',2,'Color','red')
        hold on
        plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
        hold off;
            if use_roi == 1
                name3=sprintf('density-%s-%d-con%d-con%d',roi_name,i_run,con1_count,con2_count);
            else
                name3=sprintf('density-%d-con%d-con%d',i_run,con1_count,con2_count);
            end;            
        title(name3);
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fig,name1,'-dpng','-r0')
        close(fig);

        clearvars out f1;      
    end;
% if nr_para1 > 0   
%     for ind_para = 1:nr_para1
%             for ind_run = 1:runs
%                 for count = 0:runs-1
%                     if ind_run+count <= runs
%                         %con1
%                             fprintf('...compare session %d to %d in contrast %d...\n',ind_run,ind_run+count,con1_count);
%                         for i = 1:nr_subj
%                             for j = 1:nr_subj
%                                 temp_nii_1 = [];
%                                 temp_nii_2 = [];
%                                 eval(sprintf('temp_nii_1 = FourD%d_%d_par%d(:,:,:,i);',ind_run,con1_count,ind_para));
%                                 eval(sprintf('temp_nii_2 = FourD%d_%d_par%d(:,:,:,j);',ind_run+count,con1_count,ind_para));
%                                 if use_roi == 1
%                                     temp_nii_1(~r_roi_ind) = 0;
%                                     temp_nii_2(~r_roi_ind) = 0;
%                                 end;
%                                 temp_1 = temp_nii_1(~isnan(temp_nii_1));
%                                 temp_2 = temp_nii_2(~isnan(temp_nii_2));
% 
%                                 [out.r, out.p] = corrcoef([temp_1,temp_2]);
%                                 out.r_mat(i,j) = out.r(1,2);
%                                 out.p_mat(i,j) = out.p(1,2);   
%                             end;
%                        eval(sprintf('mean_temp = mean_FourD%d_%d_par%d;',ind_run+count,con1_count,ind_para)); 
%                        [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
%                        out.r_mat(i,j+1) = out.r(1,2);
%                        out.p_mat(i,j+1) = out.p(1,2);  
%                         end;
% 
%                     cd(results_dir);
%                     roi_name=get(handles.name_roi,'String');
% 
%                     if use_roi == 1
%                         eval(sprintf('save similarity-par%d_%s_%d-%d_con%d.mat out',ind_para,roi_name,ind_run, ind_run+count,con1_count));           
%                     else
%                         eval(sprintf('save similarity-par%d_%d-%d_con%d.mat out',ind_para,ind_run, ind_run+count,con1_count));
%                     end;
% 
%                     % color matrix    
%                 f1 = figure;
%                 set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
%                 subplot(2,2,1:2);
%                 colormap('jet');
%                 imagesc(out.r_mat);
%                 caxis([-1,1]);
%                 colorbar;
%                 xlabel(sprintf('subject session %d (last column: similarity to mean image)',ind_run+count));
%                 ylabel(sprintf('subject session %d',ind_run));
%                     if use_roi == 1
%                         name1=sprintf('similarity-par%d-%s-%d-%d-con%d',ind_para,roi_name,ind_run, ind_run+count,con1_count);
%                     else
%                         name1=sprintf('similarity-par%d-%d-%d-con%d',ind_para,ind_run, ind_run+count,con1_count);
%                     end;
%                 title(name1);
% 
%                     % histograms
%                 subplot(2,2,3);
%                 triu_r_mat = triu(out.r_mat,1);
%                 vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
%                 h1 = histogram(diag(out.r_mat));
%                 hold on;
%                 h2 = histogram(vec_off_diag_r_mat);
%                 h1.Normalization = 'probability';
%                 h1.BinWidth = 0.1;
%                 h2.Normalization = 'probability';
%                 h2.BinWidth = 0.1;
%                 xlabel('similarity');
%                 ylabel('frequency in percentage');
%                     if use_roi == 1
%                         name2=sprintf('histograms-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con1_count);
%                     else
%                         name2=sprintf('histograms-%d-%d-con%d',ind_run, ind_run+count,con1_count);
%                     end;
%                 title(name2);
% 
%                     % ecdf - densitiy plots
%                 subplot(2,2,4);
% [f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
%         plot(x,f,'LineWidth',2,'Color','blue')
%         hold on
%         plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
%         hold on
%         plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
%         hold on;
%         clear x f flo fup
%         [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
%         plot(x,f,'LineWidth',2,'Color','red')
%         hold on
%         plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
%         hold on
%         plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
%         hold off;
%                     if use_roi == 1
%                         name3=sprintf('density-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con1_count);
%                     else
%                         name3=sprintf('density-%d-%d-con%d',ind_run, ind_run+count,con1_count);
%                     end;     
%                     title(name3);
%                 fig = gcf;
%                 fig.PaperPositionMode = 'auto';
%                 print(fig,name1,'-dpng','-r0')
%                 close(f1);
%                 clearvars out f1;
%                     %con2
%                            fprintf('...compare session %d to %d in contrast %d...\n',ind_run,ind_run+count,con2_count);
%                         for i = 1:nr_subj
%                             for j = 1:nr_subj
%                                 temp_nii_1 = [];
%                                 temp_nii_2 = [];
%                                 eval(sprintf('temp_nii_1 = FourD%d_%d_par%d(:,:,:,i);',ind_run,con2_count,ind_para));
%                                 eval(sprintf('temp_nii_2 = FourD%d_%d_par%d(:,:,:,j);',ind_run+count,con2_count,ind_para));
%                                 if use_roi == 1
%                                     temp_nii_1(~r_roi_ind) = 0;
%                                     temp_nii_2(~r_roi_ind) = 0;
%                                 end;
%                                 temp_1 = temp_nii_1(~isnan(temp_nii_1));
%                                 temp_2 = temp_nii_2(~isnan(temp_nii_2));
% 
%                                 [out.r, out.p] = corrcoef([temp_1,temp_2]);
%                                 out.r_mat(i,j) = out.r(1,2);
%                                 out.p_mat(i,j) = out.p(1,2);   
%                             end;
%                        eval(sprintf('mean_temp = mean_FourD%d_%d_par%d;',ind_run+count,con2_count,ind_para)); 
%                        [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
%                        out.r_mat(i,j+1) = out.r(1,2);
%                        out.p_mat(i,j+1) = out.p(1,2);                      
%                         end;
% 
%                     cd(results_dir);
%                     roi_name=get(handles.name_roi,'String');
% 
%                     if use_roi == 1
%                         eval(sprintf('save similarity-par%d_%s_%d-%d-con%d.mat out',ind_para,roi_name,ind_run, ind_run+count,con2_count));           
%                     else
%                         eval(sprintf('save similarity-par%d_%d-%d-con%d.mat out',ind_para,ind_run, ind_run+count,con2_count));
%                     end;
% 
%                     % color matrix    
%                 f1 = figure;
%                 set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
%                 subplot(2,2,1:2);
%                 colormap('jet');
%                 imagesc(out.r_mat);
%                 caxis([-1,1]);
%                 colorbar;
%                 xlabel(sprintf('subject session %d (last column: similarity to mean image)',ind_run+count));
%                 ylabel(sprintf('subject session %d',ind_run));
% 
%                     if use_roi == 1
%                         name1=sprintf('similarity-par%d-%s-%d-%d-con%d',ind_para,roi_name,ind_run, ind_run+count,con2_count);
%                     else
%                         name1=sprintf('similarity-par%d-%d-%d-con%d',ind_para,ind_run, ind_run+count,con2_count);
%                     end;
%                 title(name1);
% 
%                     % histograms
%                 subplot(2,2,3);
%                 triu_r_mat = triu(out.r_mat,1);
%                 vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
%                 h1 = histogram(diag(out.r_mat));
%                 hold on;
%                 h2 = histogram(vec_off_diag_r_mat);
%                 h1.Normalization = 'probability';
%                 h1.BinWidth = 0.1;
%                 h2.Normalization = 'probability';
%                 h2.BinWidth = 0.1;
%                 xlabel('similarity');
%                 ylabel('frequency in percentage');
%                     if use_roi == 1
%                         name2=sprintf('histograms-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con2_count);
%                     else
%                         name2=sprintf('histograms-%d-%d-con%d',ind_run, ind_run+count,con2_count);
%                     end;
%                 title(name2);
% 
%                     % ecdf - densitiy plots
%                 subplot(2,2,4);
% [f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
%         plot(x,f,'LineWidth',2,'Color','blue')
%         hold on
%         plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
%         hold on
%         plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
%         hold on;
%         clear x f flo fup
%         [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
%         plot(x,f,'LineWidth',2,'Color','red')
%         hold on
%         plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
%         hold on
%         plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
%         hold off;
%                     if use_roi == 1
%                         name3=sprintf('density-%s-%d-%d-con%d',roi_name,ind_run, ind_run+count,con2_count);
%                     else
%                         name3=sprintf('density-%d-%d-con%d',ind_run, ind_run+count,con2_count);
%                     end;     
%                     title(name3);
% 
%                 fig = gcf;
%                 fig.PaperPositionMode = 'auto';
%                 print(fig,name1,'-dpng','-r0')
%                 close(fig);  
%                 clearvars out f1;
%                     end;
% 
%                 end;
%             end;
%             % compares two contrasts within one session
%             for  i_run = 1:runs
%                 fprintf('...compare session %d in contrast %d and %d...\n',i_run,con1_count,con2_count);
%                 for i = 1:nr_subj
%                     for j = 1:nr_subj
%                         temp_nii_1 = [];
%                         temp_nii_2 = [];
%                         eval(sprintf('temp_nii_1 = FourD%d_%d_par%d(:,:,:,i);',i_run,con1_count,ind_para));
%                         eval(sprintf('temp_nii_2 = FourD%d_%d_par%d(:,:,:,j);',i_run,con2_count,ind_para));
%                         if use_roi == 1
%                             temp_nii_1(~r_roi_ind) = 0;
%                             temp_nii_2(~r_roi_ind) = 0;
%                         end;
%                         temp_1 = temp_nii_1(~isnan(temp_nii_1));
%                         temp_2 = temp_nii_2(~isnan(temp_nii_2));
% 
%                         [out.r, out.p] = corrcoef([temp_1,temp_2]);
%                         out.r_mat(i,j) = out.r(1,2);
%                         out.p_mat(i,j) = out.p(1,2);
%                         eval(sprintf('mean_temp = mean_FourD%d_%d_par%d;',i_run,con1_count,ind_para)); 
%                        [out.r, out.p] = corrcoef([temp_2,mean_temp]); 
%                        out.r_mat(j,nr_subj+2) = out.r(1,2);
%                        out.p_mat(j,nr_subj+2) = out.p(1,2); 
%                     end;
%                        eval(sprintf('mean_temp = mean_FourD%d_%d_par%d;',i_run,con2_count,ind_para)); 
%                        [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
%                        out.r_mat(i,j+1) = out.r(1,2);
%                        out.p_mat(i,j+1) = out.p(1,2);   
% 
%                 end;
% 
%                     cd(results_dir);
%                     roi_name=get(handles.name_roi,'String');
% 
%                     if use_roi == 1
%                         eval(sprintf('save similarity-par%d_%s_%d-con%d-con%d.mat out',ind_para,roi_name,i_run,con1_count,con2_count));           
%                     else
%                         eval(sprintf('save similarity-par%d_%d-con%d-con%d.mat out',ind_para,i_run,con1_count,con2_count));
%                     end;
% 
%                     % color matrix    
%                 f1 = figure;
%                 set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
%                 subplot(2,2,1:2);
%                 colormap('jet');
%                 imagesc(out.r_mat);
%                 caxis([-1,1]);
%                 colorbar;
%                 xlabel(sprintf('subject session %d (last column: similarity to mean image)',i_run));
%                 ylabel(sprintf('subject session %d',i_run));
% 
%                     if use_roi == 1
%                         name1=sprintf('similarity-par%d-%s-%d-con%d-con%d',ind_para,roi_name,i_run,con1_count,con2_count);
%                     else
%                         name1=sprintf('similarity-par%d-%d-con%d-con%d',ind_para,i_run,con1_count,con2_count);
%                     end;
%                    title(name1);
%                 % histograms
%                 subplot(2,2,3);
%                 triu_r_mat = triu(out.r_mat,1);
%                 vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
%                 h1 = histogram(diag(out.r_mat));
%                 hold on;
%                 h2 = histogram(vec_off_diag_r_mat);
%                 h1.Normalization = 'probability';
%                 h1.BinWidth = 0.1;
%                 h2.Normalization = 'probability';
%                 h2.BinWidth = 0.1;
%                 xlabel('similarity');
%                 ylabel('frequency in percentage');
%                     if use_roi == 1
%                         name2=sprintf('histograms-%s-%d-con%d-con%d',roi_name,i_run,con1_count,con2_count);
%                     else
%                         name2=sprintf('histograms-%d-con%d-con%d',i_run,con1_count,con2_count);
%                     end;
%                     title(name2);
% 
%                 % ecdf - densitiy plots
%                 subplot(2,2,4);
%  [f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
%         plot(x,f,'LineWidth',2,'Color','blue')
%         hold on
%         plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
%         hold on
%         plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
%         hold on;
%         clear x f flo fup
%         [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
%         plot(x,f,'LineWidth',2,'Color','red')
%         hold on
%         plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
%         hold on
%         plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
%         hold off;
%                     if use_roi == 1
%                         name3=sprintf('density-%s-%d-con%d-con%d',roi_name,i_run,con1_count,con2_count);
%                     else
%                         name3=sprintf('density-%d-con%d-con%d',i_run,con1_count,con2_count);
%                     end;            
%                 title(name3);
%                 fig = gcf;
%                 fig.PaperPositionMode = 'auto';
%                 print(fig,name1,'-dpng','-r0')
%                 close(fig);
% 
%                 clearvars out f1;      
%             end;
%     end;
% end;
elseif split == 1
    for  i_run = 1:runs
        fprintf('...compare splits session %d ...\n',i_run);
        for i = 1:nr_subj
            for j = 1:nr_subj
                temp_nii_1 = [];
                temp_nii_2 = [];
                eval(sprintf('temp_nii_1 = FourD%d_split1(:,:,:,i);',i_run));
                eval(sprintf('temp_nii_2 = FourD%d_split2(:,:,:,j);',i_run));
                if use_roi == 1
                    temp_nii_1(~r_roi_ind) = 0;
                    temp_nii_2(~r_roi_ind) = 0;
                end;
                temp_1 = temp_nii_1(~isnan(temp_nii_1));
                temp_2 = temp_nii_2(~isnan(temp_nii_2));

                [out.r, out.p] = corrcoef([temp_1,temp_2]);
                out.r_mat(i,j) = out.r(1,2);
                out.p_mat(i,j) = out.p(1,2);
                eval(sprintf('mean_temp = mean_FourD%d_split1(:);',i_run)); 
               [out.r, out.p] = corrcoef([temp_2,mean_temp]); 
               out.r_mat(j,nr_subj+2) = out.r(1,2);
               out.p_mat(j,nr_subj+2) = out.p(1,2); 
            end;
               eval(sprintf('mean_temp = mean_FourD%d_split2(:);',i_run)); 
               [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
               out.r_mat(i,j+1) = out.r(1,2);
               out.p_mat(i,j+1) = out.p(1,2);   
              
        end;

            cd(results_dir);
            roi_name=get(handles.name_roi,'String');

            if use_roi == 1
                eval(sprintf('save similarity_%s_%s_%d_split.mat out',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),i_run));           
            else
                eval(sprintf('save similarity_%s_%d_split.mat out',strtok(strtok(contrast_def.contrast,'.'),'.'),i_run));
            end;

            % color matrix    
        f1 = figure;
        set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
        subplot(2,2,1:2);
        colormap('jet');
        imagesc(out.r_mat);
                caxis([-1,1]);
                colorbar;
                xlabel(sprintf('subjects split 1 session %d (last column: similarity to mean image)',i_run));
        ylabel(sprintf('subjects split 2 session %d',i_run));
            if use_roi == 1
                name1=sprintf('similarity-%s-%s-%d-split',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),i_run);
            else
                name1=sprintf('similarity-%s-%d-split',strtok(strtok(contrast_def.contrast,'.'),'.'),i_run);
            end;
           title(name1);
        % histograms
        subplot(2,2,3);
        triu_r_mat = triu(out.r_mat,1);
        vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
        h1 = histogram(diag(out.r_mat));
        hold on;
        h2 = histogram(vec_off_diag_r_mat);
        h1.Normalization = 'probability';
        h1.BinWidth = 0.1;
        h2.Normalization = 'probability';
        h2.BinWidth = 0.1;
        xlabel('similarity');
        ylabel('frequency in percentage');
            if use_roi == 1
                name2=sprintf('histograms-%s-%s-%d-split',roi_name,strtok(strtok(contrast_def.contrast,'.'),'.'),i_run);
            else
                name2=sprintf('histograms-%s-%d-split',strtok(contrast_def.contrast,'.'),i_run);
            end;
        title(name2);
                
        % ecdf - densitiy plots
        subplot(2,2,4);
        [f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
        plot(x,f,'LineWidth',2,'Color','blue')
        hold on
        plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on;
        clear x f flo fup
        [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
        plot(x,f,'LineWidth',2,'Color','red')
        hold on
        plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
        hold off;
            if use_roi == 1
                name3=sprintf('density-%s-%s-%d-split',roi_name,strtok(contrast_def.contrast,'.'),i_run);
            else
                name3=sprintf('density-%s-%d-split',strtok(contrast_def.contrast,'.'),i_run);
            end;            
         title(name3);
         fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fig,name1,'-dpng','-r0')
        close(fig);
        clearvars out f1 fig;    
    end;    
if nr_para > 0
    for ind_para = 1:nr_para
        for  i_run = 1:runs
            fprintf('...compare splits session %d ...\n',i_run);
            for i = 1:nr_subj
                for j = 1:nr_subj
                    temp_nii_1 = [];
                    temp_nii_2 = [];
                    eval(sprintf('temp_nii_1 = FourD%d_split1_par%d(:,:,:,i);',i_run,ind_para));
                    eval(sprintf('temp_nii_2 = FourD%d_split2_par%d(:,:,:,j);',i_run,ind_para));
                    if use_roi == 1
                        temp_nii_1(~r_roi_ind) = 0;
                        temp_nii_2(~r_roi_ind) = 0;
                    end;
                    temp_1 = temp_nii_1(~isnan(temp_nii_1));
                    temp_2 = temp_nii_2(~isnan(temp_nii_2));

                    [out.r, out.p] = corrcoef([temp_1,temp_2]);
                    out.r_mat(i,j) = out.r(1,2);
                    out.p_mat(i,j) = out.p(1,2);
                    eval(sprintf('mean_temp = mean_FourD%d_split1_par%d;',i_run,ind_para)); 
                   [out.r, out.p] = corrcoef([temp_2,mean_temp]); 
                   out.r_mat(j,nr_subj+2) = out.r(1,2);
                   out.p_mat(j,nr_subj+2) = out.p(1,2); 
                end;
                   eval(sprintf('mean_temp = mean_FourD%d_split2_par%d;',i_run,ind_para)); 
                   [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
                   out.r_mat(i,j+1) = out.r(1,2);
                   out.p_mat(i,j+1) = out.p(1,2);   

            end;

                cd(results_dir);
                roi_name=get(handles.name_roi,'String');

                if use_roi == 1
                    eval(sprintf('save similarity-par%d_%s_%s_%d_split.mat out',ind_para,roi_name,strtok(contrast_def.contrast,'.'),i_run));           
                else
                    eval(sprintf('save similarity-par%d_%s_%d_split.mat out',ind_para,strtok(contrast_def.contrast,'.'),i_run));
                end;

                % color matrix    
            f1 = figure;
            set(gcf,'units','normalized','position',[0.25 0 0.50 1]);
            subplot(2,2,1:2);
            colormap('jet');
            imagesc(out.r_mat);
                caxis([-1,1]);
                colorbar;
                xlabel(sprintf('subjects split 1 session %d (last column: similarity to mean image)',i_run));
            ylabel(sprintf('subjects split 2 session %d',i_run));
                if use_roi == 1
                    name1=sprintf('similarity-par%d-%s-%s-%d-split',ind_para,roi_name,strtok(contrast_def.contrast,'.'),i_run);
                else
                    name1=sprintf('similarity-par%d-%s-%d-split',ind_para,strtok(contrast_def.contrast,'.'),i_run);
                end;
               title(name1);
            % histograms
            subplot(2,2,3);
            triu_r_mat = triu(out.r_mat,1);
            vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
            h1 = histogram(diag(out.r_mat));
            hold on;
            h2 = histogram(vec_off_diag_r_mat);
            h1.Normalization = 'probability';
            h1.BinWidth = 0.1;
            h2.Normalization = 'probability';
            h2.BinWidth = 0.1;
            xlabel('similarity');
            ylabel('frequency in percentage');
                if use_roi == 1
                    name2=sprintf('histograms-%s-%s-%d-split',roi_name,strtok(contrast_def.contrast,'.'),i_run);
                else
                    name2=sprintf('histograms-%d-%s-split',i_run,strtok(contrast_def.contrast,'.'));
                end;
            title(name2);

            % ecdf - densitiy plots
            subplot(2,2,4);
[f,x,flo,fup]=ecdf(diag(out.r_mat),'bounds','on');
        plot(x,f,'LineWidth',2,'Color','blue')
        hold on
        plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
        hold on;
        clear x f flo fup
        [f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
        plot(x,f,'LineWidth',2,'Color','red')
        hold on
        plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
        hold on
        plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
        hold off;
                if use_roi == 1
                    name3=sprintf('density-%s-%s-%d-split',roi_name,strtok(contrast_def.contrast,'.'),i_run);
                else
                    name3=sprintf('density-%s-%d-split',strtok(contrast_def.contrast,'.'),i_run);
                end;            
             title(name3);
             fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fig,name1,'-dpng','-r0')
            close(fig);
            clearvars out f1 fig;    
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
