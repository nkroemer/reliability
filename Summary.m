function varargout = Summary(varargin)
% SUMMARY MATLAB code for Summary.fig
%      SUMMARY, by itself, creates a new SUMMARY or raises the existing
%      singleton*.
%
%      H = SUMMARY returns the handle to a new SUMMARY or the handle to
%      the existing singleton*.
%
%      SUMMARY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUMMARY.M with the given input arguments.
%
%      SUMMARY('Property','Value',...) creates a new SUMMARY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Summary_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Summary_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Summary

% Last Modified by GUIDE v2.5 26-Jul-2017 09:54:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Summary_OpeningFcn, ...
                   'gui_OutputFcn',  @Summary_OutputFcn, ...
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


% --- Executes just before Summary is made visible.
function Summary_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Summary (see VARARGIN)

% Choose default command line output for Summary
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Summary wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Summary_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

% --- Executes on button press in icc.
function icc_Callback(hObject, eventdata, handles)
% hObject    handle to icc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of icc


% --- Executes on button press in corr.
function corr_Callback(hObject, eventdata, handles)
% hObject    handle to corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of corr
% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Starting reliability summary...');
%% define file seperator 
f = filesep;
box_path=evalin('base','box_path');

study_design = evalin('base','study_design');
contrast_def = evalin('base','contrast_def');

%% get study design info
results_dir = study_design.results_directory;
load(study_design.subject_list); % loads id list
exStats = study_design.exist_stats;
ex4D = study_design.exist_4D;
if exStats == 1
    stats_path = study_design.stats_path;
    stats_filled = sprintf(study_design.stats_directory,1);
    two_cons = contrast_def.two_contrasts;
    if two_cons == 1
        con = [contrast_def.contrast1 contrast_def.contrast_format];

    else
        con = contrast_def.contrast;
    end;
end;


nr_runs = study_design.number_sessions;


   

%%  reslice atlas and move to results directory
disp('...reslicing atlas...');
cd(box_path);
cd('atlas');
atlas_dir = pwd;
atlas_name = 'atlas.nii';
atlas_compl = [atlas_dir f atlas_name];
if exStats == 1
    matlabbatch{1}.spm.spatial.coreg.write.ref = {[stats_path f id{1} f stats_filled f con ',1']};
else 
    matlabbatch{1}.spm.spatial.coreg.write.ref = {[results_dir f 'template_3D.nii' ',1']};
end;       
    matlabbatch{1}.spm.spatial.coreg.write.source = {[atlas_compl ',1']};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

    % save batch
    save('batch_reslice_atlas','matlabbatch');

    % run batch
    spm_jobman('serial',matlabbatch);

% load resliced atlas
atlas_r = [atlas_dir f f 'r' atlas_name];
movefile (atlas_r,results_dir,'f');
cd(results_dir);
atlas = load_nii(sprintf('r%s',atlas_name));
atlas = atlas.img;

cd(atlas_dir);
load('labels.mat');
load('labels_yeo.mat');

rel_summary={};

rel_summary(2:length(Labels)+1,1) = Labels;
rel_summary(1,1) = {'atlas regions'};
rel_summary(2:length(Labels_yeo)+1,2) = Labels_yeo(:,3);
rel_summary(1,2) = {'yeo networks'};

%% get GUI input
corr = get(handles.corr,'Value');
if corr == 1
    cd(results_dir);
    %corrs = dir('CorrMap*');
    corrsZ = dir('z_Corr*');
    %list_corrs = cell(1,length(corrs));
%     for ind_corrs = 1:length(corrs)
%         tmp = load_nii(corrs(ind_corrs).name);
%         [pathstr,name,ext] = fileparts(corrs(ind_corrs).name) ;
%         eval(sprintf('%s = tmp.img;',name));
%         list_corrs{1,ind_corrs} = name;
%     end;
    %rel_summary(1,3:length(list_corrs)+2) = list_corrs;
    list_corrsZ = cell(1,length(corrsZ));
    for ind_corrsZ = 1:length(corrsZ)
        tmp = load_nii(corrsZ(ind_corrsZ).name);
        [pathstr,name,ext] = fileparts(corrsZ(ind_corrsZ).name) ;
        [a,b]=strtok(name,'z_');
        name=[a b];
        eval(sprintf('%s = tmp.img;',name));
        list_corrsZ{1,ind_corrsZ} = name;
    end;
    rel_summary(1,end+1:length(rel_summary(1,:))+length(list_corrsZ)) = list_corrsZ;
end;
icc = get(handles.icc,'Value');
if icc == 1
    cd(results_dir);
    %iccs = dir('ICC*');
    iccsZ = dir('z_ICC*');
    %list_iccs = cell(1,length(iccs));
%     for ind_iccs = 1:length(iccs)
%         tmp = load_nii(iccs(ind_iccs).name);
%         [pathstr,name,ext] = fileparts(iccs(ind_iccs).name) ;
%         eval(sprintf('%s = tmp.img;',name));
%         list_iccs{1,ind_iccs} = name;
%     end;
%     rel_summary(1,end+1:length(rel_summary(1,:))+length(list_iccs)) = list_iccs;
    list_iccsZ = cell(1,length(iccsZ));
    for ind_iccsZ = 1:length(iccsZ)
        tmp = load_nii(iccsZ(ind_iccsZ).name);
        [pathstr,name,ext] = fileparts(iccsZ(ind_iccsZ).name) ;
        [a,b]=strtok(name,'z_');
        name=[a b];
        eval(sprintf('%s = tmp.img;',name));
        list_iccsZ{1,ind_iccsZ} = name;
    end;
    rel_summary(1,end+1:length(rel_summary(1,:))+length(list_iccsZ)) = list_iccsZ;
end;
%% calculate reliability for ROIs

for ind_labels = 1:length(Labels)
    fprintf('...analyzing voxels in %s...\n',Labels{ind_labels})
    fprintf('...corresponding to yeo network %d...\n',Labels_yeo{ind_labels,2})
    [x,y,z] = ind2sub([size(atlas,1),size(atlas,2),size(atlas,3)],find(atlas==ind_labels));
    if icc == 1
        icc_vox = zeros(length(find(atlas==ind_labels)),length(list_iccsZ));
    end;
    if corr == 1
        corr_vox = zeros(length(find(atlas==ind_labels)),length(list_corrsZ));
    end;        
    for ind_vox = 1:length(find(atlas==ind_labels))
        if icc == 1
%             for ind_list = 1:length(list_iccs)
%                eval(sprintf('tmp =  %s;',list_iccs{1,ind_list}));
%                icc_vox(ind_vox,ind_list) = mean(tmp(x(ind_vox),y(ind_vox),z(ind_vox),:),'omitnan');
%             end;
            for ind_list = 1:length(list_iccsZ)
               eval(sprintf('icc_vox(ind_vox,ind_list) = mean(%s(x(ind_vox),y(ind_vox),z(ind_vox),:));',list_iccsZ{1,ind_list}));
            end;           
        end;
        if corr == 1
%             for ind_list = 1:length(list_corrs)
%                eval(sprintf('corr_vox(ind_vox,ind_list) = mean(%s(x(ind_vox),y(ind_vox),z(ind_vox),:));',list_corrs{1,ind_list}));
%             end;
            for ind_list = 1:length(list_corrsZ)
               eval(sprintf('corr_vox(ind_vox,ind_list) = mean(%s(x(ind_vox),y(ind_vox),z(ind_vox),:));',list_corrsZ{1,ind_list}));
            end;            
        end;
    end
    
    if corr == 1 && icc == 1
        for ind_list = 1:length(corr_vox(1,:))
            temp=mean(corr_vox(:,ind_list),'omitnan');
            rel_summary{ind_labels+1,ind_list+2} = (exp(2*temp)-1)./(exp(2*temp)+1);
        end;
        for ind_list = 1:length(icc_vox(1,:))
            temp = mean(icc_vox(:,ind_list),'omitnan');
            rel_summary{ind_labels+1,length(corr_vox(1,:))+ind_list+2} = (exp(2*temp)-1)./(exp(2*temp)+1);
        end;
    elseif corr == 1
        for ind_list = 1:length(corr_vox(1,:))
             temp = mean(corr_vox(:,ind_list),'omitnan');
             rel_summary{ind_labels+1,ind_list+2} =(exp(2*temp)-1)./(exp(2*temp)+1);
        end;
    elseif icc == 1
        for ind_list = 1:length(icc_vox(1,:))
            temp = mean(icc_vox(:,ind_list),'omitnan');
            rel_summary{ind_labels+1,ind_list+2} = (exp(2*temp)-1)./(exp(2*temp)+1);
        end;
    end;

end;

cd(results_dir);
save reliability_summary.mat rel_summary;
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
