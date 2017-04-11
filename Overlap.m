function varargout = Overlap(varargin)
% OVERLAP MATLAB code for Overlap.fig
%      OVERLAP, by itself, creates a new OVERLAP or raises the existing
%      singleton*.
%
%      H = OVERLAP returns the handle to a new OVERLAP or the handle to
%      the existing singleton*.
%
%      OVERLAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OVERLAP.M with the given input arguments.
%
%      OVERLAP('Property','Value',...) creates a new OVERLAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Overlap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Overlap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Overlap

% Last Modified by GUIDE v2.5 20-Mar-2017 13:59:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Overlap_OpeningFcn, ...
                   'gui_OutputFcn',  @Overlap_OutputFcn, ...
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


% --- Executes just before Overlap is made visible.
function Overlap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Overlap (see VARARGIN)

% Choose default command line output for Overlap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Overlap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Overlap_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes during object creation, after setting all properties.
function roi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi as text
%        str2double(get(hObject,'String')) returns contents of roi as a double

function roi_dir_Callback(hObject, eventdata, handles)
% hObject    handle to roi_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi_dir = spm_select(1,'dir','choose roi directory');
assignin('base','roi_dir',roi_dir);

function p_Callback(hObject, eventdata, handles)
% hObject    handle to p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p as text
%        str2double(get(hObject,'String')) returns contents of p as a double


% --- Executes during object creation, after setting all properties.
function p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reslice.
function split_Callback(hObject, eventdata, handles)
% hObject    handle to reslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reslice


% --- Executes on button press in con_def.
function con_def_Callback(hObject, eventdata, handles)
% hObject    handle to con_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

% --- Executes on button press in group.
function group_Callback(hObject, eventdata, handles)
% hObject    handle to group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of group

% --- Executes on button press in group_path.
function group_path_Callback(hObject, eventdata, handles)
% hObject    handle to group_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
group_path = spm_select(1,'dir','choose group stats directory');
assignin('base','group_path',group_path);

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% set parameters
disp('starting calulcation of overlap');
study_design=evalin('base','study_design');
contrast_def=evalin('base','contrast_def');

roi=get(handles.roi,'String');
roi_dir=evalin('base','roi_dir');

%get GUI input
p = str2double(get(handles.p,'String'));
split = get(handles.split,'value');
group = get(handles.group,'value');

%get study design information
runs=str2double(study_design.number_sessions); 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
stats_dir=study_design.stats_directory;
stats_path=study_design.stats_path;
results_dir = study_design.results_directory;
if split == 1
    split_dir = study_design.split_directory;
    nr_para = study_design.number_parametric;
end;
box_path = pwd; 
if runs == 1
    single_run = str2double(study_design.identifier_session);
end;
cd(results_dir);

%load contrast information
two_cons = contrast_def.two_contrasts;
if two_cons == 0
    con=contrast_def.contrast;
    con_count=contrast_def.number_contrast;
else
    con1=contrast_def.contrast1;
    con2=contrast_def.contrast2;
    con1_count=contrast_def.number_contrast1;
    con2_count=contrast_def.number_contrast2;
end;



% load ROI
cd(roi_dir);
roi_ful = dir(sprintf('%s*',roi));
if isstruct(roi_ful)
    if length(roi_ful)==2
        roi_ful = roi_ful(2).name;
    else
        roi_ful = roi_ful(1).name;
    end;
end;
compl = sprintf('%s\\%s',roi_dir, roi_ful);

%% reslice ROI

disp('...reslicing ROI...');
stats_filled = sprintf(stats_dir,1);
temp = sprintf('%s\\%s\\%s\\%s,1',stats_path,vp{1},stats_filled,con);
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
r_roi=dir(sprintf('r%s*',roi));
if length(r_roi)==2
    compl1 = sprintf('%s\\r%s.img',roi_dir,roi);
    movefile(compl1,results_dir,'f');
    compl2 = sprintf('%s\\r%s.hdr',roi_dir,roi);
    movefile(compl2,results_dir,'f');
    cd(results_dir);    
    r_roi = load_nii(sprintf('r%s.img',roi));
    r_roi_ind = r_roi.img==1;
else
    if ~strcmp(roi_dir,results_dir)
    compl = sprintf('%s\\r%s.nii',roi_dir,roi);
    movefile(compl,results_dir,'f');
    end;
    cd(results_dir);    
    r_roi = load_nii(sprintf('r%s.nii',roi));
    r_roi_ind = r_roi.img==1;        
end;

[roix,roiy,roiz] = ind2sub(size(r_roi_ind),find(r_roi_ind == 1));
roix=roix';
roiy=roiy';
roiz=roiz';
roi_xyz(1,:) = roix;
roi_xyz(2,:) = roiy;
roi_xyz(3,:) = roiz;


%% calculate overlap measures
if group == 0
if split == 0 && two_cons == 0
    for i = 1:length(vp)
        eval(sprintf('results_%d = [];',i));
        % load xSPM and extract suprathresholded voxels
        for j = 1:runs
            co_temp=[];

            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = sprintf('%s\\%s\\%s',stats_path,vp{i},stats_dir_filled);
            cd(box_path);
            fprintf('... create xSPM for %s in session %d...\n',vp{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,con_count);
            eval(sprintf('co_sig%d = xSPM.XYZ;',j));


            % apply ROI
            fprintf('...apply ROI...\n')
            for count_roi = 1:length(roi_xyz)
                fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                for count_sig = 1:length(eval(sprintf('co_sig%d(1,:)',j)))
                    eval(sprintf('temp=co_sig%d(:,count_sig);',j));
                    if isequal(temp,roi_xyz(:,count_roi))
                        co_temp(:,end+1)=temp;
                        fprintf('...significant...\n');
                    end;
                end;
            end;
            eval(sprintf('co_sig%d= co_temp;',j));
            eval(sprintf('nr_sig%d = size(co_sig%d,2);',j,j));
            clear co_temp;
        end;

        % compare suprathresholded voxels
        for k = runs:-1:2
            for ind = 1:k-1
                fprintf('...compare %d to %d...\n',k,k-ind);
                if eval(sprintf('nr_sig%d',k)) > eval(sprintf('nr_sig%d',k-ind))
                    count = 0;
                        for l = 1:eval(sprintf('nr_sig%d',k))
                            for m = 1:eval(sprintf('nr_sig%d',k-ind))
                                comp = isequal(eval(sprintf('co_sig%d(:,l)',k)),eval(sprintf('co_sig%d(:,m)',k-ind)));
                                count = count + comp;
                            end;
                        end;
                else
                    count = 0;
                        for l = 1:eval(sprintf('nr_sig%d',k-ind))
                            for m = 1:eval(sprintf('nr_sig%d',k))
                                comp = isequal(eval(sprintf('co_sig%d(:,m)',k)),eval(sprintf('co_sig%d(:,l)',k-ind)));
                                count = count + comp;
                            end;
                        end;
                end;
                eval(sprintf('overlap_%d_%d = count;',k-ind,k));
                if eval(sprintf('overlap_%d_%d',k-ind,k)) == 0
                    eval(sprintf('dice_%d_%d = 0;',k-ind,k));
                    eval(sprintf('jaccard_%d_%d = 0;',k-ind,k));
                else
                    eval(sprintf('dice_%d_%d = (2.*overlap_%d_%d)./(nr_sig%d+nr_sig%d);',k-ind,k,k-ind,k,k-ind,k));
                    eval(sprintf('jaccard_%d_%d = overlap_%d_%d./(nr_sig%d+nr_sig%d-overlap_%d_%d);',k-ind,k,k-ind,k,k-ind,k,k-ind,k));
                end;


            end;


        end;
           for k1 = runs:-1:2
            for ind1 = 1:k1-1 
            % subject result table
                eval(sprintf('results_%d(1,end+1)=[dice_%d_%d];',i,k1-ind1,k1));
                eval(sprintf('results_%d(1,end+1)=[jaccard_%d_%d];',i,k1-ind1,k1));
            end;
           end;
        eval(sprintf('results(i,:)=results_%d;',i));
        eval(sprintf('clear results_%d',i));

    end;

 %create results table with names
 row = {};   
 row = vp;
    %generate column names
    cols = {};
    for k2 = runs:-1:2
            for ind2 = 1:k2-1 
            % subject result table
                cols{1,end+1}=(sprintf('dice_%d_%d;',k2-ind2,k2));
                cols{1,end+1}=(sprintf('jaccard_%d_%d;',k2-ind2,k2));
            end;
    end;
elseif split == 1
    
    for i = 1:length(vp)
        eval(sprintf('results_%d = [];',i));
        % load xSPM and extract suprathresholded voxels
        for j = 1:runs
            if runs == 1
                j = single_run;
            end;
        %split1
            co_temp=[];
            stats_dir_filled = sprintf(split_dir,j);
            SPM_path = sprintf('%s\\%d\\%s',stats_path,vp{i},stats_dir_filled);
            cd(box_path);
            fprintf('... create xSPM for %d split 1 in session %d...\n',vp(i),j)

            xSPM = create_xSPM(SPM_path,box_path,p,1);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);

            % apply ROI
            fprintf('...apply ROI...\n')
            for count_roi = 1:length(roi_xyz)
                for count_sig = 1:length(eval(sprintf('co_sig%d(1,:)',j)))
                    if isequal(eval(sprintf('co_sig%d(:,count_sig)',j)),roi_xyz(:,count_roi))
                        str = sprintf('co_temp(:,end+1)=co_sig%d(:,count_sig);',j);
                        eval(str);
                    end;
                end;
            end;
            evalstr2 = sprintf('co_sig_split1_%d= co_temp;',j);
            eval(evalstr2);  
            evalstr = sprintf('nr_sig_split1_%d = size(co_sig_split1_%d,2);',j,j);
            eval(evalstr);        
            clear co_temp;
      %split2
            co_temp=[];
            stats_dir_filled = sprintf(split_dir,j);
            SPM_path = sprintf('%s\\%d\\%s',stats_path,vp{i},stats_dir_filled);
            cd(box_path);
            fprintf('... create xSPM for %d split 1 in session %d...\n',vp(i),j)

            xSPM = create_xSPM(SPM_path,box_path,p,2+nr_para);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);

            % apply ROI
            fprintf('...apply ROI...\n')
            for count_roi = 1:length(roi_xyz)
                for count_sig = 1:length(eval(sprintf('co_sig%d(1,:)',j)))
                    if isequal(eval(sprintf('co_sig%d(:,count_sig)',j)),roi_xyz(:,count_roi))
                        str = sprintf('co_temp(:,end+1)=co_sig%d(:,count_sig);',j);
                        eval(str);
                    end;
                end;
            end;
            evalstr2 = sprintf('co_sig_split2_%d= co_temp;',j);
            eval(evalstr2);  
            evalstr = sprintf('nr_sig_split2_%d = size(co_sig_split2_%d,2);',j,j);
            eval(evalstr);        
            clear co_temp;
        end;

        % compare suprathresholded voxels
        for k = 1:runs
            if runs == 1
                k = single_run;
            end;
                fprintf('...compare split1 to split2 in session %d...\n',k);
                if eval(sprintf('nr_sig_split1_%d',k)) > eval(sprintf('nr_sig_split2_%d',k))
                    count = 0;
                        for l = 1:eval(sprintf('nr_sig_split1_%d',k))
                            for m = 1:eval(sprintf('nr_sig_split2_%d',k))
                                comp = isequal(eval(sprintf('co_sig_split1_%d(:,l)',k)),eval(sprintf('co_sig_split2_%d(:,m)',k)));
                                count = count + comp;
                            end;
                        end;
                else
                    count = 0;
                        for l = 1:eval(sprintf('nr_sig_split2_%d',k))
                            for m = 1:eval(sprintf('nr_sig_split1_%d',k))
                                comp = isequal(eval(sprintf('co_sig_split1_%d(:,m)',k)),eval(sprintf('co_sig_split2_%d(:,l)',k)));
                                count = count + comp;
                            end;
                        end;
                end;
                eval(sprintf('overlap_%d = count;',k));
                if eval(sprintf('overlap_%d',k)) == 0
                    eval(sprintf('dice_%d = 0;',k));
                    eval(sprintf('jaccard_%d = 0;',k));
                else
                    eval(sprintf('dice_%d = (2.*overlap_%d)./(nr_sig_split1_%d+nr_sig_split2_%d);',k,k,k,k));
                    eval(sprintf('jaccard_%d = overlap_%d./(nr_sig_split1_%d+nr_sig_split2_%d-overlap_%d);',k,k,k,k,k));
                end;
                % subject result table
                eval(sprintf('results_%d(1,end+1)=[dice_%d];',i,k));
                eval(sprintf('results_%d(1,end+1)=[jaccard_%d];',i,k));
            end;

           end;
        eval(sprintf('results(i,:)=results_%d;',i));
        eval(sprintf('clear results_%d',i)); 

 %create results table with names
 row = {};   
 row = vp;
    %generate column names
    cols = {};
    for k = 1:runs
        if runs == 1
            k = single_run;
        end;
           cols{1,end+1}=(sprintf('dice_%d;',k));
           cols{1,end+1}=(sprintf('jaccard_%d;',k));
    end;

elseif two_cons == 1
    for i = 1:length(vp)
        eval(sprintf('results_%d = [];',i));
        % load xSPM and extract suprathresholded voxels
        for j = 1:runs
            if runs == 1
                j = single_run;
            end;
        %con1
            co_temp=[];
            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = sprintf('%s\\%s\\%s',stats_path,vp{i},stats_dir_filled);
            cd(box_path);
            fprintf('... create xSPM for %s con 1 in session %d...\n',vp{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,con1_count);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);

            % apply ROI
            fprintf('...apply ROI...\n')
            for count_roi = 1:length(roi_xyz)
                for count_sig = 1:length(eval(sprintf('co_sig%d(1,:)',j)))
                    if isequal(eval(sprintf('co_sig%d(:,count_sig)',j)),roi_xyz(:,count_roi))
                        str = sprintf('co_temp(:,end+1)=co_sig%d(:,count_sig);',j);
                        eval(str);
                    end;
                end;
            end;
            evalstr2 = sprintf('co_sig_con1_%d= co_temp;',j);
            eval(evalstr2);  
            evalstr = sprintf('nr_sig_con1_%d = size(co_sig_con1_%d,2);',j,j);
            eval(evalstr);        
            clear co_temp;
      %con2
            co_temp=[];
            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = sprintf('%s\\%s\\%s',stats_path,vp{i},stats_dir_filled);
            cd(box_path);
            fprintf('... create xSPM for %s con 2 in session %d...\n',vp{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,con2_count);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);

            % apply ROI
            fprintf('...apply ROI...\n')
            for count_roi = 1:length(roi_xyz)
                for count_sig = 1:length(eval(sprintf('co_sig%d(1,:)',j)))
                    if isequal(eval(sprintf('co_sig%d(:,count_sig)',j)),roi_xyz(:,count_roi))
                        str = sprintf('co_temp(:,end+1)=co_sig%d(:,count_sig);',j);
                        eval(str);
                    end;
                end;
            end;
            evalstr2 = sprintf('co_sig_con2_%d= co_temp;',j);
            eval(evalstr2);  
            evalstr = sprintf('nr_sig_con2_%d = size(co_sig_con2_%d,2);',j,j);
            eval(evalstr);        
            clear co_temp;
        end;

        % compare suprathresholded voxels
        for k = 1:runs
            if runs == 1
                k = single_run;
            end;
                fprintf('...compare con1 to con2 in session %d...\n',k);
                if eval(sprintf('nr_sig_con1_%d',k)) > eval(sprintf('nr_sig_con2_%d',k))
                    count = 0;
                        for l = 1:eval(sprintf('nr_sig_con1_%d',k))
                            for m = 1:eval(sprintf('nr_sig_con2_%d',k))
                                comp = isequal(eval(sprintf('co_sig_con1_%d(:,l)',k)),eval(sprintf('co_sig_con2_%d(:,m)',k)));
                                count = count + comp;
                            end;
                        end;
                else
                    count = 0;
                        for l = 1:eval(sprintf('nr_sig_con2_%d',k))
                            for m = 1:eval(sprintf('nr_sig_con1_%d',k))
                                comp = isequal(eval(sprintf('co_sig_con1_%d(:,m)',k)),eval(sprintf('co_sig_con2_%d(:,l)',k)));
                                count = count + comp;
                            end;
                        end;
                end;
                eval(sprintf('overlap_%d = count;',k));
                if eval(sprintf('overlap_%d',k)) == 0
                    eval(sprintf('dice_%d = 0;',k));
                    eval(sprintf('jaccard_%d = 0;',k));
                else
                    eval(sprintf('dice_%d = (2.*overlap_%d)./(nr_sig_con1_%d+nr_sig_con2_%d);',k,k,k,k));
                    eval(sprintf('jaccard_%d = overlap_%d./(nr_sig_con1_%d+nr_sig_con2_%d-overlap_%d);',k,k,k,k,k));
                end;
                % subject result table
                eval(sprintf('results_%d(1,end+1)=[dice_%d];',i,k));
                eval(sprintf('results_%d(1,end+1)=[jaccard_%d];',i,k));
        end;
        eval(sprintf('results(i,:)=results_%d;',i));
        eval(sprintf('clear results_%d',i)); 
           end;


 %create results table with names
    %generate column names
     row = {};   
     row = vp;
    cols = {};
    for k = 1:runs
        if runs == 1
            k = single_run;
        end;
           cols{1,end+1}=(sprintf('dice_%d;',k));
           cols{1,end+1}=(sprintf('jaccard_%d;',k));
    end;
    
end;
assignin('base','results',results);
assignin('base','cols',cols);
assignin('base','row',row);


results_overlap_subj=dataset({results,cols{:}}, ...
                'obsnames', row);

            
cd(results_dir);             
save results_overlap_subjectwise.mat results_overlap_subj         
disp('DONE');
cd(box_path);
else
    group_path = evalin('base','group_path');
    % group level overlap measures
    con = spm_input;
    for j = 1:runs
        co_temp=[];
        SPM_path = group_path;
        disp('...please type in last contrast before sessions, e.g. 1 if forst contrast of interest is 2...');
        conj = con + j;
        cd(box_path)
        xSPM = create_xSPM(SPM_path,box_path,p,conj);
        eval(sprintf('co_sig%d = xSPM.XYZ;',j));
        % apply ROI
        fprintf('...apply ROI...\n')
        for count_roi = 1:length(roi_xyz)
            fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
            for count_sig = 1:length(eval(sprintf('co_sig%d(1,:)',j)))
                if isequal(eval(sprintf('co_sig%d(:,count_sig)',j)),roi_xyz(:,count_roi))
                    eval(sprintf('co_temp(:,end+1)=co_sig%d(:,count_sig);',j));
                    fprintf('...significant...\n');
                end;
            end;
        end;
        eval(sprintf('co_sig%d= co_temp;',j));
        eval(sprintf('nr_sig%d = size(co_sig%d,2);',j,j));
        clear co_temp;        
    end;
    
    % compare suprathresholded voxels
    for k = runs:-1:2
        for ind = 1:k-1
            fprintf('...compare %d to %d...\n',k,k-ind);
            if eval(sprintf('nr_sig%d',k)) > eval(sprintf('nr_sig%d',k-ind))
                count = 0;
                    for l = 1:eval(sprintf('nr_sig%d',k))
                        for m = 1:eval(sprintf('nr_sig%d',k-ind))
                            comp = isequal(eval(sprintf('co_sig%d(:,l)',k)),eval(sprintf('co_sig%d(:,m)',k-ind)));
                            count = count + comp;
                        end;
                    end;
            else
                count = 0;
                    for l = 1:eval(sprintf('nr_sig%d',k-ind))
                        for m = 1:eval(sprintf('nr_sig%d',k))
                            comp = isequal(eval(sprintf('co_sig%d(:,m)',k)),eval(sprintf('co_sig%d(:,l)',k-ind)));
                            count = count + comp;
                        end;
                    end;
            end;
            eval(sprintf('overlap_%d_%d = count;',k-ind,k));
            if eval(sprintf('overlap_%d_%d',k-ind,k)) == 0
                eval(sprintf('dice_%d_%d = 0;',k-ind,k));
                eval(sprintf('jaccard_%d_%d = 0;',k-ind,k));
            else
                eval(sprintf('dice_%d_%d = (2.*overlap_%d_%d)./(nr_sig%d+nr_sig%d);',k-ind,k,k-ind,k,k-ind,k));
                eval(sprintf('jaccard_%d_%d = overlap_%d_%d./(nr_sig%d+nr_sig%d-overlap_%d_%d);',k-ind,k,k-ind,k,k-ind,k,k-ind,k));
            end;


        end;


    end;
    
    %create results table
    results_overlap = [];
    for k1 = runs:-1:2
        for ind1 = 1:k1-1 
            eval(sprintf('results_overlap(1,end+1)=[dice_%d_%d];',k1-ind1,k1));
            eval(sprintf('results_overlap(1,end+1)=[jaccard_%d_%d];',k1-ind1,k1));
        end;
    end;
    %generate column names
    cols = {};
    for k2 = runs:-1:2
            for ind2 = 1:k2-1 
            % subject result table
                cols{1,end+1}=(sprintf('dice_%d_%d;',k2-ind2,k2));
                cols{1,end+1}=(sprintf('jaccard_%d_%d;',k2-ind2,k2));
            end;
    end;   
    results_overlap_group=dataset({results_overlap,cols{:}});
    cd(results_dir);             
    save results_overlap_grouplevel.mat results_overlap_group  
    disp('DONE');
end;
