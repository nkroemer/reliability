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

% Last Modified by GUIDE v2.5 25-Jul-2017 16:44:47

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
set(handles.roi,'TooltipString','Type in name of ROI file WITHOUT suffix');


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

% --- Executes on button press in use_roi.
function use_roi_Callback(hObject, eventdata, handles)
% hObject    handle to use_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_roi

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% define file seperator 
f = filesep;
box_path=evalin('base','box_path');

%% set parameters
disp('Starting calulcation of overlap...');
study_design=evalin('base','study_design');
contrast_def=evalin('base','contrast_def');
%get GUI input
p = str2double(get(handles.p,'String'));
split = get(handles.split,'value');
group = get(handles.group,'value');
use_roi = get(handles.use_roi,'value');

if use_roi == 1
    roi=get(handles.roi,'String');
    roi_dir=evalin('base','roi_dir');
    str = 'ROI';
else
    str = '';
end;

%% get study design information
runs=study_design.number_sessions; 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
stats_dir=study_design.stats_directory;
stats_path=study_design.stats_path;
results_dir = study_design.results_directory;
if split == 1
    split_dir = study_design.split_directory;
end;
if runs == 1
    single_run = str2double(study_design.identifier_session);
end;
cd(results_dir);

%% load contrast information
two_cons = contrast_def.two_contrasts;
if two_cons == 0
    con=contrast_def.contrast;
    %con_count=contrast_def.number_regressor;
    con_count = contrast_def.contrast_number;
    nr_para = study_design.number_parametric;

else
    con=[contrast_def.contrast1 contrast_def.contrast_format];
    con1=contrast_def.contrast1;
    con2=contrast_def.contrast2;
    con1_count=contrast_def.number_regressor1;
    con2_count=contrast_def.number_regressor2;
   
    nr_para1 = study_design.number_parametric1;
    nr_para2 = study_design.number_parametric2;

end;

%% load ROI
if use_roi == 1
    cd(roi_dir);
    roi_ful = dir(sprintf('%s*',roi));
    if isstruct(roi_ful)
        if length(roi_ful)==2
            roi_ful = roi_ful(2).name;
        else
            roi_ful = roi_ful(1).name;
        end;
    end;
compl = [roi_dir f roi_ful];


%% reslice ROI

disp('...reslicing ROI...');
stats_filled = sprintf(stats_dir,1);
temp = [stats_path f id{1} f stats_filled f con ',1'];
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
    if ~strcmp(roi_dir,results_dir)
        compl1 = [roi_dir f 'r' roi '.img'];
        movefile(compl1,results_dir,'f');
        compl2 = [roi_dir f 'r' roi '.hdr'];
        movefile(compl2,results_dir,'f');
    end;
    cd(results_dir);    
    r_roi = load_untouch_nii(sprintf('r%s.img',roi));
    r_roi_ind = r_roi.img==1;
else
    if ~strcmp(roi_dir,results_dir)
    compl = [roi_dir f 'r' roi '.nii'];
    movefile(compl,results_dir,'f');
    end;
    cd(results_dir);    
    r_roi = load_untouch_nii(sprintf('r%s.nii',roi));
    r_roi_ind = r_roi.img==1;        
end;

[roix,roiy,roiz] = ind2sub(size(r_roi_ind),find(r_roi_ind == 1));
roix=roix';
roiy=roiy';
roiz=roiz';
roi_xyz(1,:) = roix;
roi_xyz(2,:) = roiy;
roi_xyz(3,:) = roiz;
end;

%% calculate overlap measures
if group == 0
if split == 0 && two_cons == 0
    for i = 1:length(id)
        eval(sprintf('results_%d = [];',i));
        % load xSPM and extract suprathresholded voxels
        for j = 1:runs
            co_temp=[];

            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = [stats_path f id{i} f stats_dir_filled];
            cd(box_path);
            fprintf('... create xSPM for %s in session %d...\n',id{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,con_count);
            eval(sprintf('co_sig%d = xSPM.XYZ;',j));

            if use_roi == 1
                % apply ROI
                for count_roi = 1:length(roi_xyz)
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                    if idx>0
                        eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                    end;
                end;
                eval(sprintf('co_sig_%d= co_temp;',j));
            else
                eval(sprintf('co_sig_%d= co_sig%d;',j,j));
            end;
            eval(sprintf('nr_sig_%d = size(co_sig_%d,2);',j,j));
            clear co_temp;           
        end;

        % compare suprathresholded voxels
        for k = runs:-1:2
            for ind = 1:k-1
                fprintf('...compare %d to %d...\n',k,k-ind);
                count = 0;
                for l = 1:eval(sprintf('nr_sig_%d',k))
                   if ~isempty(eval(sprintf('co_sig_%d',k-ind)))
                        eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_%d,co_sig_%d(:,l))));',k-ind,k));
                        count = count + length(idx);
                   end;
                end;
                
                eval(sprintf('overlap_%d_%d = count;',k-ind,k));
                if eval(sprintf('overlap_%d_%d',k-ind,k)) == 0
                    eval(sprintf('dice_%d_%d = 0;',k-ind,k));
                    eval(sprintf('jaccard_%d_%d = 0;',k-ind,k));
                else
                    eval(sprintf('dice_%d_%d = (2.*overlap_%d_%d)./(nr_sig_%d+nr_sig_%d);',k-ind,k,k-ind,k,k-ind,k));
                    eval(sprintf('jaccard_%d_%d = overlap_%d_%d./(nr_sig_%d+nr_sig_%d-overlap_%d_%d);',k-ind,k,k-ind,k,k-ind,k,k-ind,k));
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
 row = id;
 
 %generate column names
 cols = {};
 for k2 = runs:-1:2
    for ind2 = 1:k2-1 
        % subject result table
        cols{1,end+1}=(sprintf('dice_%d_%d;',k2-ind2,k2));
        cols{1,end+1}=(sprintf('jaccard_%d_%d;',k2-ind2,k2));
    end;
 end;
assignin('base','results',results);
assignin('base','cols',cols);
assignin('base','row',row);

cd(results_dir);             

results_overlap_subj=dataset({results,cols{:}}, ...
                'obsnames', row);
f2 = figure;            
imagesc(results,[0,1]);
colormap('inferno');
colorbar;
xlabel(sprintf('%s',cols{:}));
ylabel('subject ID');
title('results overlap subjectwise');
print(f2,sprintf('results_overlap%s_subjectwise',str),'-dpng','-r1200');
            
eval(sprintf('save results%s_overlap_subjectwise_%g.mat results_overlap_subj',str,p)); 
clear results;

if nr_para>0
    for ind_para = 1:nr_para
         for i = 1:length(id)
                eval(sprintf('results_%d = [];',i));
                % load xSPM and extract suprathresholded voxels
                for j = 1:runs
                    co_temp=[];

                    stats_dir_filled = sprintf(stats_dir,j);
                    SPM_path = [stats_path f id{i} f stats_dir_filled];
                    cd(box_path);
                    fprintf('... create xSPM for %s in session %d...\n',id{i},j)

                    xSPM = create_xSPM(SPM_path,box_path,p,con_count+ind_para);
                    eval(sprintf('co_sig%d_par = xSPM.XYZ;',j));

                    if use_roi == 1
                        % apply ROI
                        for count_roi = 1:length(roi_xyz)
                            %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                            eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d_par,roi_xyz(:,count_roi))));',j));
                            if idx>0
                                eval(sprintf('co_temp(:,end+1)=co_sig%d_par(:,idx);',j));
                               % fprintf('...significant...\n');
                            end;
                        end;
                        eval(sprintf('co_sig_%d_par= co_temp;',j));
                    else
                        eval(sprintf('co_sig_%d_par= co_sig%d_par;',j,j));
                    end;
                    eval(sprintf('nr_sig_%d_par = size(co_sig_%d_par,2);',j,j));
                    clear co_temp;           
                 end;

            % compare suprathresholded voxels
            for k = runs:-1:2
                for ind = 1:k-1
                    fprintf('...compare %d to %d...\n',k,k-ind);
                    count = 0;
                    for l = 1:eval(sprintf('nr_sig_%d_par',k))
                       if ~isempty(eval(sprintf('co_sig_%d_par',k-ind)))
                           eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_%d_par,co_sig_%d_par(:,l))));',k-ind,k));
                           count = count + length(idx);
                       end;
                    end;
                    eval(sprintf('overlap_%d_%d = count;',k-ind,k));
                    if eval(sprintf('overlap_%d_%d',k-ind,k)) == 0
                        eval(sprintf('dice_%d_%d = 0;',k-ind,k));
                        eval(sprintf('jaccard_%d_%d = 0;',k-ind,k));
                    else
                        eval(sprintf('dice_%d_%d = (2.*overlap_%d_%d)./(nr_sig_%d_par+nr_sig_%d_par);',k-ind,k,k-ind,k,k-ind,k));
                        eval(sprintf('jaccard_%d_%d = overlap_%d_%d./(nr_sig_%d_par+nr_sig_%d_par-overlap_%d_%d);',k-ind,k,k-ind,k,k-ind,k,k-ind,k));
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
        row = id;
        %generate column names
        cols = {};
        for k2 = runs:-1:2
            for ind2 = 1:k2-1 
                % subject result table
                cols{1,end+1}=(sprintf('dice_%d_%d;',k2-ind2,k2));
                cols{1,end+1}=(sprintf('jaccard_%d_%d;',k2-ind2,k2));
            end;
        end;
        assignin('base','results',results);
        assignin('base','cols',cols);
        assignin('base','row',row);

        cd(results_dir);             

        results_overlap_subj_par=dataset({results,cols{:}}, ...
                        'obsnames', row);
        f2 = figure;            
        imagesc(results,[0,1]);
        colormap('inferno');
        colorbar;
        xlabel(sprintf('%s',cols{:}));
        ylabel('subject ID');
        title('results overlap subjectwise parametric');
        print(f2,sprintf('results_overlap%s_subjectwise_par%d',str,ind_para),'-dpng','-r1200');

        eval(sprintf('save results_overlap%s_subjectwise_par%d_%g.mat results_overlap_subj_par',str,ind_para,p));   
        clear results
    end;
end;
    
disp('DONE');
cd(box_path);    

elseif split == 1
    for i = 1:length(id)
        eval(sprintf('results_%d = [];',i));
        % load xSPM and extract suprathresholded voxels
        for j = 1:runs
            if runs == 1
                j = single_run;
            end;
            
            %split1
            co_temp=[];
            stats_dir_filled = sprintf(split_dir,j);
            SPM_path = [stats_path f id{i} f stats_dir_filled];
            cd(box_path);
            fprintf('... create xSPM for %s split 1 in session %d...\n',id{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,1);
            eval(sprintf('co_sig%d = xSPM.XYZ;',j));

            if use_roi == 1    
                % apply ROI
                fprintf('...apply ROI...\n')
                for count_roi = 1:length(roi_xyz)
                    %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                    if idx>0
                        eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                       % fprintf('...significant...\n');
                    end;
                end;
                eval(sprintf('co_sig_split1_%d= co_temp;',j));
            else
                eval(sprintf('co_sig_split1_%d= co_sig%d;',j,j));
            end;
                
            eval(sprintf('nr_sig_split1_%d = size(co_sig_split1_%d,2);',j,j));
            clear co_temp;
            
            %split2
            co_temp=[];
            stats_dir_filled = sprintf(split_dir,j);
            SPM_path = [stats_path f id{i} f stats_dir_filled];
            cd(box_path);
            fprintf('... create xSPM for %s split 2 in session %d...\n',id{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,2+nr_para);
            eval(sprintf('co_sig%d = xSPM.XYZ;',j));
            
            if use_roi == 1
                % apply ROI
                for count_roi = 1:length(roi_xyz)
                    %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                    if idx>0
                        eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                       % fprintf('...significant...\n');
                    end;
                end;
                eval(sprintf('co_sig_split2_%d= co_temp;',j));
            else
                eval(sprintf('co_sig_split2_%d= co_sig%d;',j,j));
            end;                
            eval(sprintf('nr_sig_split2_%d = size(co_sig_split2_%d,2);',j,j));
            clear co_temp;
        end;

        % compare suprathresholded voxels  
        for k = 1:runs
            if runs == 1
                k = single_run;
            end;
            fprintf('...compare split1 to split2 in session %d...\n',k);
            count = 0;
            if  eval(sprintf('~isempty(co_sig_split2_%d)',k)) && eval(sprintf('~isempty(co_sig_split1_%d)',k))
                for l = 1:eval(sprintf('nr_sig_split1_%d',k))
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_split2_%d,co_sig_split1_%d(:,l))));',k,k));
                    count = count + length(idx);
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
                eval(sprintf('results_%d(1,end+1)=[jaccard_%d];',i,k))
            else
                eval(sprintf('overlap_%d = 0;',k));
                eval(sprintf('dice_%d = 0;',k));
                eval(sprintf('jaccard_%d = 0;',k));
                % subject result table
                eval(sprintf('results_%d(1,end+1)=[dice_%d];',i,k));
                eval(sprintf('results_%d(1,end+1)=[jaccard_%d];',i,k));
            end;
         end;
         eval(sprintf('results(i,:)=results_%d;',i));
         eval(sprintf('clear results_%d',i)); 
    end;
%create results table with names
row = id;
%generate column names
cols = {};
for k = 1:runs
    if runs == 1
        k = single_run;
    end;
    cols{1,end+1}=(sprintf('dice_%d;',k));
    cols{1,end+1}=(sprintf('jaccard_%d;',k));
end;
assignin('base','results',results);
assignin('base','cols',cols);
assignin('base','row',row);

results_overlap_subj_split=dataset({results,cols{:}}, ...
                'obsnames', row);
f2 = figure;            
imagesc(results,[0,1]);
colormap('inferno');
colorbar;
xlabel(sprintf('%s',cols{:}));
ylabel('subject ID');
title('results overlap subjectwise');
print(f2,['results_overlap_subjectwise_split' str],'-dpng','-r1200');

            
cd(results_dir);             
eval(sprintf('save results_overlap%s_subjectwise_split_%g.mat results_overlap_subj_split',str,p));       

if nr_para>0
    for ind_para=1:nr_para
        for i = 1:length(id)
            eval(sprintf('results_%d = [];',i));
            % load xSPM and extract suprathresholded voxels
            for j = 1:runs
                if runs == 1
                    j = single_run;
                end;
                %split1
                co_temp=[];
                stats_dir_filled = sprintf(split_dir,j);
                SPM_path = [stats_path f id{i} f stats_dir_filled];
                cd(box_path);
                fprintf('... create xSPM for %s split 1 in session %d...\n',id{i},j)

                xSPM = create_xSPM(SPM_path,box_path,p,1+ind_para);
                eval(sprintf('co_sig%d = xSPM.XYZ;',j));
            
                if use_roi == 1
                    % apply ROI
                    fprintf('...apply ROI...\n')
                    for count_roi = 1:length(roi_xyz)
                        %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                        eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                        if idx>0
                            eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                           % fprintf('...significant...\n');
                        end;
                    end;
                    eval(sprintf('co_sig_split1_%d= co_temp;',j));
                else
                    eval(sprintf('co_sig_split1_%d= co_sig%d;',j,j));               
                end;
                eval(sprintf('nr_sig_split1_%d = size(co_sig_split1_%d,2);',j,j));
                clear co_temp;
            
                %split2
                co_temp=[];
                stats_dir_filled = sprintf(split_dir,j);
                SPM_path = [stats_path f id{i} f stats_dir_filled];
                cd(box_path);
                fprintf('... create xSPM for %s split 2 in session %d...\n',id{i},j)

                xSPM = create_xSPM(SPM_path,box_path,p,2+ind_para+nr_para);
                eval(sprintf('co_sig%d = xSPM.XYZ;',j));

                if use_roi == 1
                    % apply ROI
                    for count_roi = 1:length(roi_xyz)
                        %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                        eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                        if idx>0
                            eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                           % fprintf('...significant...\n');
                        end;
                    end;
                    eval(sprintf('co_sig_split2_%d= co_temp;',j));
                else
                    eval(sprintf('co_sig_split2_%d= co_sig%d;',j,j));                
                end;
                eval(sprintf('nr_sig_split2_%d = size(co_sig_split2_%d,2);',j,j));
                clear co_temp;
            end;

            % compare suprathresholded voxels
        for k = 1:runs
            if runs == 1
                k = single_run;
            end;
            fprintf('...compare split1 to split2 in session %d...\n',k);
            count = 0;
            if  eval(sprintf('~isempty(co_sig_split2_%d)',k)) && eval(sprintf('~isempty(co_sig_split1_%d)',k))
                for l = 1:eval(sprintf('nr_sig_split1_%d',k))
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_split2_%d,co_sig_split1_%d(:,l))));',k,k));
                    count = count + length(idx);
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
            else
                eval(sprintf('overlap_%d = 0;',k));
                eval(sprintf('dice_%d = 0;',k));
                eval(sprintf('jaccard_%d = 0;',k));
                % subject result table
                eval(sprintf('results_%d(1,end+1)=[dice_%d];',i,k));
                eval(sprintf('results_%d(1,end+1)=[jaccard_%d];',i,k));
            end;
         end;
        eval(sprintf('results(i,:)=results_%d;',i));
        eval(sprintf('clear results_%d',i)); 
    end;

     %create results table with names
     row = id;
    %generate column names
    cols = {};
    for k = 1:runs
        if runs == 1
            k = single_run;
        end;
           cols{1,end+1}=(sprintf('dice_%d;',k));
           cols{1,end+1}=(sprintf('jaccard_%d;',k));
    end;
    assignin('base','results',results);
    assignin('base','cols',cols);
    assignin('base','row',row);


    results_overlap_subj_split=dataset({results,cols{:}}, ...
                    'obsnames', row);
    f2 = figure;            
    imagesc(results,[0,1]);
    colormap('inferno');
    colorbar;
    xlabel(sprintf('%s',cols{:}));
    ylabel('subject ID');
    title('results overlap subjectwise');
    print(f2,sprintf('results_overlap%s_subjectwise_split_par%d',str,ind_para),'-dpng','-r1200');

    cd(results_dir);             
    eval(sprintf('save results_overlap%s_subjectwise_split_par%d_%g.mat results_overlap_subj_split',str,ind_para,p));     
    end;
end;
disp('DONE');
cd(box_path);  

elseif two_cons == 1
    for i = 1:length(id)
        eval(sprintf('results_%d = [];',i));
        % load xSPM and extract suprathresholded voxels
        for j = 1:runs
            if runs == 1
                j = single_run;
            end;
            %con1
            co_temp=[];
            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = [stats_path f id{i} f stats_dir_filled];
            cd(box_path);
            fprintf('... create xSPM for %s con 1 in session %d...\n',id{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,con1_count);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);
            
            if use_roi == 1
                % apply ROI
                for count_roi = 1:length(roi_xyz)
                    %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                    if idx>0
                        eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                       % fprintf('...significant...\n');
                    end;
                end;
                eval(sprintf('co_sig_con1_%d= co_temp;',j));
            else
                eval(sprintf('co_sig_con1_%d= co_sig%d;',j,j));
            end;
                
            eval(sprintf('nr_sig_con1_%d = size(co_sig_con1_%d,2);',j,j));
            clear co_temp;    
      %con2
            co_temp=[];
            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = [stats_path f id{i} f stats_dir_filled];
            cd(box_path);
            fprintf('... create xSPM for %s con 2 in session %d...\n',id{i},j)

            xSPM = create_xSPM(SPM_path,box_path,p,con2_count);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);
            
            if use_roi == 1
                % apply ROI
                disp('...apply ROI...');
                for count_roi = 1:length(roi_xyz)
                    %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                    if idx>0
                        eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                       % fprintf('...significant...\n');
                    end;
                end;
                eval(sprintf('co_sig_con2_%d= co_temp;',j));
            else
                eval(sprintf('co_sig_con2_%d= co_sig%d;',j,j));                
            end;
                
            eval(sprintf('nr_sig_con2_%d = size(co_sig_con2_%d,2);',j,j));
            clear co_temp;  
        end;

        % compare suprathresholded voxels betweens cons in one session
        for k = 1:runs
            if runs == 1
                k = single_run;
            end;
                fprintf('...compare con1 to con2 in session %d...\n',k);
                count = 0;
                if eval(sprintf('nr_sig_con1_%d',k)) > 0 && eval(sprintf('nr_sig_con2_%d',k)) > 0
                for l = 1:eval(sprintf('nr_sig_con1_%d',k))
                   eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_con2_%d,co_sig_con1_%d(:,l))));',k,k));
                   count = count + length(idx);
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
        for ind_con = 1:2
            eval(sprintf('con = con%d;',ind_con));
            eval(sprintf('con_count = con%d_count;',ind_con));
        for k = runs:-1:2
                for ind = 1:k-1
                    fprintf('...compare %d to %d in %s...\n',k,k-ind, con);
                        count = 0;
                        if eval(sprintf('nr_sig_con%d_%d',ind_con,k)) > 0 && eval(sprintf('nr_sig_con%d_%d',ind_con,k-ind)) > 0
                    for l = 1:eval(sprintf('nr_sig_con%d_%d',ind_con,k))
                       eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_con%d_%d,co_sig_con%d_%d(:,l))));',ind_con,k-ind,ind_con,k));
                       count = count + length(idx);
                    end;
                        end;
                    eval(sprintf('%s_overlap_%d_%d = count;',con,k-ind,k));
                    if eval(sprintf('%s_overlap_%d_%d',con,k-ind,k)) == 0
                        eval(sprintf('%s_dice_%d_%d = 0;',con,k-ind,k));
                        eval(sprintf('%s_jaccard_%d_%d = 0;',con,k-ind,k));
                    else
                        eval(sprintf('%s_dice_%d_%d = (2.*%s_overlap_%d_%d)./(nr_sig_con%d_%d+nr_sig_con%d_%d);',con,k-ind,k,con,k-ind,k,ind_con,k-ind,ind_con,k));
                        eval(sprintf('%s_jaccard_%d_%d = %s_overlap_%d_%d./(nr_sig_con%d_%d+nr_sig_con%d_%d-%s_overlap_%d_%d);',con,k-ind,k,con,k-ind,k,ind_con,k-ind,ind_con,k,con,k-ind,k));
                    end;
                end;

        end;
        end;
        for ind_con = 1:2
            eval(sprintf('con = con%d;',ind_con));
            for k1 = runs:-1:2
                for ind1 = 1:k1-1 
                    % subject result table
                    eval(sprintf('results_%d(1,end+1)=[%s_dice_%d_%d];',i,con,k1-ind1,k1));
                    eval(sprintf('results_%d(1,end+1)=[%s_jaccard_%d_%d];',i,con,k1-ind1,k1));
                end;
            end;
        end;
            
            
        eval(sprintf('results(i,:)=results_%d;',i));
        eval(sprintf('clear results_%d',i)); 
    end;


 %create results table with names
    %id in rows
  
     row = id;
    cols = {};
    for k = 1:runs
        if runs == 1
            k = single_run;
        end;
           cols{1,end+1}=(sprintf('dice_betweenCons_%d;',k));
           cols{1,end+1}=(sprintf('jaccard_betweenCons_%d;',k));
    end;
       for ind_con = 1:2
            eval(sprintf('con = con%d;',ind_con));          
            for k1 = runs:-1:2
                for ind1 = 1:k1-1 
                    cols{1,end+1}=(sprintf('dice_%s_%d_%d;',con,k1,ind1));
                    cols{1,end+1}=(sprintf('jaccard_%s_%d_%d;',con,k1,ind1));    
                end;
            end;
        end;
    
    assignin('base','results',results);
    assignin('base','cols',cols);
    assignin('base','row',row);

    cd(results_dir);             

    results_overlap_subj=dataset({results,cols{:}}, ...
                    'obsnames', row);
    f2 = figure;            
    imagesc(results,[0,1]);
    colormap('inferno');
    colorbar;
    xlabel(sprintf('%s',cols{:}));
    ylabel('subject ID');
    title('results overlap subjectwise');
    print(f2,sprintf('results_overlap%s_%g_subjectwise',str,p),'-dpng','-r1200');

    eval(sprintf('save results%s_overlap_subjectwise_%g.mat results_overlap_subj',str,p)); 
    
    if nr_para1>0 && nr_para2>0
    for ind_para = 1:nr_para1
    for i = 1:length(id)
        eval(sprintf('results_%d = [];',i));
        % load xSPM and extract suprathresholded voxels
        for j = 1:runs
            if runs == 1
                j = single_run;
            end;
            %con1
            co_temp=[];
            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = [stats_path f id{i} f stats_dir_filled];
            cd(box_path);
            fprintf('... create xSPM for %s con 1 parametric %d in session %d...\n',id{i},ind_para,j)

            xSPM = create_xSPM(SPM_path,box_path,p,con1_count+ind_para);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);
            
            if use_roi == 1
                % apply ROI
                for count_roi = 1:length(roi_xyz)
                    %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                    if idx>0
                        eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                       % fprintf('...significant...\n');
                    end;
                end;
                eval(sprintf('co_sig_con1_%d= co_temp;',j));
            else
                eval(sprintf('co_sig_con1_%d= co_sig%d;',j,j));
            end;
                
            eval(sprintf('nr_sig_con1_%d = size(co_sig_con1_%d,2);',j,j));
            clear co_temp;    
      %con2
            co_temp=[];
            stats_dir_filled = sprintf(stats_dir,j);
            SPM_path = [stats_path f id{i} f stats_dir_filled];
            cd(box_path);
            fprintf('... create xSPM for %s con 2 parametric %d in session %d...\n',id{i},ind_para, j)

            xSPM = create_xSPM(SPM_path,box_path,p,con2_count+ind_para);
            evalstr1 = sprintf('co_sig%d = xSPM.XYZ;',j);
            eval(evalstr1);
            
            if use_roi == 1
                % apply ROI
                disp('...apply ROI...');
                for count_roi = 1:length(roi_xyz)
                    eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                    if idx>0
                        eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                    end;
                end;
                eval(sprintf('co_sig_con2_%d= co_temp;',j));
            else
                eval(sprintf('co_sig_con2_%d= co_sig%d;',j,j));                
            end;
                
            eval(sprintf('nr_sig_con2_%d = size(co_sig_con2_%d,2);',j,j));
            clear co_temp;  
        end;

        % compare suprathresholded voxels betweens cons in one session
        for k = 1:runs
            if runs == 1
                k = single_run;
            end;
                fprintf('...compare con1 to con2 parametric %d in session %d...\n',ind_para,k);
                count = 0;
                if eval(sprintf('nr_sig_con1_%d',k)) > 0 && eval(sprintf('nr_sig_con2_%d',k)) > 0
                for l = 1:eval(sprintf('nr_sig_con1_%d',k))
                   eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_con2_%d,co_sig_con1_%d(:,l))));',k,k));
                   count = count + length(idx);
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
        for ind_con = 1:2
            eval(sprintf('con = con%d;',ind_con));
            eval(sprintf('con_count = con%d_count;',ind_con));
        for k = runs:-1:2
                for ind = 1:k-1
                    fprintf('...compare %d to %d in %s...\n',k,k-ind, con);
                        count = 0;
                        if eval(sprintf('nr_sig_con%d_%d',ind_con,k)) > 0 && eval(sprintf('nr_sig_con%d_%d',ind_con,k-ind)) > 0
                    for l = 1:eval(sprintf('nr_sig_con%d_%d',ind_con,k))
                       eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig_con%d_%d,co_sig_con%d_%d(:,l))));',ind_con,k-ind,ind_con,k));
                       count = count + length(idx);
                    end;
                        end;
                    eval(sprintf('%s_overlap_%d_%d = count;',con,k-ind,k));
                    if eval(sprintf('%s_overlap_%d_%d',con,k-ind,k)) == 0
                        eval(sprintf('%s_dice_%d_%d = 0;',con,k-ind,k));
                        eval(sprintf('%s_jaccard_%d_%d = 0;',con,k-ind,k));
                    else
                        eval(sprintf('%s_dice_%d_%d = (2.*%s_overlap_%d_%d)./(nr_sig_con%d_%d+nr_sig_con%d_%d);',con,k-ind,k,con,k-ind,k,ind_con,k-ind,ind_con,k));
                        eval(sprintf('%s_jaccard_%d_%d = %s_overlap_%d_%d./(nr_sig_con%d_%d+nr_sig_con%d_%d-%s_overlap_%d_%d);',con,k-ind,k,con,k-ind,k,ind_con,k-ind,ind_con,k,con,k-ind,k));
                    end;
                end;

        end;
        end;
        for ind_con = 1:2
            eval(sprintf('con = con%d;',ind_con));
            for k1 = runs:-1:2
                for ind1 = 1:k1-1 
                    % subject result table
                    eval(sprintf('results_%d(1,end+1)=[%s_dice_%d_%d];',i,con,k1-ind1,k1));
                    eval(sprintf('results_%d(1,end+1)=[%s_jaccard_%d_%d];',i,con,k1-ind1,k1));
                end;
            end;
        end;
            
            
        eval(sprintf('results(i,:)=results_%d;',i));
        eval(sprintf('clear results_%d',i)); 
    end;


 %create results table with names
    %id in rows
  
     row = id;
    cols = {};
    for k = 1:runs
        if runs == 1
            k = single_run;
        end;
           cols{1,end+1}=(sprintf('dice_betweenCons_%d_par%d;',k,ind_para));
           cols{1,end+1}=(sprintf('jaccard_betweenCons_%d_par%d;',k,ind_para));
    end;
       for ind_con = 1:2
            eval(sprintf('con = con%d;',ind_con));          
            for k1 = runs:-1:2
                for ind1 = 1:k1-1 
                    cols{1,end+1}=(sprintf('dice_%s_%d_%d_par%d;',con,k1,ind1,ind_para));
                    cols{1,end+1}=(sprintf('jaccard_%s_%d_%d_par%d;',con,k1,ind1,ind_para));    
                end;
            end;
        end;
    
    assignin('base','results',results);
    assignin('base','cols',cols);
    assignin('base','row',row);

    cd(results_dir);             

    results_overlap_subj=dataset({results,cols{:}}, ...
                    'obsnames', row);
    f2 = figure;            
    imagesc(results,[0,1]);
    colormap('inferno');
    colorbar;
    xlabel(sprintf('%s',cols{:}));
    ylabel('subject ID');
    title('results overlap subjectwise');
    print(f2,sprintf('results_overlap%s_%g_subjectwise_par%d',str,p,ind_para),'-dpng','-r1200');
    end;
    end;
    eval(sprintf('save results%s_overlap_subjectwise_%g_par%d.mat results_overlap_subj',str,p,ind_para)); 
    clear results;
    disp('DONE');
    cd(box_path);    
end;
    

else
    disp('...calculate overlap on group level...');
    if two_cons == 0 && split == 0
    for j = 1:runs
        co_temp=[];
        cd(results_dir)
        % calculate T statistics via 4D images
        tmp = load_untouch_nii(sprintf('4D_%d.nii',j));
        tmp = tmp.img;
        if ~isa(tmp,'double')
            tmp = double(tmp);
        end
        % Get standard error of beta maps
            Beta_sd = std(tmp,0,4);  
            Beta_se = Beta_sd / sqrt(nr_subj);
        % Get mean over subjects of contrast beta maps
            Beta_mean = mean(tmp,4);
        % Get new T map
            Tmap = Beta_mean ./ Beta_se;
        % Get p-values
            Pmap = 1-tcdf(Tmap,nr_subj-1);
        [x,y,z] = ind2sub([size(Pmap,1),size(Pmap,2),size(Pmap,3)],find(Pmap<=p)); 
        x = x';y=y';z=z';
        eval(sprintf('co_sig%d = [x;y;z];',j));
        if use_roi == 1
        % apply ROI
        fprintf('...apply ROI...\n')
        for count_roi = 1:length(roi_xyz)
            %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
            eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
            if idx>0
                eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
            end;
        end;
        eval(sprintf('co_sig%d= co_temp;',j));
        clear co_temp;        
        end;
        eval(sprintf('nr_sig%d = size(co_sig%d,2);',j,j));
        
    end;
    
    % compare suprathresholded voxels
    for k = runs:-1:2
        for ind = 1:k-1
            fprintf('...compare %d to %d...\n',k,k-ind);
                count = 0;
                    for l = 1:eval(sprintf('nr_sig%d',k))
                           eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d(:,l),co_sig%d)));',k,k-ind));
                           count = count + length(idx);
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
    eval(sprintf('save results_overlap%s_grouplevel_%g.mat results_overlap_group',str,p)); 
    clear results_overlap_group
    
    if nr_para>0
        for ind_para=1:nr_para
            disp('...calculate overlap on group level...');

            for j = 1:runs
                co_temp=[];
                eval(sprintf('co_sig%d = [];',j));
                cd(results_dir)
                % calculate T statistics via 4D images
                tmp = load_untouch_nii(sprintf('4D_par%d_%d.nii',ind_para,j));
                tmp = tmp.img;
                if ~isa(tmp,'double')
                    tmp = double(tmp);
                end
                % Get standard error of beta maps
                    Beta_sd = std(tmp,0,4);  
                    Beta_se = Beta_sd / sqrt(nr_subj);
                % Get mean over subjects of contrast beta maps
                    Beta_mean = mean(tmp,4);
                % Get new T map
                    Tmap = Beta_mean ./ Beta_se;
                % Get p-values
                    Pmap = 1-tcdf(Tmap,nr_subj-1);
                [x,y,z] = ind2sub([size(Pmap,1),size(Pmap,2),size(Pmap,3)],find(Pmap<=p)); 
                x = x';y=y';z=z';
                eval(sprintf('co_sig%d = [x;y;z];',j));
            
                        % apply ROI
            if use_roi == 1
            fprintf('...apply ROI...\n')
            for count_roi = 1:length(roi_xyz)
                %fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d,roi_xyz(:,count_roi))));',j));
                if idx>0
                    eval(sprintf('co_temp(:,end+1)=co_sig%d(:,idx);',j));
                end;
            end;
            eval(sprintf('co_sig%d= co_temp;',j));
            clear co_temp;        
            end;
                eval(sprintf('nr_sig%d = size(co_sig%d,2);',j,j));
            end;

            % compare suprathresholded voxels
            for k = runs:-1:2
                for ind = 1:k-1
                    fprintf('...compare %d to %d...\n',k,k-ind);
                    count = 0;
                    for l = 1:eval(sprintf('nr_sig%d',k))
                        if eval(sprintf('nr_sig%d',k))> 0 && eval(sprintf('nr_sig%d',k-ind))> 0
                           eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d(:,l),co_sig%d)));',k,k-ind));
                           count = count + length(idx);
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
            eval(sprintf('save results_overlap%s_grouplevel_par%d_%g.mat results_overlap_group',str,ind_para,p));     
            clear results_overlap_group

        end;
    end;
    
    elseif split == 1
       for j = 1:runs
           for ind_split = 1:2
        co_temp=[];
        cd(results_dir)
        % calculate T statistics via 4D images
        tmp = load_untouch_nii(sprintf('4D_split%d_%d.nii',ind_split,j));
        tmp = tmp.img;
        if ~isa(tmp,'double')
            tmp = double(tmp);
        end
        % Get standard error of beta maps
            Beta_sd = std(tmp,0,4);  
            Beta_se = Beta_sd / sqrt(nr_subj);
        % Get mean over subjects of contrast beta maps
            Beta_mean = mean(tmp,4);
        % Get new T map
            Tmap = Beta_mean ./ Beta_se;
        % Get p-values
            Pmap = 1-tcdf(Tmap,nr_subj-1);
        [x,y,z] = ind2sub([size(Pmap,1),size(Pmap,2),size(Pmap,3)],find(Pmap<=p)); 
        x = x';y=y';z=z';
        eval(sprintf('co_sig%d_split%d = [x;y;z];',j,ind_split));
        % apply ROI
        if use_roi == 1
                    fprintf('...apply ROI...\n')

        for count_roi = 1:length(roi_xyz)
            fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
            for count_sig = 1:length(eval(sprintf('co_sig%d_split%d(1,:)',j,ind_split)))
                if isequal(eval(sprintf('co_sig%d_split%d(:,count_sig)',j,ind_split)),roi_xyz(:,count_roi))
                    eval(sprintf('co_temp(:,end+1)=co_sig%d_split%d(:,count_sig);',j,ind_split));
                    fprintf('...significant...\n');
                end;
            end;
        end;
        eval(sprintf('co_sig%d_split%d= co_temp;',j,ind_split));
        end;
        eval(sprintf('nr_sig%d_split%d = size(co_sig%d_split%d,2);',j,ind_split,j,ind_split));
        clear co_temp;        
            end;
       end;
    
    % compare suprathresholded voxels
    for k = 1:runs
            fprintf('...compare splits in session %d ...\n',k);
            count = 0;
                    for l = 1:eval(sprintf('nr_sig%d_split1',k))
                           eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d_split1(:,l),co_sig%d_split2)));',k,k));
                           count = count + length(idx);
                    end;
            eval(sprintf('overlap_split%d = count;',k));
            if eval(sprintf('overlap_split%d',k)) == 0
                eval(sprintf('dice_split%d = 0;',k));
                eval(sprintf('jaccard_split%d = 0;',k));
            else
                eval(sprintf('dice_split%d = (2.*overlap_split%d)./(nr_sig%d_split1+nr_sig%d_split2);',k,k,k,k));
                eval(sprintf('jaccard_split%d = overlap_split%d./(nr_sig%d_split1+nr_sig%d_split2-overlap_split%d);',k,k,k,k,k));
            end;
    end;
    
    %create results table
    results_overlap = [];
    for k1 = 1:runs
            eval(sprintf('results_overlap(1,end+1)=[dice_split%d];',k1));
            eval(sprintf('results_overlap(1,end+1)=[jaccard_split%d];',k1));
    end;
    %generate column names
    cols = {};
    for k2 = 1:runs
            
            % subject result table
                cols{1,end+1}=(sprintf('dice_split%d;',k2));
                cols{1,end+1}=(sprintf('jaccard_split%d;',k2));
    end;   
    results_overlap_group=dataset({results_overlap,cols{:}});
    cd(results_dir);             
    eval(sprintf('save results_overlap%s_grouplevel_split_%g.mat results_overlap_group',str,p)); 
    
    if nr_para>0
        for ind_para=1:nr_para
            disp('...calculate overlap on group level...');

            for j = 1:runs
                for ind_split = 1:2
                co_temp=[];
                cd(results_dir)
                % calculate T statistics via 4D images
                tmp = load_untouch_nii(sprintf('4D_split%d_par%d_%d.nii',ind_split,ind_para,j));
                tmp = tmp.img;
                if ~isa(tmp,'double')
                    tmp = double(tmp);
                end
                % Get standard error of beta maps
                    Beta_sd = std(tmp,0,4);  
                    Beta_se = Beta_sd / sqrt(nr_subj);
                % Get mean over subjects of contrast beta maps
                    Beta_mean = mean(tmp,4);
                % Get new T map
                    Tmap = Beta_mean ./ Beta_se;
                % Get p-values
                    Pmap = 1-tcdf(Tmap,nr_subj-1);
                [x,y,z] = ind2sub([size(Pmap,1),size(Pmap,2),size(Pmap,3)],find(Pmap<=p)); 
                x = x';y=y';z=z';
                eval(sprintf('co_sig%d_split%d = [x;y;z];',j,ind_split));
                % apply ROI
                if use_roi == 1
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
                eval(sprintf('co_sig%d_split%d= co_temp;',j,ind_split));
                end;
                eval(sprintf('nr_sig%d_split%d = size(co_sig%d_split%d,2);',j,ind_split,j,ind_split));
                clear co_temp;        
            end;

            % compare suprathresholded voxels
            for k = 1:runs
                    fprintf('...compare splits session %d ...\n',k);
                    count = 0;
                    for l = 1:eval(sprintf('nr_sig%d_split1',k))
                           eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d_split1(:,l),co_sig%d_split2)));',k,k));
                           count = count + length(idx);
                    end;
                    eval(sprintf('overlap_split%d = count;',k));
                    if eval(sprintf('overlap_split%d',k)) == 0
                        eval(sprintf('dice_split%d = 0;',k));
                        eval(sprintf('jaccard_split%d = 0;',k));
                    else
                        eval(sprintf('dice_split%d = (2.*overlap_split%d)./(nr_sig%d_split1+nr_sig%d_split2);',k,k,k,k));
                        eval(sprintf('jaccard_split%d = overlap_split%d./(nr_sig%d_split1+nr_sig%d_split2-overlap_split%d);',k,k,k,k,k));
                    end;


                end;


            end;

            %create results table
            results_overlap = [];
            for k1 = 1:3
                    eval(sprintf('results_overlap(1,end+1)=[dice_split%d];',k1));
                    eval(sprintf('results_overlap(1,end+1)=[jaccard_split%d];',k1));
            end;
            %generate column names
            cols = {};
            for k2 = 1:runs
                    % subject result table
                        cols{1,end+1}=(sprintf('dice_split%d;',k2));
                        cols{1,end+1}=(sprintf('jaccard_split%d;',k2));
            end;   
            results_overlap_group=dataset({results_overlap,cols{:}});
            cd(results_dir);             
            eval(sprintf('save results_overlap%s_grouplevel_split_par%d_%g.mat results_overlap_group',str,ind_para,p));     
        
        end;
    end; 
    elseif two_cons == 1
        for j = 1:runs
            for ind_con = 1:2
                eval(sprintf('con_temp = con%d;',ind_con));
                co_temp=[];
                cd(results_dir)
                % calculate T statistics via 4D images
                tmp = load_untouch_nii(sprintf('4D_%s_%d.nii',con_temp,j));
                tmp = tmp.img;
                if ~isa(tmp,'double')
                    tmp = double(tmp);
                end
                % Get standard error of beta maps
                    Beta_sd = std(tmp,0,4);  
                    Beta_se = Beta_sd / sqrt(nr_subj);
                % Get mean over subjects of contrast beta maps
                    Beta_mean = mean(tmp,4);
                % Get new T map
                    Tmap = Beta_mean ./ Beta_se;
                % Get p-values
                    Pmap = 1-tcdf(Tmap,nr_subj-1);
                [x,y,z] = ind2sub([size(Pmap,1),size(Pmap,2),size(Pmap,3)],find(Pmap<=p)); 
                x = x';y=y';z=z';
        eval(sprintf('co_sig%d_split%d = [x;y;z];',j,ind_con));
        % apply ROI
        if use_roi == 1
                    fprintf('...apply ROI...\n')

        for count_roi = 1:length(roi_xyz)
            fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
            for count_sig = 1:length(eval(sprintf('co_sig%d_split%d(1,:)',j,ind_con)))
                if isequal(eval(sprintf('co_sig%d_split%d(:,count_sig)',j,ind_con)),roi_xyz(:,count_roi))
                    eval(sprintf('co_temp(:,end+1)=co_sig%d_split%d(:,count_sig);',j,ind_con));
                    fprintf('...significant...\n');
                end;
            end;
        end;
        eval(sprintf('co_sig%d_split%d= co_temp;',j,ind_con));
        end;
        eval(sprintf('nr_sig%d_split%d = size(co_sig%d_split%d,2);',j,ind_con,j,ind_con));
        clear co_temp;        
            end;
       end;
    
    % compare suprathresholded voxels
    for k = 1:runs
            fprintf('...compare contrasts in session %d ...\n',k);
            count = 0;
                    for l = 1:eval(sprintf('nr_sig%d_split1',k))
                           eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d_split1(:,l),co_sig%d_split2)));',k,k));
                           count = count + length(idx);
                    end;
            eval(sprintf('overlap_between_%d = count;',k));
            if eval(sprintf('overlap_between_%d',k)) == 0
                eval(sprintf('dice_between_%d = 0;',k));
                eval(sprintf('jaccard_between_%d = 0;',k));
            else
                eval(sprintf('dice_between_%d = (2.*overlap_between_%d)./(nr_sig%d_split1+nr_sig%d_split2);',k,k,k,k));
                eval(sprintf('jaccard_between_%d = overlap_between_%d./(nr_sig%d_split1+nr_sig%d_split2-overlap_between_%d);',k,k,k,k,k));
            end;
    end;
    
    %create results table
    results_overlap = [];
    for k1 = 1:runs
            eval(sprintf('results_overlap(1,end+1)=[dice_between_%d];',k1));
            eval(sprintf('results_overlap(1,end+1)=[jaccard_between_%d];',k1));
    end;
    %generate column names
    cols = {};
    for k2 = 1:runs
            
            % subject result table
                cols{1,end+1}=(sprintf('dice_between_%d;',k2));
                cols{1,end+1}=(sprintf('jaccard_between_%d;',k2));
    end;   
    
     % compare suprathresholded voxels between sessions same contrast
     for ind_con = 1:2
            eval(sprintf('con_temp = con%d;',ind_con));
         
            for k = runs:-1:2
                for ind = 1:k-1
                    fprintf('...compare %d to %d...\n',k,k-ind);
                    count = 0;
                    for l = 1:eval(sprintf('nr_sig%d_split%d',k,ind_con))
                           eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d_split%d(:,l),co_sig%d_split%d)));',k,ind_con,k-ind,ind_con));
                           count = count + length(idx);
                    end;
                    eval(sprintf('overlap_%d_%d_%s = count;',k-ind,k,con_temp));
                    if eval(sprintf('overlap_%d_%d_%s',k-ind,k,con_temp)) == 0
                        eval(sprintf('dice_%d_%d_%s = 0;',k-ind,k,con_temp));
                        eval(sprintf('jaccard_%d_%d_%s = 0;',k-ind,k,con_temp));
                    else
                        eval(sprintf('dice_%d_%d_%s = (2.*overlap_%d_%d_%s)./(nr_sig%d_split%d+nr_sig%d_split%d);',k-ind,k,con_temp,k-ind,k,con_temp,k-ind,ind_con,k,ind_con));
                        eval(sprintf('jaccard_%d_%d_%s = overlap_%d_%d_%s./(nr_sig%d_split%d+nr_sig%d_split%d-overlap_%d_%d_%s);',k-ind,k,con_temp,k-ind,k,con_temp,k-ind,ind_con,k,ind_con,k-ind,k,con_temp));
                    end;


                end;


            end;

            %create results table
            for k1 = runs:-1:2
                for ind1 = 1:k1-1 
                    eval(sprintf('results_overlap(1,end+1)=[dice_%d_%d_%s];',k1-ind1,k1,con_temp));
                    eval(sprintf('results_overlap(1,end+1)=[jaccard_%d_%d_%s];',k1-ind1,k1,con_temp));
                end;
            end;
            %generate column names
            for k2 = runs:-1:2
                    for ind2 = 1:k2-1 
                    % subject result table
                        cols{1,end+1}=(sprintf('dice_%d_%d_%s;',k2-ind2,k2,con_temp));
                        cols{1,end+1}=(sprintf('jaccard_%d_%d_%s;',k2-ind2,k2,con_temp));
                    end;
            end;   
    
     end;
   
    
    results_overlap_group=dataset({results_overlap,cols{:}});
    cd(results_dir);             
    eval(sprintf('save results_overlap%s_grouplevel_%g.mat results_overlap_group',str,p)); 
    
    if nr_para1>0 && nr_para2>0
        for ind_para = 1:nr_para1
           for j = 1:runs
                for ind_con = 1:2
                    eval(sprintf('con_temp = con%d;',ind_con));
                    co_temp=[];
                    cd(results_dir)
                    % calculate T statistics via 4D images
                    tmp = load_untouch_nii(sprintf('4D_%s_par%d_%d.nii',con_temp,ind_para,j));
                    tmp = tmp.img;
                    if ~isa(tmp,'double')
                        tmp = double(tmp);
                    end
                    % Get standard error of beta maps
                        Beta_sd = std(tmp,0,4);  
                        Beta_se = Beta_sd / sqrt(nr_subj);
                    % Get mean over subjects of contrast beta maps
                        Beta_mean = mean(tmp,4);
                    % Get new T map
                        Tmap = Beta_mean ./ Beta_se;
                    % Get p-values
                        Pmap = 1-tcdf(Tmap,nr_subj-1);
                    [x,y,z] = ind2sub([size(Pmap,1),size(Pmap,2),size(Pmap,3)],find(Pmap<=p)); 
                    x = x';y=y';z=z';
            eval(sprintf('co_sig%d_split%d = [x;y;z];',j,ind_con));
            % apply ROI
            if use_roi == 1
                        fprintf('...apply ROI...\n')

            for count_roi = 1:length(roi_xyz)
                fprintf('...voxel x = %d .. y = %d .. z = %d ...\n',roi_xyz(1,count_roi),roi_xyz(2,count_roi),roi_xyz(3,count_roi));
                for count_sig = 1:length(eval(sprintf('co_sig%d_split%d(1,:)',j,ind_con)))
                    if isequal(eval(sprintf('co_sig%d_split%d(:,count_sig)',j,ind_con)),roi_xyz(:,count_roi))
                        eval(sprintf('co_temp(:,end+1)=co_sig%d_split%d(:,count_sig);',j,ind_con));
                        fprintf('...significant...\n');
                    end;
                end;
            end;
            eval(sprintf('co_sig%d_split%d= co_temp;',j,ind_con));
            end;
            eval(sprintf('nr_sig%d_split%d = size(co_sig%d_split%d,2);',j,ind_con,j,ind_con));
            clear co_temp;        
                end;
           end;

        % compare suprathresholded voxels
        for k = 1:runs
                fprintf('...compare contrasts in session %d ...\n',k);
                count = 0;
                        for l = 1:eval(sprintf('nr_sig%d_split1',k))
                               eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d_split1(:,l),co_sig%d_split2)));',k,k));
                               count = count + length(idx);
                        end;
                eval(sprintf('overlap_between_%d = count;',k));
                if eval(sprintf('overlap_between_%d',k)) == 0
                    eval(sprintf('dice_between_%d = 0;',k));
                    eval(sprintf('jaccard_between_%d = 0;',k));
                else
                    eval(sprintf('dice_between_%d = (2.*overlap_between_%d)./(nr_sig%d_split1+nr_sig%d_split2);',k,k,k,k));
                    eval(sprintf('jaccard_between_%d = overlap_between_%d./(nr_sig%d_split1+nr_sig%d_split2-overlap_between_%d);',k,k,k,k,k));
                end;
        end;

        %create results table
        results_overlap = [];
        for k1 = 1:runs
                eval(sprintf('results_overlap(1,end+1)=[dice_between_%d];',k1));
                eval(sprintf('results_overlap(1,end+1)=[jaccard_between_%d];',k1));
        end;
        %generate column names
        cols = {};
        for k2 = 1:runs

                % subject result table
                    cols{1,end+1}=(sprintf('dice_between_%d;',k2));
                    cols{1,end+1}=(sprintf('jaccard_between_%d;',k2));
        end;   

         % compare suprathresholded voxels between sessions same contrast
         for ind_con = 1:2
                eval(sprintf('con_temp = con%d;',ind_con));

                for k = runs:-1:2
                    for ind = 1:k-1
                        fprintf('...compare %d to %d...\n',k,k-ind);
                        count = 0;
                        for l = 1:eval(sprintf('nr_sig%d_split%d',k,ind_con))
                               eval(sprintf('idx=find(~any(bsxfun(@minus,co_sig%d_split%d(:,l),co_sig%d_split%d)));',k,ind_con,k-ind,ind_con));
                               count = count + length(idx);
                        end;
                        eval(sprintf('overlap_%d_%d_%s = count;',k-ind,k,con_temp));
                        if eval(sprintf('overlap_%d_%d_%s',k-ind,k,con_temp)) == 0
                            eval(sprintf('dice_%d_%d_%s = 0;',k-ind,k,con_temp));
                            eval(sprintf('jaccard_%d_%d_%s = 0;',k-ind,k,con_temp));
                        else
                            eval(sprintf('dice_%d_%d_%s = (2.*overlap_%d_%d_%s)./(nr_sig%d_split%d+nr_sig%d_split%d);',k-ind,k,con_temp,k-ind,k,con_temp,k-ind,ind_con,k,ind_con));
                            eval(sprintf('jaccard_%d_%d_%s = overlap_%d_%d_%s./(nr_sig%d_split%d+nr_sig%d_split%d-overlap_%d_%d_%s);',k-ind,k,con_temp,k-ind,k,con_temp,k-ind,ind_con,k,ind_con,k-ind,k,con_temp));
                        end;


                    end;


                end;

                %create results table
                for k1 = runs:-1:2
                    for ind1 = 1:k1-1 
                        eval(sprintf('results_overlap(1,end+1)=[dice_%d_%d_%s];',k1-ind1,k1,con_temp));
                        eval(sprintf('results_overlap(1,end+1)=[jaccard_%d_%d_%s];',k1-ind1,k1,con_temp));
                    end;
                end;
                %generate column names
                for k2 = runs:-1:2
                        for ind2 = 1:k2-1 
                        % subject result table
                            cols{1,end+1}=(sprintf('dice_%d_%d_%s;',k2-ind2,k2,con_temp));
                            cols{1,end+1}=(sprintf('jaccard_%d_%d_%s;',k2-ind2,k2,con_temp));
                        end;
                end;   

         end;


        results_overlap_group=dataset({results_overlap,cols{:}});
        cd(results_dir);             
        eval(sprintf('save results_overlap%s_grouplevel_%g_par%d.mat results_overlap_group',str,p,ind_para)); 

        end;
    end;
    
    
    end;
        disp('DONE');
    cd(box_path);

end;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
axes(hObject);
imshow('logo.png');
