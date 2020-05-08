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

% Last Modified by GUIDE v2.5 24-Sep-2018 14:10:35

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
runs = study_design.number_sessions;
nr_subj = str2double(study_design.number_subjects);
exStats = study_design.exist_stats;
ex4D = study_design.exist_4D;
if exStats == 1
    load(study_design.subject_list);
    stats = study_design.stats_directory;
    path = study_design.stats_path;
elseif ex4D == 1
    nr_cond = contrast_def.number_conditions;
    if nr_cond > 1
        conditions = contrast_def.conditions;
    end;
end;

%% get GUI input
sim2mean = get(handles.sim2mean,'Value');
split = get(handles.split,'Value');
use_roi = get(handles.use_roi,'Value');
if use_roi == 1
    str = 'ROI';
else
    str = '';
end;

%% load contrast information
if exStats == 1
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
        nr_para1 = study_design.number_parametric1;
        nr_para2 = study_design.number_parametric2;

    end;   
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
    if exStats == 1
        stats_filled = sprintf(stats,1);
        temp = [path f id{1} f stats_filled f con ',1'];
    elseif ex4D == 1
        cd(results_dir)
        abk_4Dto3D([results_dir f '4D_1.nii'],1)
        temp = [results_dir f 'template_3D.nii' ',1'];
    end;
        
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
        r_roi = load_untouch_nii(sprintf('r%s.img',roi_pur));
        r_roi_ind = r_roi.img>0.0001;
    else
        if ~strcmp(roi_dir,results_dir)
        compl = [roi_dir f  'r' roi_pur '.nii'];
        movefile(compl,results_dir,'f');
        end;
        cd(results_dir);    
        r_roi = load_untouch_nii(sprintf('r%s',roi_name));
        r_roi_ind = r_roi.img>0.0001;        
    end;
else
    str = '';
    if exStats == 1
        stats_filled = sprintf(stats,1);
        temp = [path f id{1} f stats_filled f con];
        temp = load_untouch_nii(temp);
    elseif ex4D == 1
        temp = [results_dir f 'template_3D.nii'];
        temp = load_untouch_nii(temp);
    end;
    [x,y,z] = size(temp.img);
    r_roi_ind = zeros(x,y,z);
    r_roi_ind = r_roi_ind==0;
end;


%% load 4D images
if exStats == 1
if two_cons == 0 && split == 0
    for ind_run = 1:runs
        fprintf('...load 4D image for run %d...\n',ind_run);
        file = [results_dir f '4D_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file;',ind_run));
        eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d = FourD%d.img;',ind_run,ind_run));
        eval(sprintf('dims = size(FourD%d(:,:,:,1));',ind_run));
        x = dims(1);
        y = dims(2);
        z = dims(3);
    end
    if ~isa(FourD1,'double')
        for ind_run = 1:runs
            eval(sprintf('FourD%d = double(FourD%d);',ind_run,ind_run));
        end    
    end
    % load parametric modulator
    if nr_para > 0
        for ind_para = 1:nr_para
            for ind_run = 1:runs
                fprintf('...load 4D parametric image for run %d...\n',ind_run);
                file = [results_dir f '4D_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d_par = file;',ind_run));
                eval(sprintf('FourD%d_par = load_untouch_nii(FourD%d_par);',ind_run,ind_run));
                eval(sprintf('FourD%d_par%d = FourD%d_par.img;',ind_run,ind_para,ind_run));

            end
        end
        if ~isa(FourD1_par1,'double')
            for ind_para = 1:nr_para
                for ind_run = 1:runs
                    eval(sprintf('FourD%d_par%d = double(FourD%d_par%d);',ind_run,ind_para,ind_run,ind_para));
                end    
            end
        end       
    end
    
elseif two_cons == 1
     for ind_run = 1:runs
        fprintf('...load 4D image for run %d and contrast %s...\n',ind_run,con1);
        file1 = [results_dir f '4D_' con1 '_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file1;',ind_run));
        eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d_%d = FourD%d.img;',ind_run,con1_count,ind_run));
        eval(sprintf('dims = size(FourD%d_%d(:,:,:,1));',ind_run,con1_count));
        x = dims(1);
        y = dims(2);
        z = dims(3);

        fprintf('...load 4D image for run %d and contrast %s...\n',ind_run,con2);
        file2 = [results_dir f '4D_' con2 '_' num2str(ind_run) '.nii'];
        eval(sprintf('FourD%d = file2;',ind_run));
        eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d_%d = FourD%d.img;',ind_run,con2_count,ind_run));
     end
    eval(sprintf('temp = FourD1_%d;',con2_count));    
    if ~isa(temp,'double')
        for ind_run = 1:runs
            eval(sprintf('FourD%d_%d = double(FourD%d_%d);',ind_run,con1_count,ind_run,con1_count));
            eval(sprintf('FourD%d_%d = double(FourD%d_%d);',ind_run,con2_count,ind_run,con2_count));
        end    
    end
    if nr_para1 > 0
        for ind_para = 1:nr_para1
            for ind_run = 1:runs
                fprintf('...load 4D image parametric for run %d and contrast %s...\n',ind_run,con1);
                file1 = [results_dir f '4D_' con1 '_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file1;',ind_run));
                eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d_%d_par%d = FourD%d.img;',ind_run,con1_count,ind_para,ind_run));
            end
        end
        eval(sprintf('temp = FourD1_%d_par1;',con1_count));    
        if ~isa(temp,'double')
            for ind_para = 1:nr_para1
                for ind_run = 1:runs
                    eval(sprintf('FourD%d_%d_par%d = double(FourD%d_%d_par%d);',ind_run,con1_count,ind_para,ind_run,con1_count,ind_para));
                end    
            end
        end
    end
    
    if nr_para2 > 0
        for ind_para = 1:nr_para2
            for ind_run = 1:runs    
                fprintf('...load 4D image parametric for run %d and contrast %s...\n',ind_run,con2);
                file2 = [results_dir f '4D_' con2 '_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file2;',ind_run));
                eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d_%d_par%d = FourD%d.img;',ind_run,con2_count,ind_para,ind_run));
            end
        end
        eval(sprintf('temp = FourD1_%d_par1;',con2_count));    
        if ~isa(temp,'double')
            for ind_para = 1:nr_para2
                for ind_run = 1:runs
                    eval(sprintf('FourD%d_%d_par%d = double(FourD%d_%d_par%d);',ind_run,con2_count,ind_para,ind_run,con2_count,ind_para));
                end    
            end
        end        
    end
elseif split == 1
    for ind_run = 1:runs
            fprintf('...load 4D image for run %d and split 1...\n',ind_run);
            file1 = [results_dir f '4D_split1_' num2str(ind_run) '.nii'];
            eval(sprintf('FourD%d = file1;',ind_run));
            eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
            eval(sprintf('FourD%d_split1 = FourD%d.img;',ind_run,ind_run));
            eval(sprintf('dims = size(FourD%d_split1(:,:,:,1));',ind_run));
            x = dims(1);
            y = dims(2);
            z = dims(3);

            fprintf('...load 4D image for run %d and split 2...\n',ind_run);
            file2 = [results_dir f '4D_split2_' num2str(ind_run) '.nii'];
            eval(sprintf('FourD%d = file2;',ind_run));
            eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
            eval(sprintf('FourD%d_split2 = FourD%d.img;',ind_run,ind_run));

    end
    if ~isa(FourD1_split1,'double')
        for ind_run = 1:runs
            eval(sprintf('FourD%d_split1 = double(FourD%d_split1);',ind_run,ind_run));
            eval(sprintf('FourD%d_split2 = double(FourD%d_split2);',ind_run,ind_run));
        end    
    end
if nr_para > 0 
    for ind_para = 1:nr_para
        for ind_run = 1:runs
                fprintf('...load 4D image parametric for run %d and split 1...\n',ind_run);
                file1 = [results_dir f '4D_split1_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];               
                eval(sprintf('FourD%d = file1;',ind_run));
                eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d_split1_par%d = FourD%d.img;',ind_run,ind_para,ind_run));

                fprintf('...load 4D image parametric for run %d and split 2...\n',ind_run);
                file2 = [results_dir f '4D_split2_par' num2str(ind_para) '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file2;',ind_run));
                eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d_split2_par%d = FourD%d.img;',ind_run,ind_para,ind_run));
                
        end
    end
    if ~isa(FourD1_split1_par1,'double')
        for ind_para = 1:nr_para
            for ind_run = 1:runs
                eval(sprintf('FourD%d_split1_par%d = double(FourD%d_split1_par%d);',ind_run,ind_para,ind_run,ind_para));
                eval(sprintf('FourD%d_split2_par%d = double(FourD%d_split2_par%d);',ind_run,ind_para,ind_run,ind_para));
            end    
        end
    end
end
end
elseif ex4D == 1
    for ind_cond = 1:nr_cond
        for ind_run = 1:runs
            fprintf('...load 4D image for run %d...\n',ind_run);
            if nr_cond > 1
                fprintf('...condition %s...\n',conditions{ind_cond,1})
            end
            
            if nr_cond == 1
                file = [results_dir f '4D_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file;',ind_run));
                eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d = FourD%d.img;',ind_run,ind_run));

                %calculation of mean activation 
                eval(sprintf('dims = size(FourD%d(:,:,:,1));',ind_run));
                x = dims(1);
                y = dims(2);
                z = dims(3);

            else
                file = [results_dir f '4D_' conditions{ind_cond,1} '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file;',ind_run));
                eval(sprintf('FourD%d = load_untouch_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('check_double = FourD%d.img;',ind_run))
                if ~isa(check_double,'double')
                    eval(sprintf('FourD%d.img = double(FourD%d.img)',ind_run,ind_run))
                end                
                eval(sprintf('FourD%s%d = FourD%d.img;',conditions{ind_cond,1},ind_run,ind_run));
             end
        end
    end
end

%% calculation of similarity
if exStats == 1
if two_cons == 0 && split == 0
for ind_run = 1:runs
    for count = 0:runs-1
        if ind_run+count <= runs
                fprintf('...compare session %d to %d...\n',ind_run,ind_run+count);
                eval(sprintf('TempFourD1 = FourD%d;',ind_run));
                eval(sprintf('TempFourD2 = FourD%d;',ind_run+count));
                out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
        
                cd(results_dir);
                eval(sprintf('save similarity-%s-%s-%d-%d.mat out',str,con,ind_run, ind_run+count));           
                
               [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
               fig.PaperPositionMode = 'auto';
               tempc = strsplit(con,'.');
               print(fig,sprintf('similarity-%s-%s-%d-%d',str,tempc{1,1},ind_run, ind_run+count),'-dpng','-r1200')

                close(fig);
                clearvars out fig TempFourD1 TempFourD2 ;
        end;
    end;
end;
if nr_para > 0 
    for ind_para = 1:nr_para
        for ind_run = 1:runs
            for count = 0:runs-1
                if ind_run+count <= runs
                        fprintf('...compare parametric session %d to %d...\n',ind_run,ind_run+count);
                        eval(sprintf('TempFourD1 = FourD%d_par%d;',ind_run,ind_para));
                        eval(sprintf('TempFourD2 = FourD%d_par%d;',ind_run+count,ind_para));
                        out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
                        
                        cd(results_dir);
                        eval(sprintf('save similarity-par%d-%s-%s-%d-%d.mat out',ind_para,str,con,ind_run, ind_run+count));           

                        [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
                        fig.PaperPositionMode = 'auto';
                        tempc = strsplit(con,'.');
                        print(fig,sprintf('similarity-%s-%s-par%d-%d-%d',str,tempc{1,1},ind_para,ind_run, ind_run+count),'-dpng','-r1200')

                        close(fig);
                        clearvars out fig TempFourD1 TempFourD2 ;

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
                eval(sprintf('TempFourD1 = FourD%d_%d;',ind_run,con1_count));
                eval(sprintf('TempFourD2 = FourD%d_%d;',ind_run+count,con1_count));
                out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
                
                cd(results_dir);
                eval(sprintf('save similarity-%s-%s-%d-%d.mat out',str,con1,ind_run, ind_run+count));           
                
                [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
                fig.PaperPositionMode = 'auto';
                tempc = strsplit(con1,'.');
                print(fig,sprintf('similarity-%s-%s-%d-%d',str,tempc{1,1},ind_run, ind_run+count),'-dpng','-r1200')

                close(fig);
                clearvars out fig TempFourD1 TempFourD2 ;
                
                %con1 parametric
                if nr_para1 > 0
                    for ind_para = 1:nr_para1
                        fprintf('...compare session %d to %d in parametric modulator %d of contrast %d...\n',ind_run,ind_run+count,ind_para,con1_count);
                        eval(sprintf('TempFourD1 = FourD%d_%d_par%d;',ind_run,con1_count,ind_para));
                        eval(sprintf('TempFourD2 = FourD%d_%d_par%d;',ind_run+count,con1_count,ind_para));
                        out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);

                        cd(results_dir);
                        eval(sprintf('save similarity-%s-%s-par%d-%d-%d.mat out',str,con1,ind_para,ind_run, ind_run+count));           

                        [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
                        fig.PaperPositionMode = 'auto';
                        tempc = strsplit(con1,'.');
                        print(fig,sprintf('similarity-%s-%s-par%-d-%d-%d',str,tempc{1,1},ind_para,ind_run, ind_run+count),'-dpng','-r1200')

                        close(fig);
                        clearvars out fig TempFourD1 TempFourD2 ;    
                    end
                end

                %con2
                fprintf('...compare session %d to %d in contrast %d...\n',ind_run,ind_run+count,con2_count);
                eval(sprintf('TempFourD1 = FourD%d_%d;',ind_run,con2_count));
                eval(sprintf('TempFourD2 = FourD%d_%d;',ind_run+count,con2_count));
                out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
                
                cd(results_dir);
                eval(sprintf('save similarity-%s-%s-%d-%d.mat out',str,con2,ind_run, ind_run+count));           
                
                [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
                fig.PaperPositionMode = 'auto';
                tempc = strsplit(con2,'.');
                print(fig,sprintf('similarity-%s-%s-%d-%d',str,tempc{1,1},ind_run, ind_run+count),'-dpng','-r1200')

                close(fig);
                clearvars out fig TempFourD1 TempFourD2 ;

                %con2 parametric
                if nr_para2 > 0
                    for ind_para = 1:nr_para2
                        fprintf('...compare session %d to %d in parametric modulator %d of contrast %d...\n',ind_run,ind_run+count,ind_para,con2_count);
                        eval(sprintf('TempFourD1 = FourD%d_%d_par%d;',ind_run,con2_count,ind_para));
                        eval(sprintf('TempFourD2 = FourD%d_%d_par%d;',ind_run+count,con2_count,ind_para));
                        out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);

                        cd(results_dir);
                        eval(sprintf('save similarity-%s-%s-par%d-%d-%d.mat out',str,con2,ind_para,ind_run, ind_run+count));           

                        [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
                        fig.PaperPositionMode = 'auto';
                        tempc = strsplit(con2,'.');
                        print(fig,sprintf('similarity-%s-%s-par%-d-%d-%d',str,tempc{1,1},ind_para,ind_run, ind_run+count),'-dpng','-r1200')

                        close(fig);
                        clearvars out fig TempFourD1 TempFourD2 ;    
                    end;
                end;
            end;
        end;
    end;

    % compares two contrasts within one session
    for  i_run = 1:runs
        fprintf('...compare session %d in contrast %d and %d...\n',i_run,con1_count,con2_count);
        eval(sprintf('TempFourD1 = FourD%d_%d;',i_run,con1_count));
        eval(sprintf('TempFourD2 = FourD%d_%d;',i_run,con2_count));
        out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
                
        cd(results_dir);
        eval(sprintf('save similarity-%s-%s-%s-%d.mat out',str,con1,con2,ind_run));           
                
        [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run,str);
        fig.PaperPositionMode = 'auto';
        tempc1 = strsplit(con1,'.');
        tempc2 = strsplit(con2,'.');
        print(fig,sprintf('similarity-%s-%s-%s-%d',str,tempc1{1,1},tempc2{1,1},ind_run),'-dpng','-r1200')

        close(fig);
        clearvars out fig TempFourD1 TempFourD2 ;

        if nr_para1 > 0 && nr_para2 > 0
            for ind_para = 1:nr_para1
                eval(sprintf('TempFourD1 = FourD%d_%d_par%d;',i_run,con1_count,ind_para));
                eval(sprintf('TempFourD2 = FourD%d_%d_par%d;',i_run,con2_count,ind_para));
                out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);

                cd(results_dir);
                eval(sprintf('save similarity-%s-%s-%s-par%d-%d.mat out',str,con1,con2,ind_para,ind_run));           

                [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run,str);
                fig.PaperPositionMode = 'auto';
                tempc1 = strsplit(con1,'.');
                tempc2 = strsplit(con2,'.');                
                print(fig,sprintf('similarity-%s-%s-%s-par%d-%d',str,tempc1{1,1},tempc2{1,1},ind_para,ind_run),'-dpng','-r1200')

                close(fig);
                clearvars out fig TempFourD1 TempFourD2 ;
            end;          
        end;
                
  
    end;

elseif split == 1
    for  i_run = 1:runs
        fprintf('...compare splits session %d ...\n',i_run);
        eval(sprintf('TempFourD1 = FourD%d_split1;',i_run));
        eval(sprintf('TempFourD2 = FourD%d_split2;',i_run));
        out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
        
        cd(results_dir);
        eval(sprintf('save similarity-%s-%s-%d-split.mat out',str,con,i_run));           
        
        [fig] = similarity_figure(out.r_mat,sim2mean,i_run,i_run,str);
        fig.PaperPositionMode = 'auto';
        tempc = strsplit(con,'.');                
        print(fig,sprintf('similarity-%s-%s-%d-split',str,tempc{1,1},i_run),'-dpng','-r1200')

        close(fig);
        clearvars out fig TempFourD1 TempFourD2 ;
        
        if nr_para > 0
            for ind_para = 1:nr_para
                fprintf('...compare splits parametric modulator %d session %d ...\n',ind_para,i_run);
                eval(sprintf('TempFourD1 = FourD%d_split1_par%d(:,:,:,i);',i_run,ind_para));
                eval(sprintf('TempFourD2 = FourD%d_split2_par%d(:,:,:,j);',i_run,ind_para));
                out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
                cd(results_dir);
                eval(sprintf('save similarity-par%d-%s-%s-%d-split.mat out',ind_para,str,con,i_run));           
                
                [fig] = similarity_figure(out.r_mat,sim2mean,i_run,i_run,str);
                fig.PaperPositionMode = 'auto';
                tempc = strsplit(con,'.');
                print(fig,sprintf('similarity-%s-%s-par%d-%d-split',str,tempc{1,1},ind_para,i_run),'-dpng','-r1200')

                close(fig);
                clearvars out fig TempFourD1 TempFourD2 ;
            end;
         end;           
    end;
end;
elseif ex4D == 1
if nr_cond == 1
for ind_run = 1:runs
    for count = 0:runs-1
        if ind_run+count <= runs
                fprintf('...compare session %d to %d...\n',ind_run,ind_run+count);
                eval(sprintf('TempFourD1 = FourD%d;',ind_run));
                eval(sprintf('TempFourD2 = FourD%d;',ind_run+count));
                out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
        
                cd(results_dir);
                eval(sprintf('save similarity-%s-%d-%d.mat out',str,ind_run, ind_run+count));           
                
               [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
               fig.PaperPositionMode = 'auto';
               
               print(fig,sprintf('similarity-%s-%d-%d',str,ind_run, ind_run+count),'-dpng','-r1200')

                close(fig);
                clearvars out fig TempFourD1 TempFourD2 ;
            
        end;
    end;
end;
else
    for ind_cond = 1:nr_cond
    for ind_run = 1:runs
        for count = 0:runs-1
            if ind_run+count <= runs
                    fprintf('...compare session %d to %d in condition %s...\n',ind_run,ind_run+count,conditions{ind_cond,1});
                    eval(sprintf('TempFourD1 = FourD%s%d;',conditions{ind_cond,1},ind_run));
                    eval(sprintf('TempFourD2 = FourD%s%d;',conditions{ind_cond,1},ind_run+count));
                    out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
                    cd(results_dir);
                    eval(sprintf('save similarity-%s-%d-%d-%s.mat out',str,ind_run, ind_run+count,conditions{ind_cond,1}));           
                    [fig] = similarity_figure(out.r_mat,sim2mean,ind_run,ind_run+count,str);
                    fig.PaperPositionMode = 'auto';
                    print(fig,sprintf('similarity-%s-%d-%d-%s',str,ind_run, ind_run+count,conditions{ind_cond,1}),'-dpng','-r1200')

                    close(fig);
                    clearvars out fig TempFourD1 TempFourD2 ;

             end;
        end;
    end;
    end;
            
    % compares two conditions within one session
    for  i_run = 1:runs
        fprintf('...compare session %d in condition %s and %s...\n',i_run,conditions{1,1},conditions{2,1});
        eval(sprintf('TempFourD1 = FourD%s%d(:,:,:,i);',conditions{1,1},i_run));
        eval(sprintf('TempFourD2 = FourD%s%d(:,:,:,j);',conditions{2,1},i_run));
        out = similarity_subjectwise(nr_subj,TempFourD1,TempFourD2,use_roi,sim2mean,r_roi_ind);
        cd(results_dir);
        eval(sprintf('save similarity-%s-%d-%s-%s.mat out',str,i_run,conditions{1,1},conditions{2,1}));           
        [fig] = similarity_figure(out.r_mat,sim2mean,i_run,ind_run+count,str);
        fig.PaperPositionMode = 'auto';
        print(fig,sprintf('similarity-%s-%d-%s-%s',str,i_run,conditions{1,1},conditions{2,1}),'-dpng','-r1200')

        close(fig);
        clearvars out fig TempFourD1 TempFourD2 ;
        
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


% --- Executes on button press in sim2mean.
function sim2mean_Callback(hObject, eventdata, handles)
% hObject    handle to sim2mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim2mean
