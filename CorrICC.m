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

% Last Modified by GUIDE v2.5 11-Sep-2017 15:46:11

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
set(handles.roi_name,'TooltipString','Type in name of ROI file WITHOUT suffix');

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
box_path=evalin('base','box_path');

%% set parameters
%get information out of study_design
study_design=evalin('base','study_design');
contrast_def=evalin('base','contrast_def');
runs=study_design.number_sessions; 
nr_subj=study_design.number_subjects;
exStats = study_design.exist_stats;
ex4D = study_design.exist_4D;
if exStats == 1
    load(study_design.subject_list);
    stats_dir=study_design.stats_directory;
    stats_path=study_design.stats_path;
end;
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
    end;
elseif ex4D == 1
    nr_cond = contrast_def.number_conditions;
    if nr_cond > 1
        conditions = contrast_def.conditions;
    end;
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
    compl = [roi_dir f roi_ful];

    % reslice ROI
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
    r_roi=dir(sprintf('r%s*',roi_name));
    if length(r_roi)==2
        if ~strcmp(roi_dir,dir_results)
        compl1 = [roi_dir f 'r' roi_name '.img'];
        movefile(compl1,dir_results,'f');
        compl2 = [roi_dir f 'r' roi_name '.hdr'];
        movefile(compl2,dir_results,'f');
        end;
        cd(dir_results);    
        r_roi = load_nii(sprintf('r%s.img',roi_name));
        r_roi_ind = r_roi.img==1;
    else
        compl = [roi_dir f 'r' roi_name '.nii'];
        if ~strcmp(roi_dir,dir_results)
        movefile(compl,dir_results,'f');
        end;
        cd(dir_results);    
        r_roi = load_nii(sprintf('r%s.nii',roi_name));
        r_roi_ind = r_roi.img==1;        
    end;
    cd(dir_results);
    clear r_roi roi_ful;
end;    
%% load 4D images whole-brain, without ROI
if exStats == 1
    disp('...loads image dimensions..');
    stats_temp =sprintf(stats_dir,1);
    temp_img = [stats_path f id{1} f stats_temp f con];
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
        clear img;
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
                    img2 = load_nii(sprintf('4D_split2_par%d_%d.nii',ind_para,i));
                    eval(sprintf('img_%d_par%d_split1 = img1.img;',i,ind_para));
                    eval(sprintf('img_%d_par%d_split2 = img2.img;',i,ind_para)); 
                end;
            end;
        end;     
        clear img1 img2
    elseif two_cons == 1
        for i = 1:runs    
            if runs == 1
                i = single_run;
            end;
            img1 = load_nii(sprintf('4D_%s_%d.nii',con1,i));
            img2 = load_nii(sprintf('4D_%s_%d.nii',con2,i));
            eval(sprintf('img_%d_con1 = img1.img;',i));
            eval(sprintf('img_%d_con2 = img2.img;',i));     
%             if nr_para > 0
%                 for ind_para = 1:nr_para
%                     img = load_nii(sprintf('4D_%s_par%d_%d.nii',con1,ind_para,i));
%                     eval(sprintf('img_con1_par%d_%d = img.img;',ind_para,i));
%                     img = load_nii(sprintf('4D_%s_par%d_%d.nii',con2,ind_para,i));
%                     eval(sprintf('img_con2_par%d_%d = img.img;',ind_para,i));                    
%                 end;
%             end;
        end; 
        clear img1 img2 img
    end;
elseif ex4D == 1
    for ind_cond = 1:nr_cond
        for ind_run = 1:runs
            fprintf('...load 4D image for run %d...\n',ind_run);
            if nr_cond > 1
                fprintf('...condition %s...\n',conditions{ind_cond,1})
            end;
            
            if nr_cond == 1
                file = [dir_results f '4D_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file;',ind_run));
                cd(dir_results);
                temp_img = load_nii('temp_img.nii');
                eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%d = FourD%d.img;',ind_run,ind_run));

                %image dimensions 
                eval(sprintf('dims = size(FourD%d(:,:,:,1));',ind_run));
                x = dims(1);
                y = dims(2);
                z = dims(3);
            else
                file = [dir_results f '4D_' conditions{ind_cond,1} '_' num2str(ind_run) '.nii'];
                eval(sprintf('FourD%d = file;',ind_run));
                cd(dir_results);
                temp_img = load_nii('temp_img.nii');
                eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
                eval(sprintf('FourD%s%d = FourD%d.img;',conditions{ind_cond,1},ind_run,ind_run));

                %image dimensions
                eval(sprintf('dims = size(FourD%s%d(:,:,:,1));',conditions{ind_cond,1},ind_run));
                x = dims(1);
                y = dims(2);
                z = dims(3);                           
            end;
        end;
    end;
end;
nr_vox = x*y*z;

%% create correlation maps 
cd (dir_results); 
    % for naming of outputs
    if roi == 1
        str = 'ROI';
    else
        str = '';
    end;

    
    if roi == 1
        disp('...calculate correlation of mean ROI activity...')
        if exStats == 1
            if split == 0 && two_cons == 0
                for i = 1:runs
                    eval(sprintf('img_%d_ROI = NaN(nr_subj,1);',i));
                    for i_subj = 1:nr_subj
                        img_subj = [];
                        eval(sprintf('img_subj = img_%d(:,:,:,%d);',i,i_subj));
                        eval(sprintf('img_%d_ROI(i_subj,1) = mean(img_subj(r_roi_ind);',i));
                    end;
                    if nr_para > 0
                        for ind_para = 1:nr_para
                            eval(sprintf('img_par%d_%d_ROI = NaN(nr_subj,1);',ind_para,i));
                            for i_subj = 1:nr_subj
                                img_subj = [];
                                eval(sprintf('img_subj = img_par%d_%d(:,:,:,%d);',ind_para,i,i_subj));
                                eval(sprintf('img_par%d_%d_ROI(i_subj,1) = mean(img_subj(r_roi_ind);',ind_para,i));
                            end;
                        end;
                    end;
                    clear img_subj;
                end;
            elseif split == 1
                for i = 1:runs
                    if runs == 1
                        i = single_run;
                    end;
                    eval(sprintf('img_%d_ROI_split1 = NaN(nr_subj,1);',i));
                    eval(sprintf('img_%d_ROI_split2 = NaN(nr_subj,1);',i));
                    for i_subj = 1:nr_subj
                        img_subj = [];
                        eval(sprintf('img_subj = img_%d_split1(:,:,:,%d);',i,i_subj));
                        eval(sprintf('img_%d_ROI_split1(i_subj,1) = mean(img_subj(r_roi_ind);',i));
                        img_subj = [];
                        eval(sprintf('img_subj = img_%d_split2(:,:,:,%d);',i,i_subj));
                        eval(sprintf('img_%d_ROI_split2(i_subj,1) = mean(img_subj(r_roi_ind);',i));
                    end;
                    if nr_para > 0
                        for ind_para = 1:nr_para
                            eval(sprintf('img_par%d_%d_ROI_split1 = NaN(nr_subj,1);',ind_para,i));
                            eval(sprintf('img_par%d_%d_ROI_split2 = NaN(nr_subj,1);',ind_para,i));
                            for i_subj = 1:nr_subj
                                img_subj = [];
                                eval(sprintf('img_subj = img_%d_par%d_split1(:,:,:,%d);',i,ind_para,i_subj));
                                eval(sprintf('img_par%d_%d_ROI_split1(i_subj,1) = mean(img_subj(r_roi_ind);',ind_para,i));
                                img_subj = [];
                                eval(sprintf('img_subj = img_%d_par%d_split2(:,:,:,%d);',i,ind_para,i_subj));
                                eval(sprintf('img_par%d_%d_ROI_split2(i_subj,1) = mean(img_subj(r_roi_ind);',ind_para,i));
                            end;
                        end;
                    end;
                end;
                
            elseif two_cons == 1
                for i = 1:runs
                    if runs == 1
                        i = single_run;
                    end;
                    eval(sprintf('img_%d_ROI_con1 = NaN(nr_subj,1);',i));
                    eval(sprintf('img_%d_ROI_con2 = NaN(nr_subj,1);',i));
                    for i_subj = 1:nr_subj
                        img_subj = [];
                        eval(sprintf('img_subj = img_%d_con1(:,:,:,%d);',i,i_subj));
                        eval(sprintf('img_%d_ROI_con1(i_subj,1) = mean(img_subj(r_roi_ind);',i));
                        img_subj = [];
                        eval(sprintf('img_subj = img_%d_con2(:,:,:,%d);',i,i_subj));
                        eval(sprintf('img_%d_ROI_con2(i_subj,1) = mean(img_subj(r_roi_ind);',i));
                    end;
                end;
            end;
            % initiate summary
            cols = {};
            summary = [];
            if split == 0 && two_cons == 0
                for i = 1:runs
                    for j = (runs-1):1
                        if i+j <= runs
                            if pear == 1 && spea == 1
                                eval(sprintf('first = img_%d_ROI;',i))
                                eval(sprintf('second = img_%d_ROI;',i+j))
                                r = corrcoef(first,second, 'rows', 'pairwise');
                                if isnan (r(1,2))
                                    r1(i_vox,1) = 0;
                                    zr1(i_vox,1) = 0;
                                else
                                    r1(i_vox,1) = r(1,2);
                                    zr1(i_vox,1) = atanh(r(1,2));
                                end;
                                r = corr(first,second, 'Type','Spearman', 'rows', 'pairwise');
                                r2(i_vox,1) = r;  
                                zr2(i_vox,1) = atanh(r);                            
                                
                                summary(1,end+1)=r1;
                                cols{1,end+1} = sprintf('Pearson_%d_%d',i,i+j);   
                                summary(1,end+1)=zr1;
                                cols{1,end+1} = sprintf('zPearson_%d_%d',i,i+j); 
                                summary(1,end+1)=r2;
                                cols{1,end+1} = sprintf('Spearman_%d_%d',i,i+j);   
                                summary(1,end+1)=zr2;
                                cols{1,end+1} = sprintf('zSpearman_%d_%d',i,i+j);
                                
                                clear first second
                                if nr_para > 0
                                    for ind_p = 1:nr_para
                                        eval(sprintf('first = img_par%d_%d_ROI;',ind_p,i))
                                        eval(sprintf('second = img_par%d_%d_ROI;',ind_p,i+j))
                                        r = corrcoef(first,second, 'rows', 'pairwise');
                                        if isnan (r(1,2))
                                            r1(i_vox,1) = 0;
                                            zr1(i_vox,1) = 0;
                                        else
                                            r1(i_vox,1) = r(1,2);
                                            zr1(i_vox,1) = atanh(r(1,2));
                                        end;
                                        r = corr(first,second, 'Type','Spearman', 'rows', 'pairwise');
                                        r2(i_vox,1) = r;  
                                        zr2(i_vox,1) = atanh(r);                            

                                        summary(1,end+1)=r1;
                                        cols{1,end+1} = sprintf('Pearson_%d_%d_par%d',i,i+j,ind_p);   
                                        summary(1,end+1)=zr1;
                                        cols{1,end+1} = sprintf('zPearson_%d_%d_par%d',i,i+j,ind_p); 
                                        summary(1,end+1)=r2;
                                        cols{1,end+1} = sprintf('Spearman_%d_%d_par%d',i,i+j,ind_p);   
                                        summary(1,end+1)=zr2;
                                        cols{1,end+1} = sprintf('zSpearman_%d_%d_par%d',i,i+j,ind_p);

                                        clear first second                                        
                                    end;
                                end;
                                
                                
                            elseif pear == 1 && spea == 0
                                eval(sprintf('first = img_%d_ROI;',i))
                                eval(sprintf('second = img_%d_ROI;',i+j))
                                r = corrcoef(first,second, 'rows', 'pairwise');
                                if isnan (r(1,2))
                                    r1(i_vox,1) = 0;
                                    zr1(i_vox,1) = 0;
                                else
                                    r1(i_vox,1) = r(1,2);
                                    zr1(i_vox,1) = atanh(r(1,2));
                                end;
                                summary(1,end+1)=r1;
                                cols{1,end+1} = sprintf('Pearson_%d_%d',i,i+j);   
                                summary(1,end+1)=zr1;
                                cols{1,end+1} = sprintf('zPearson_%d_%d',i,i+j); 
                                
                                clear first second 
                                if nr_para > 0
                                    for ind_p = 1:nr_para
                                        eval(sprintf('first = img_par%d_%d_ROI;',ind_p,i))
                                        eval(sprintf('second = img_par%d_%d_ROI;',ind_p,i+j))
                                        r = corrcoef(first,second, 'rows', 'pairwise');
                                        if isnan (r(1,2))
                                            r1(i_vox,1) = 0;
                                            zr1(i_vox,1) = 0;
                                        else
                                            r1(i_vox,1) = r(1,2);
                                            zr1(i_vox,1) = atanh(r(1,2));
                                        end;
                                        summary(1,end+1)=r1;
                                        cols{1,end+1} = sprintf('Pearson_%d_%d_par%d',i,i+j,ind_p);   
                                        summary(1,end+1)=zr1;
                                        cols{1,end+1} = sprintf('zPearson_%d_%d_par%d',i,i+j,ind_p); 

                                        clear first second 
                                    end;
                                end;
                            elseif pear == 0 && spea == 1
                                eval(sprintf('first = img_%d_ROI;',i))
                                eval(sprintf('second = img_%d_ROI;',i+j))
                                r = corr(first,second, 'Type','Spearman', 'rows', 'pairwise');
                                r2(i_vox,1) = r;  
                                zr2(i_vox,1) = atanh(r); 
                                summary(1,end+1)=r2;
                                cols{1,end+1} = sprintf('Spearman_%d_%d',i,i+j);   
                                summary(1,end+1)=zr2;
                                cols{1,end+1} = sprintf('zSpearman_%d_%d',i,i+j);
                                
                                clear first second
                                if nr_para > 0
                                    for ind_p = 1:nr_para
                                        eval(sprintf('first = img_par%d_%d_ROI;',ind_p,i))
                                        eval(sprintf('second = img_par%d_%d_ROI;',ind_p,i+j))
                                        r = corr(first,second, 'Type','Spearman', 'rows', 'pairwise');
                                        r2(i_vox,1) = r;  
                                        zr2(i_vox,1) = atanh(r); 
                                        summary(1,end+1)=r2;
                                        cols{1,end+1} = sprintf('Spearman_%d_%d_par%d',i,i+j,ind_p);   
                                        summary(1,end+1)=zr2;
                                        cols{1,end+1} = sprintf('zSpearman_%d_%d_par%d',i,i+j,ind_p);
                                        clear first second 
                                    end;
                                end;
                            end;    
                        end;
                    end;
                end;
            elseif two_cons == 1
            elseif split == 1
            end;
        elseif ex4D == 1
            for ind_cond = 1:nr_cond
                for i = 1:runs                      
                    if nr_cond == 1
                        eval(sprintf('img_%d_ROI = NaN(nr_subj,1);',i));
                        for i_subj = 1:nr_subj
                            img_subj = [];
                            eval(sprintf('img_subj = FourD%d(:,:,:,%d);',i,i_subj));
                            eval(sprintf('img_%d_ROI(i_subj,1) = mean(img_subj(r_roi_ind);',i));
                            clear img_subj
                        end;
                    else
                        if runs == 1
                            i = single_run;
                        end;
                        eval(sprintf('img_%d_ROI_%s = NaN(nr_subj,1);',i,conditions{ind_cond,1}));
                        for i_subj = 1:nr_subj
                            img_subj = [];
                            eval(sprintf('img_subj = FourD%s%d(:,:,:,%d);',conditions{ind_cond,1},i,i_subj));
                            eval(sprintf('img_%d_ROI_%s(i_subj,1) = mean(img_subj(r_roi_ind);',i,conditions{ind_cond,1}));
                            clear img_subj
                        end;

                    end;
                end;
            end;            
        end;
    end;

if pear == 1 || spea == 1
% initiate summary
cols = {};
summary = [];
if exStats == 1
if split == 0 && two_cons == 0 && runs > 1    
    count_comp = 0;
    % create data matrix with data for all voxels in cols, subjects in rows
    % and runs in third dimension
    data = zeros(nr_subj,nr_vox,runs);
    for ind_runs = 1:runs
        for ind_subj = 1:nr_subj
            eval(sprintf('temp = img_%d(:,:,:,ind_subj);',ind_runs));
            data(ind_subj,:,ind_runs)=temp(:);
        end;
    end;
    clear temp;
    
    %correlations for each comparison
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                count_comp = count_comp +1;
                fprintf('...creates correlation maps for session %d and session %d...\n',i_run,i_run+i_sec);
                %load 4D images
                eval(sprintf('one=data(:,:,%d);',i_run));
                eval(sprintf('second=data(:,:,%d);',i_run+i_sec));
                
                r1 = zeros(nr_vox,1);
                zr1 = zeros(nr_vox,1);
                r2 = zeros(nr_vox,1);
                zr2 = zeros(nr_vox,1);
                for i_vox = 1:nr_vox
                    if pear == 1 
                        r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
                            if isnan (r(1,2))
                                r1(i_vox,1) = 0;
                                zr1(i_vox,1) = 0;
                            else
                                r1(i_vox,1) = r(1,2);
                                zr1(i_vox,1) = atanh(r(1,2));
                            end;
                    end
                    if spea == 1
                            r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                            r2(i_vox,1) = r;  
                            zr2(i_vox,1) = atanh(r);  
                    end;
                    
                end;
                clear one second
                
                % create correlation matrices for image 
                if pear == 1
                r_vec_pear_1_2 = reshape(r1,x,y,z);
                z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
                clear r1 zr1
                end;
                
                if spea == 1
                r_vec_spea_1_2 = reshape(r2,x,y,z);
                z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
                clear r2 zr2 
                end;
             

            % save correlation maps
            if pear == 1
            target_img = temp_img;

            file = sprintf('CorrMaps%s_pear_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_pear_1_2;
                else
                    r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_pear_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_pear_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_pear_1_2;
                else
                    z_r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_pear_1_2;
                end;            
            save_nii(target_img,target_img.fileprefix);      
            end;
            
            if spea == 1
            target_img = temp_img;
            file = sprintf('CorrMaps%s_spea_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_spea_1_2;
                else
                    r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_spea_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_spea_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_spea_1_2;
                else
                    z_r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_spea_1_2;
                end;               
            save_nii(target_img,target_img.fileprefix);                      
            end;
            
            %compute z mean and inverse Fisher's transformation
            if pear == 1
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear%s_%d_%d',str,i_run,i_run+i_sec);
            end;
            if spea == 1
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea%s_%d_%d',str,i_run,i_run+i_sec);    
            end;
            end; 
        end;
    end;
    clear r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2 data

    % calculate average correlation for all comparisons
    disp('...creates average correlation maps ...');   
    avg_corr_4D_pear = zeros(nr_vox,count_comp);
    avg_corr_4D_spea = zeros(nr_vox,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                ind_comp = ind_comp+1;
                if pear == 1
                map = load_nii(sprintf('z_CorrMaps%s_pear_%d_%d.nii',str,i_run,i_run+i_sec));
                avg_corr_4D_pear(:,ind_comp) = map.img(:);
                clear map
                end;
                if spea == 1
                map = load_nii(sprintf('z_CorrMaps%s_spea_%d_%d.nii',str,i_run,i_run+i_sec));
                avg_corr_4D_spea(:,ind_comp) = map.img(:);
                clear map
                end;                
            end; 
        end;
    end;
    clear map
    
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                

    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear.nii',str);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea.nii',str);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
clear avg_spea avg_pear avg_corr_4D_spea avg_corr_4D_pear z_r_vec_spea_1_2 

if nr_para > 0    
     count_comp = 0;
     for i_par = 1:nr_para
        data = zeros(nr_subj,nr_vox,runs);
        for ind_runs = 1:runs
            for ind_subj = 1:nr_subj
                eval(sprintf('temp = img_par%d_%d(:,:,:,ind_subj);',i_par,ind_runs));
                data(ind_subj,:,ind_runs)=temp(:);
            end;
        end;
        clear temp
    
    %correlations for each comparison
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                count_comp = count_comp +1;
                fprintf('...creates correlation maps for parametric contrast in session %d and session %d...\n',i_run,i_run+i_sec);
                %load 4D images
                eval(sprintf('one=data(:,:,%d);',i_run));
                eval(sprintf('second=data(:,:,%d);',i_run+i_sec));
                
                r1 = zeros(nr_vox,1);
                zr1 = zeros(nr_vox,1);
                r2 = zeros(nr_vox,1);
                zr2 = zeros(nr_vox,1);
                for i_vox = 1:nr_vox
                    if pear == 1 
                        r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
                            if isnan (r(1,2))
                                r1(i_vox,1) = 0;
                                zr1(i_vox,1) = 0;
                            else
                                r1(i_vox,1) = r(1,2);
                                zr1(i_vox,1) = atanh(r(1,2));
                            end;
                    end
                    if spea == 1
                            r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                            r2(i_vox,1) = r;  
                            zr2(i_vox,1) = atanh(r);  
                    end;
                    
                end;
                clear one second
                
                % create correlation matrices for image 
                if pear == 1
                    r_vec_pear_1_2 = reshape(r1,x,y,z);
                    z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
                    clear r1 zr1
                end;
                if spea == 1
                    r_vec_spea_1_2 = reshape(r2,x,y,z);
                    z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
                    clear r2 zr2;
                end;
             

            % save correlation maps
            if pear == 1
            target_img = temp_img;

            file = sprintf('CorrMaps%s_pear_par%d_%d_%d.nii',str,i_par,i_run,i_run+i_sec);
            target_img.fileprefix = file;
            if roi == 0
                target_img.img = r_vec_pear_1_2;
            else
                r_vec_pear_1_2(~r_roi_ind)=0;
                target_img.img = r_vec_pear_1_2;
            end;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_pear_par%d_%d_%d.nii',str,i_par,i_run,i_run+i_sec);
            target_img.fileprefix = file;
            if roi == 0
                target_img.img = z_r_vec_pear_1_2;
            else
                z_r_vec_pear_1_2(~r_roi_ind)=0;
                target_img.img = z_r_vec_pear_1_2;
            end;            
            save_nii(target_img,target_img.fileprefix);      
            end;
            
            if spea == 1
            target_img = temp_img;
            file = sprintf('CorrMaps%s_spea_par%d_%d_%d.nii',str,i_par,i_run,i_run+i_sec);
            target_img.fileprefix = file;
            if roi == 0
                target_img.img = r_vec_spea_1_2;
            else
                r_vec_spea_1_2(~r_roi_ind)=0;
                target_img.img = r_vec_spea_1_2;
            end;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_spea_par%d_%d_%d.nii',str,i_par,i_run,i_run+i_sec);
            target_img.fileprefix = file;
            if roi == 0
                target_img.img = z_r_vec_spea_1_2;
            else
                z_r_vec_spea_1_2(~r_roi_ind)=0;
                target_img.img = z_r_vec_spea_1_2;
            end;               
            save_nii(target_img,target_img.fileprefix);                      
            end;
            
            %compute z mean and inverse Fisher's transformation
            if pear == 1
            mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_pear%s_par%d_%d_%d',str,i_par,i_run,i_run+i_sec);
            end;
            if spea == 1
            mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
            mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
            summary(1,end+1)=mean_r;
            cols{1,end+1} = sprintf('mean_spea%s_par%d_%d_%d',str,i_par,i_run,i_run+i_sec);    
            end;
            clear z_r_vec_spea_1_2 r_vec_spea_1_2 r_vec_pear_1_2 z_r_vec_pear_1_2
            end; 
        end;
    end;
    clear data 
    % calculate average correlation for all comparisons
    disp('...creates average correlation maps for parametric contrast...');   
    avg_corr_4D_pear = zeros(nr_vox,count_comp);
    avg_corr_4D_spea = zeros(nr_vox,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                ind_comp = ind_comp+1;
                if pear == 1
                map = load_nii(sprintf('z_CorrMaps%s_pear_par%d_%d_%d.nii',str,i_par,i_run,i_run+i_sec));
                avg_corr_4D_pear(:,ind_comp) = map.img(:);
                end;
                if spea == 1
                map = load_nii(sprintf('z_CorrMaps%s_spea_par%d_%d_%d.nii',str,i_par,i_run,i_run+i_sec));
                avg_corr_4D_spea(:,ind_comp) = map.img(:);
                end;                
            end; 
        end;
    end;
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                

    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear_par%d.nii',str,i_par);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea_par%d.nii',str,i_par);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
clear avg_spea avg_pear avg_corr_4D_spea avg_corr_4D_pear z_r_vec_spea_1_2 
     end;
end;
results_corr=dataset({summary(1,:),cols{:}});
assignin('base','results_corr',results_corr);
    
cd(dir_results);             
eval(sprintf('save results_corr%s.mat results_corr',str)); 
%% based on split half
elseif split == 1
    count_comp = 0;
  for i = 1:runs
        if runs == 1
            i = single_run;
        end;
        % create data matrix with data for all voxels in cols, subjects in rows
        % and 2 splits in third dimension
        data1 = zeros(nr_subj,nr_vox);
        data2 = zeros(nr_subj,nr_vox);
            for ind_subj = 1:nr_subj
                eval(sprintf('temp = img_%d_split1(:,:,:,ind_subj);',i));
                data1(ind_subj,:)=temp(:);
                eval(sprintf('temp = img_%d_split2(:,:,:,ind_subj);',i));
                data2(ind_subj,:)=temp(:);                
            end;
        clear temp
        eval(sprintf('clear img_%d_split1 img_%d_split2;',i,i));
    
    %correlations for each split
    count_comp = count_comp +1;
    fprintf('...creates correlation maps for splitted session %d ...\n',i);
    %load 4D images
%     one=data(:,:,1);
%     second=data(:,:,2);

    r1 = zeros(nr_vox,1);
    zr1 = zeros(nr_vox,1);
    r2 = zeros(nr_vox,1);
    zr2 = zeros(nr_vox,1);
    for i_vox = 1:nr_vox
        if pear == 1 
            r = corrcoef(data1(:,i_vox),data2(:,i_vox), 'rows', 'pairwise');
                if isnan (r(1,2))
                    r1(i_vox,1) = 0;
                    zr1(i_vox,1) = 0;
                else
                    r1(i_vox,1) = r(1,2);
                    zr1(i_vox,1) = atanh(r(1,2));
                end;
        end
        if spea == 1
                r = corr(data1(:,i_vox),data2(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                r2(i_vox,1) = r;  
                zr2(i_vox,1) = atanh(r);  
        end;
    end;
%     clear one second
                
    % create correlation matrices for image 
    r_vec_pear_1_2 = reshape(r1,x,y,z);
    r_vec_spea_1_2 = reshape(r2,x,y,z);
    z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
    z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
    clear r1 zr1 r2 zr2;
             
    % save correlation maps
    if pear == 1
        target_img = temp_img;
        file = sprintf('CorrMaps%s_pear_%d_split.nii',str,i);
        target_img.fileprefix = file;
        if roi == 0
            target_img.img = r_vec_pear_1_2;
        else
            r_vec_pear_1_2(~r_roi_ind)=0;
            target_img.img = r_vec_pear_1_2;
        end;
        save_nii(target_img,target_img.fileprefix);  

        target_img = temp_img;
        file = sprintf('z_CorrMaps%s_pear_%d_split.nii',str,i);
        target_img.fileprefix = file;
        if roi == 0
            target_img.img = z_r_vec_pear_1_2;
        else
            z_r_vec_pear_1_2(~r_roi_ind)=0;
            target_img.img = z_r_vec_pear_1_2;
        end;            
        save_nii(target_img,target_img.fileprefix);      
    end;

    if spea == 1
    target_img = temp_img;
    file = sprintf('CorrMaps%s_spea_%d_split.nii',str,i);
    target_img.fileprefix = file;
    if roi == 0
        target_img.img = r_vec_spea_1_2;
    else
        r_vec_spea_1_2(~r_roi_ind)=0;
        target_img.img = r_vec_spea_1_2;
    end;
    save_nii(target_img,target_img.fileprefix);

    target_img = temp_img;
    file = sprintf('z_CorrMaps%s_spea_%d_split.nii',str,i);
    target_img.fileprefix = file;
    if roi == 0
        target_img.img = z_r_vec_spea_1_2;
    else
        z_r_vec_spea_1_2(~r_roi_ind)=0;
        target_img.img = z_r_vec_spea_1_2;
    end;               
    save_nii(target_img,target_img.fileprefix);                      
    end;

    %compute z mean and inverse Fisher's transformation
    if pear == 1
        mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_pear%s_%d_split',str,i);
    end;
    if spea == 1
        mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_spea%s_%d_split',str,i);    
    end;
    clear data1 data 2
  end; 
clear z_r_vec_spea_1_2 r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2

    % calculate average correlation for all comparisons
    disp('...creates average correlation maps over all splits ...');   
    avg_corr_4D_pear = zeros(nr_vox,count_comp);
    avg_corr_4D_spea = zeros(nr_vox,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        ind_comp = ind_comp+1;
        if pear == 1
        map = load_nii(sprintf('z_CorrMaps%s_pear_%d_split.nii',str,i));
        avg_corr_4D_pear(:,ind_comp) = map.img(:);
        clear map
        end;
        if spea == 1
        map = load_nii(sprintf('z_CorrMaps%s_spea_%d_split.nii',str,i));
        avg_corr_4D_spea(:,ind_comp) = map.img(:);
        clear map
        end;                
    end;
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                
    clear avg_corr_4D_pear avg_corr_4D_spea
    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear_split.nii',str);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea_split.nii',str);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
    clear avg_spea avg_pear data
%calculation for splitted parametric contrasts
if nr_para > 0    
     count_comp = 0;
     for i_par = 1:nr_para
        for i = 1:runs
            if runs == 1
                i = single_run;
            end;
            % create data matrix with data for all voxels in cols, subjects in rows
            % and 2 splits in third dimension
            data = zeros(nr_subj,nr_vox,2);
            for ind_split = 1:2
                for ind_subj = 1:nr_subj
                    eval(sprintf('temp = img_%d_par%d_split%d(:,:,:,ind_subj);',i,i_par,ind_split));
                    data(ind_subj,:,ind_split)=temp(:);
                end;
            end;
            clear temp
            eval(sprintf('clear img_%d_par%d_split1 img_%d_par%d_split2;',i,i_par,i,i_par));
            %correlations for each split
            count_comp = count_comp +1;
            fprintf('...creates correlation maps for splitted session %d ...\n',i);
            %load 4D images
             one=data(:,:,1);
             second=data(:,:,2);

            r1 = zeros(nr_vox,1);
            zr1 = zeros(nr_vox,1);
            r2 = zeros(nr_vox,1);
            zr2 = zeros(nr_vox,1);
            for i_vox = 1:nr_vox
                if pear == 1 
                    r = corrcoef(one(:,i_vox),second(:,i_vox), 'rows', 'pairwise');
                        if isnan (r(1,2))
                            r1(i_vox,1) = 0;
                            zr1(i_vox,1) = 0;
                        else
                            r1(i_vox,1) = r(1,2);
                            zr1(i_vox,1) = atanh(r(1,2));
                        end;
                end
                if spea == 1
                        r = corr(one(:,i_vox),second(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                        r2(i_vox,1) = r;  
                        zr2(i_vox,1) = atanh(r);  
                end;
            end;
             clear one second
                
            % create correlation matrices for image 
            r_vec_pear_1_2 = reshape(r1,x,y,z);
            r_vec_spea_1_2 = reshape(r2,x,y,z);
            z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
            z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
            clear r1 zr1 r2 zr2;
             
            % save correlation maps
            if pear == 1
                target_img = temp_img;
                file = sprintf('CorrMaps%s_pear_%d_par%d_split.nii',str,i,i_par);
                target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_pear_1_2;
                else
                    r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_pear_1_2;
                end;
                save_nii(target_img,target_img.fileprefix);  

                target_img = temp_img;
                file = sprintf('z_CorrMaps%s_pear_%d_par%d_split.nii',str,i,i_par);
                target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_pear_1_2;
                else
                    z_r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_pear_1_2;
                end;            
                save_nii(target_img,target_img.fileprefix);      
            end;

            if spea == 1
            target_img = temp_img;
            file = sprintf('CorrMaps%s_spea_%d_par%d_split.nii',str,i,i_par);
            target_img.fileprefix = file;
            if roi == 0
                target_img.img = r_vec_spea_1_2;
            else
                r_vec_spea_1_2(~r_roi_ind)=0;
                target_img.img = r_vec_spea_1_2;
            end;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_spea_%d_par%d_split.nii',str,i,i_par);
            target_img.fileprefix = file;
            if roi == 0
                target_img.img = z_r_vec_spea_1_2;
            else
                z_r_vec_spea_1_2(~r_roi_ind)=0;
                target_img.img = z_r_vec_spea_1_2;
            end;               
            save_nii(target_img,target_img.fileprefix);                      
            end;

            %compute z mean and inverse Fisher's transformation
            if pear == 1
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear%s_%d_par%d_split',str,i,i_par);
            end;
            if spea == 1
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea%s_%d_par%d_split',str,i,i_par);    
            end;
        end; 
        clear z_r_vec_spea_1_2 r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2

    % calculate average correlation for all comparisons
    disp('...creates average correlation maps over all splits of parametric contrast...');   
    avg_corr_4D_pear = zeros(nr_vox,count_comp);
    avg_corr_4D_spea = zeros(nr_vox,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        ind_comp = ind_comp+1;
        if pear == 1
        map = load_nii(sprintf('z_CorrMaps%s_pear_%d_par%d_split.nii',str,i,i_par));
        avg_corr_4D_pear(:,ind_comp) = map.img(:);
        clear map
        end;
        if spea == 1
        map = load_nii(sprintf('z_CorrMaps%s_spea_%d_par%d_split.nii',str,i,i_par));
        avg_corr_4D_spea(:,ind_comp) = map.img(:);
        clear map
        end;                
    end;
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                
    clear avg_corr_4D_pear avg_corr_4D_spea
    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear_par%d_split.nii',str,i_par);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea_par%d_split.nii',str,i_par);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
         
  
     end;
end;  
results_corr=dataset({summary(1,:),cols{:}});
assignin('base','results_corr',results_corr);
    
cd(dir_results);             
eval(sprintf('save results_corr%s_split.mat results_corr',str)); 
%% two contrasts out of one statistic    
elseif two_cons == 1
     % compare contrasts within sessions
     for i_run = 1:runs
            if runs == 1
                i_run = single_run;
            end;
        % create data matrix with data for all voxels in cols, subjects in rows
        % and 2 splits in third dimension
        data = zeros(nr_subj,nr_vox,2);
        for ind_split = 1:2
            for ind_subj = 1:nr_subj
                eval(sprintf('temp = img_%d_con%d(:,:,:,ind_subj);',i_run,ind_split));
                data(ind_subj,:,ind_split)=temp(:);
            end;
        end;
        clear temp
    
    %correlations for each split
    fprintf('...creates correlation maps for two contrasts in session %d ...\n',i_run);
    %load 4D images
    one=data(:,:,1);
    second=data(:,:,2);

    r1 = zeros(nr_vox,1);
    zr1 = zeros(nr_vox,1);
    r2 = zeros(nr_vox,1);
    zr2 = zeros(nr_vox,1);
    for i_vox = 1:nr_vox
        if pear == 1 
            r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
                if isnan (r(1,2))
                    r1(i_vox,1) = 0;
                    zr1(i_vox,1) = 0;
                else
                    r1(i_vox,1) = r(1,2);
                    zr1(i_vox,1) = atanh(r(1,2));
                end;
        end
        if spea == 1
                r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                r2(i_vox,1) = r;  
                zr2(i_vox,1) = atanh(r);  
        end;
    end;
    clear one second
                
    % create correlation matrices for image 
    r_vec_pear_1_2 = reshape(r1,x,y,z);
    r_vec_spea_1_2 = reshape(r2,x,y,z);
    z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
    z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
    clear r1 zr1 r2 zr2;
             
    % save correlation maps
    if pear == 1
        target_img = temp_img;
        file = sprintf('CorrMaps%s_pear_%d_between_contrasts.nii',str,i_run);
        target_img.fileprefix = file;
        if roi == 0
            target_img.img = r_vec_pear_1_2;
        else
            r_vec_pear_1_2(~r_roi_ind)=0;
            target_img.img = r_vec_pear_1_2;
        end;
        save_nii(target_img,target_img.fileprefix);  

        target_img = temp_img;
        file = sprintf('z_CorrMaps%s_pear_%d_between_contrasts.nii',str,i_run);
        target_img.fileprefix = file;
        if roi == 0
            target_img.img = z_r_vec_pear_1_2;
        else
            z_r_vec_pear_1_2(~r_roi_ind)=0;
            target_img.img = z_r_vec_pear_1_2;
        end;            
        save_nii(target_img,target_img.fileprefix);      
    end;

    if spea == 1
    target_img = temp_img;
    file = sprintf('CorrMaps%s_spea_%d_between_contrasts.nii',str,i_run);
    target_img.fileprefix = file;
    if roi == 0
        target_img.img = r_vec_spea_1_2;
    else
        r_vec_spea_1_2(~r_roi_ind)=0;
        target_img.img = r_vec_spea_1_2;
    end;
    save_nii(target_img,target_img.fileprefix);

    target_img = temp_img;
    file = sprintf('z_CorrMaps%s_spea_%d_between_contrasts.nii',str,i_run);
    target_img.fileprefix = file;
    if roi == 0
        target_img.img = z_r_vec_spea_1_2;
    else
        z_r_vec_spea_1_2(~r_roi_ind)=0;
        target_img.img = z_r_vec_spea_1_2;
    end;               
    save_nii(target_img,target_img.fileprefix);                      
    end;

    %compute z mean and inverse Fisher's transformation
    if pear == 1
        mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_pear%s_%d_between_contrasts',str,i_run);
    end;
    if spea == 1
        mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_spea%s_%d_between_contrasts',str,i_run);    
    end;
clear z_r_vec_spea_1_2 r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2
     end;
     
    % calculate average correlation for all comparisons
    disp('...creates average correlation maps over all between contrast correlations ...');   
    avg_corr_4D_pear = zeros(nr_vox,runs);
    avg_corr_4D_spea = zeros(nr_vox,runs);
    ind_comp = 0;
    for i_run = 1:runs
        ind_comp = ind_comp+1;
        if pear == 1
        map = load_nii(sprintf('z_CorrMaps%s_pear_%d_between_contrasts.nii',str,i));
        avg_corr_4D_pear(:,ind_comp) = map.img(:);
        clear map
        end;
        if spea == 1
        map = load_nii(sprintf('z_CorrMaps%s_spea_%d_between_contrasts.nii',str,i));
        avg_corr_4D_spea(:,ind_comp) = map.img(:);
        clear map
        end;                
    end;
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                
    clear avg_corr_4D_pear avg_corr_4D_spea
    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear_between_contrasts.nii',str);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea_between_contrasts.nii',str);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
%calculation for parametric regressors of both contrasts
% if nr_para > 0    
%      count_comp = 0;
%      for i_par = 1:nr_para
%         for i = 1:runs
%             if runs == 1
%                 i = single_run;
%             end;
%             % create data matrix with data for all voxels in cols, subjects in rows
%             % and 2 splits in third dimension
%             data = zeros(nr_subj,nr_vox,2);
%             for ind_split = 1:2
%                 for ind_subj = 1:nr_subj
%                     eval(sprintf('temp = img_%d_par%d_con%d(:,:,:,ind_subj);',i,i_par,ind_split));
%                     data(ind_subj,:,ind_split)=temp(:);
%                 end;
%             end;
%             clear temp
%     
%             %correlations for each split
%             count_comp = count_comp +1;
%             fprintf('...creates correlation maps for parametric regressors of both contrasts in session %d ...\n',i);
%             %load 4D images
%             one=data(:,:,1);
%             second=data(:,:,2);
% 
%             r1 = zeros(nr_vox,1);
%             zr1 = zeros(nr_vox,1);
%             r2 = zeros(nr_vox,1);
%             zr2 = zeros(nr_vox,1);
%             for i_vox = 1:nr_vox
%                 if pear == 1 
%                     r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
%                         if isnan (r(1,2))
%                             r1(i_vox,1) = 0;
%                             zr1(i_vox,1) = 0;
%                         else
%                             r1(i_vox,1) = r(1,2);
%                             zr1(i_vox,1) = atanh(r(1,2));
%                         end;
%                 end
%                 if spea == 1
%                         r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
%                         r2(i_vox,1) = r;  
%                         zr2(i_vox,1) = atanh(r);  
%                 end;
%             end;
%             clear one second
%                 
%             % create correlation matrices for image 
%             r_vec_pear_1_2 = reshape(r1,x,y,z);
%             r_vec_spea_1_2 = reshape(r2,x,y,z);
%             z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
%             z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
%             clear r1 zr1 r2 zr2;
%              
%             % save correlation maps
%             if pear == 1
%                 target_img = temp_img;
%                 file = sprintf('CorrMaps%s_pear_%d_par%d_between_contrasts.nii',str,i,i_par);
%                 target_img.fileprefix = file;
%                 if roi == 0
%                     target_img.img = r_vec_pear_1_2;
%                 else
%                     r_vec_pear_1_2(~r_roi_ind)=0;
%                     target_img.img = r_vec_pear_1_2;
%                 end;
%                 save_nii(target_img,target_img.fileprefix);  
% 
%                 target_img = temp_img;
%                 file = sprintf('z_CorrMaps%s_pear_%d_par%d_between_contrasts.nii',str,i,i_par);
%                 target_img.fileprefix = file;
%                 if roi == 0
%                     target_img.img = z_r_vec_pear_1_2;
%                 else
%                     z_r_vec_pear_1_2(~r_roi_ind)=0;
%                     target_img.img = z_r_vec_pear_1_2;
%                 end;            
%                 save_nii(target_img,target_img.fileprefix);      
%             end;
% 
%             if spea == 1
%             target_img = temp_img;
%             file = sprintf('CorrMaps%s_spea_%d_par%d_between_contrasts.nii',str,i,i_par);
%             target_img.fileprefix = file;
%             if roi == 0
%                 target_img.img = r_vec_spea_1_2;
%             else
%                 r_vec_spea_1_2(~r_roi_ind)=0;
%                 target_img.img = r_vec_spea_1_2;
%             end;
%             save_nii(target_img,target_img.fileprefix);
% 
%             target_img = temp_img;
%             file = sprintf('z_CorrMaps%s_spea_%d_par%d_between_contrasts.nii',str,i,i_par);
%             target_img.fileprefix = file;
%             if roi == 0
%                 target_img.img = z_r_vec_spea_1_2;
%             else
%                 z_r_vec_spea_1_2(~r_roi_ind)=0;
%                 target_img.img = z_r_vec_spea_1_2;
%             end;               
%             save_nii(target_img,target_img.fileprefix);                      
%             end;
% 
%             %compute z mean and inverse Fisher's transformation
%             if pear == 1
%                 mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
%                 mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%                 summary(1,end+1)=mean_r;
%                 cols{1,end+1} = sprintf('mean_pear%s_%d_par%d_between_contrasts',str,i,i_par);
%             end;
%             if spea == 1
%                 mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
%                 mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%                 summary(1,end+1)=mean_r;
%                 cols{1,end+1} = sprintf('mean_spea%s_%d_par%d_between_contrasts',str,i,i_par);    
%             end;
%         end; 
%         clear z_r_vec_spea_1_2 r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2
% 
%     % calculate average correlation for all comparisons
%     disp('...creates average correlation maps over all parametric regressors of both contrasts...');   
%     avg_corr_4D_pear = zeros(nr_vox,count_comp);
%     avg_corr_4D_spea = zeros(nr_vox,count_comp);
%     ind_comp = 0;
%     for i_run = 1:runs
%         ind_comp = ind_comp+1;
%         if pear == 1
%         map = load_nii(sprintf('z_CorrMaps%s_pear_%d_par%d_between_contrasts.nii',str,i,i_par));
%         avg_corr_4D_pear(:,ind_comp) = map.img(:);
%         clear map
%         end;
%         if spea == 1
%         map = load_nii(sprintf('z_CorrMaps%s_spea_%d_par%d_between_contrasts.nii',str,i,i_par));
%         avg_corr_4D_spea(:,ind_comp) = map.img(:);
%         clear map
%         end;                
%     end;
%     avg_pear = zeros(nr_vox,1);
%     avg_spea = zeros(nr_vox,1);
%     if pear == 1
%         for ind_vox = 1:nr_vox
%             avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
%             avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
%         end;
%     end;
%     if spea == 1
%         for ind_vox = 1:nr_vox      
%             avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
%             avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
%         end;
%     end;                
%     clear avg_corr_4D_pear avg_corr_4D_spea
%     if pear == 1
%         target_img = temp_img;
%         target_img.fileprefix = sprintf('CorrMaps%s_pear_par%d_between_contrasts.nii',str,i_par);
%         target_img.img = avg_pear;
%         save_nii(target_img,target_img.fileprefix); 
%     end;
%     if spea == 1
%         target_img = temp_img;
%         target_img.fileprefix = sprintf('CorrMaps%s_spea_par%d_between_contrasts.nii',str,i_par);
%         target_img.img = avg_spea;
%         save_nii(target_img,target_img.fileprefix);        
%     end;
%          
%   
%      end;
% end;       

% compare contrasts between sessions

for i_con = 1:2
    eval(sprintf('con=con%d;',i_con));
    eval(sprintf('con_count = con%d_count;',i_con));
    count_comp = 0;
    % create data matrix with data for all voxels in cols, subjects in rows
    % and runs in third dimension
    data = zeros(nr_subj,nr_vox,runs);
    for ind_runs = 1:runs
        for ind_subj = 1:nr_subj
            eval(sprintf('temp = img_%d_con%d(:,:,:,ind_subj);',ind_runs,i_con));
            data(ind_subj,:,ind_runs)=temp(:);
        end;
    end;
    clear temp;
    
    %correlations for each comparison
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                count_comp = count_comp +1;
                fprintf('...creates correlation maps for session %d and session %d...\n',i_run,i_run+i_sec);
                %load 4D images
                eval(sprintf('one=data(:,:,%d);',i_run));
                eval(sprintf('second=data(:,:,%d);',i_run+i_sec));
                
                r1 = zeros(nr_vox,1);
                zr1 = zeros(nr_vox,1);
                r2 = zeros(nr_vox,1);
                zr2 = zeros(nr_vox,1);
                for i_vox = 1:nr_vox
                    if pear == 1 
                        r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
                            if isnan (r(1,2))
                                r1(i_vox,1) = 0;
                                zr1(i_vox,1) = 0;
                            else
                                r1(i_vox,1) = r(1,2);
                                zr1(i_vox,1) = atanh(r(1,2));
                            end;
                    end
                    if spea == 1
                            r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                            r2(i_vox,1) = r;  
                            zr2(i_vox,1) = atanh(r);  
                    end;
                    
                end;
                clear one second
                
                % create correlation matrices for image 
                if pear == 1
                r_vec_pear_1_2 = reshape(r1,x,y,z);
                z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
                clear r1 zr1
                end;
                
                if spea == 1
                r_vec_spea_1_2 = reshape(r2,x,y,z);
                z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
                clear r2 zr2 
                end;
             

            % save correlation maps
            if pear == 1
            target_img = temp_img;

            file = sprintf('CorrMaps%s_pear_%d_%d_con%d.nii',str,i_run,i_run+i_sec,con_count);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_pear_1_2;
                else
                    r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_pear_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_pear_%d_%d_con%d.nii',str,i_run,i_run+i_sec,con_count);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_pear_1_2;
                else
                    z_r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_pear_1_2;
                end;            
            save_nii(target_img,target_img.fileprefix);      
            end;
            
            if spea == 1
            target_img = temp_img;
            file = sprintf('CorrMaps%s_spea_%d_%d_con%d.nii',str,i_run,i_run+i_sec,con_count);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_spea_1_2;
                else
                    r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_spea_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_spea_%d_%d_con%d.nii',str,i_run,i_run+i_sec,con_count);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_spea_1_2;
                else
                    z_r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_spea_1_2;
                end;               
            save_nii(target_img,target_img.fileprefix);                      
            end;
            
            %compute z mean and inverse Fisher's transformation
            if pear == 1
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear%s_%d_%d_con%d',str,i_run,i_run+i_sec,con_count);
            end;
            if spea == 1
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea%s_%d_%d_con%d',str,i_run,i_run+i_sec,con_count);    
            end;
            end; 
        end;
    end;
    clear r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2

    % calculate average correlation for all comparisons
    disp('...creates average correlation maps ...');   
    avg_corr_4D_pear = zeros(nr_vox,count_comp);
    avg_corr_4D_spea = zeros(nr_vox,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                ind_comp = ind_comp+1;
                if pear == 1
                map = load_nii(sprintf('z_CorrMaps%s_pear_%d_%d_con%d.nii',str,i_run,i_run+i_sec,con_count));
                avg_corr_4D_pear(:,ind_comp) = map.img(:);
                end;
                if spea == 1
                map = load_nii(sprintf('z_CorrMaps%s_spea_%d_%d_con%d.nii',str,i_run,i_run+i_sec,con_count));
                avg_corr_4D_spea(:,ind_comp) = map.img(:);
                end;                
            end; 
        end;
    end;
    clear map
    
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                

    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear_con%d.nii',str,con_count);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea_con%d.nii',str,con_count);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
clear avg_spea avg_pear avg_corr_4D_spea avg_corr_4D_pear z_r_vec_spea_1_2 

% if nr_para > 0    
%      count_comp = 0;
%      for i_par = 1:nr_para
%         data = zeros(nr_subj,nr_vox,runs);
%         for ind_runs = 1:runs
%             for ind_subj = 1:nr_subj
%                 eval(sprintf('temp = img_%s_par%d_%d(:,:,:,ind_subj);',con,i_par,ind_runs));
%                 data(ind_subj,:,ind_runs)=temp(:);
%             end;
%         end;
%         clear temp
%     
%     %correlations for each comparison
%     for i_run = 1:runs
%         for i_sec = 1:runs-i_run
%             if i_run+i_sec <= runs
%                 count_comp = count_comp +1;
%                 fprintf('...creates correlation maps for parametric contrast in session %d and session %d...\n',i_run,i_run+i_sec);
%                 %load 4D images
%                 eval(sprintf('one=data(:,:,%d);',i_run));
%                 eval(sprintf('second=data(:,:,%d);',i_run+i_sec));
%                 
%                 r1 = zeros(nr_vox,1);
%                 zr1 = zeros(nr_vox,1);
%                 r2 = zeros(nr_vox,1);
%                 zr2 = zeros(nr_vox,1);
%                 for i_vox = 1:nr_vox
%                     if pear == 1 
%                         r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
%                             if isnan (r(1,2))
%                                 r1(i_vox,1) = 0;
%                                 zr1(i_vox,1) = 0;
%                             else
%                                 r1(i_vox,1) = r(1,2);
%                                 zr1(i_vox,1) = atanh(r(1,2));
%                             end;
%                     end
%                     if spea == 1
%                             r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
%                             r2(i_vox,1) = r;  
%                             zr2(i_vox,1) = atanh(r);  
%                     end;
%                     
%                 end;
%                 clear one second
%                 
%                 % create correlation matrices for image 
%                 if pear == 1
%                     r_vec_pear_1_2 = reshape(r1,x,y,z);
%                     r_vec_spea_1_2 = reshape(r2,x,y,z);
%                     clear r1 zr1
%                 end;
%                 if spea == 1
%                     z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
%                     z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
%                     clear r2 zr2;
%                 end;
%              
% 
%             % save correlation maps
%             if pear == 1
%             target_img = temp_img;
% 
%             file = sprintf('CorrMaps%s_pear_par%d_%d_%d_con%d.nii',str,i_par,i_run,i_run+i_sec,con_count);
%             target_img.fileprefix = file;
%             if roi == 0
%                 target_img.img = r_vec_pear_1_2;
%             else
%                 r_vec_pear_1_2(~r_roi_ind)=0;
%                 target_img.img = r_vec_pear_1_2;
%             end;
%             save_nii(target_img,target_img.fileprefix);  
% 
%             target_img = temp_img;
%             file = sprintf('z_CorrMaps%s_pear_par%d_%d_%d_con%d.nii',str,i_par,i_run,i_run+i_sec,con_count);
%             target_img.fileprefix = file;
%             if roi == 0
%                 target_img.img = z_r_vec_pear_1_2;
%             else
%                 z_r_vec_pear_1_2(~r_roi_ind)=0;
%                 target_img.img = z_r_vec_pear_1_2;
%             end;            
%             save_nii(target_img,target_img.fileprefix);      
%             end;
%             
%             if spea == 1
%             target_img = temp_img;
%             file = sprintf('CorrMaps%s_spea_par%d_%d_%d_con%d.nii',str,i_par,i_run,i_run+i_sec,con_count);
%             target_img.fileprefix = file;
%             if roi == 0
%                 target_img.img = r_vec_spea_1_2;
%             else
%                 r_vec_spea_1_2(~r_roi_ind)=0;
%                 target_img.img = r_vec_spea_1_2;
%             end;
%             save_nii(target_img,target_img.fileprefix);
% 
%             target_img = temp_img;
%             file = sprintf('z_CorrMaps%s_spea_par%d_%d_%d_con%d.nii',str,i_par,i_run,i_run+i_sec,con_count);
%             target_img.fileprefix = file;
%             if roi == 0
%                 target_img.img = z_r_vec_spea_1_2;
%             else
%                 z_r_vec_spea_1_2(~r_roi_ind)=0;
%                 target_img.img = z_r_vec_spea_1_2;
%             end;               
%             save_nii(target_img,target_img.fileprefix);                      
%             end;
%             
%             %compute z mean and inverse Fisher's transformation
%             if pear == 1
%             mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
%             mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%             summary(1,end+1)=mean_r;
%             cols{1,end+1} = sprintf('mean_pear%s_par%d_%d_%d_con%d',str,i_par,i_run,i_run+i_sec,con_count);
%             end;
%             if spea == 1
%             mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
%             mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%             summary(1,end+1)=mean_r;
%             cols{1,end+1} = sprintf('mean_spea%s_par%d_%d_%d_con%d',str,i_par,i_run,i_run+i_sec,con_count);    
%             end;
%             clear z_r_vec_spea_1_2 r_vec_spea_1_2 r_vec_pear_1_2 z_r_vec_pear_1_2
%             end; 
%         end;
%     end;
%     
%     % calculate average correlation for all comparisons
%     disp('...creates average correlation maps for parametric contrast.../n');   
%     avg_corr_4D_pear = zeros(nr_vox,count_comp);
%     avg_corr_4D_spea = zeros(nr_vox,count_comp);
%     ind_comp = 0;
%     for i_run = 1:runs
%         for i_sec = 1:runs-i_run
%             if i_run+i_sec <= runs
%                 ind_comp = ind_comp+1;
%                 if pear == 1
%                 map = load_nii(sprintf('z_CorrMaps%s_pear_par%d_%d_%d_con%d.nii',str,i_par,i_run,i_run+i_sec,con_count));
%                 avg_corr_4D_pear(:,ind_comp) = map.img(:);
%                 end;
%                 if spea == 1
%                 map = load_nii(sprintf('z_CorrMaps%s_spea_par%d_%d_%d_con%d.nii',str,i_par,i_run,i_run+i_sec,con_count));
%                 avg_corr_4D_spea(:,ind_comp) = map.img(:);
%                 end;                
%             end; 
%         end;
%     end;
%     avg_pear = zeros(nr_vox,1);
%     avg_spea = zeros(nr_vox,1);
%     if pear == 1
%         for ind_vox = 1:nr_vox
%             avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
%             avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
%         end;
%     end;
%     if spea == 1
%         for ind_vox = 1:nr_vox      
%             avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
%             avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
%         end;
%     end;                
% 
%     if pear == 1
%         target_img = temp_img;
%         target_img.fileprefix = sprintf('CorrMaps%s_pear_par%d_con%d.nii',str,i_par,con_count);
%         target_img.img = avg_pear;
%         save_nii(target_img,target_img.fileprefix); 
%     end;
%     if spea == 1
%         target_img = temp_img;
%         target_img.fileprefix = sprintf('CorrMaps%s_spea_par%d_con%d.nii',str,i_par,con_count);
%         target_img.img = avg_spea;
%         save_nii(target_img,target_img.fileprefix);        
%     end;
%     clear avg_spea avg_pear avg_corr_4D_spea avg_corr_4D_pear z_r_vec_spea_1_2 
%      end;
% end; 
results_corr=dataset({summary(1,:),cols{:}});
assignin('base','results_corr',results_corr);
    
cd(dir_results);             
eval(sprintf('save results_corr%s_2cons.mat results_corr',str)); 
end;
           
end; 
elseif ex4D == 1
    if nr_cond == 1
    count_comp = 0;
    % create data matrix with data for all voxels in cols, subjects in rows
    % and runs in third dimension
    data = zeros(nr_subj,nr_vox,runs);
    for ind_runs = 1:runs
        for ind_subj = 1:nr_subj
            eval(sprintf('temp = FourD%d(:,:,:,ind_subj);',ind_runs));
            data(ind_subj,:,ind_runs)=temp(:);
        end;
    end;
    clear temp;
    
    %correlations for each comparison
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                count_comp = count_comp +1;
                fprintf('...creates correlation maps for session %d and session %d...\n',i_run,i_run+i_sec);
                %load 4D images
                eval(sprintf('one=data(:,:,%d);',i_run));
                eval(sprintf('second=data(:,:,%d);',i_run+i_sec));
                
                r1 = zeros(nr_vox,1);
                zr1 = zeros(nr_vox,1);
                r2 = zeros(nr_vox,1);
                zr2 = zeros(nr_vox,1);
                for i_vox = 1:nr_vox
                    if pear == 1 
                        r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
                            if isnan (r(1,2))
                                r1(i_vox,1) = 0;
                                zr1(i_vox,1) = 0;
                            else
                                r1(i_vox,1) = r(1,2);
                                zr1(i_vox,1) = atanh(r(1,2));
                            end;
                    end
                    if spea == 1
                            r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                            r2(i_vox,1) = r;  
                            zr2(i_vox,1) = atanh(r);  
                    end;
                    
                end;
                clear one second
                
                % create correlation matrices for image 
                if pear == 1
                r_vec_pear_1_2 = reshape(r1,x,y,z);
                z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
                clear r1 zr1
                end;
                
                if spea == 1
                r_vec_spea_1_2 = reshape(r2,x,y,z);
                z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
                clear r2 zr2 
                end;
             

            % save correlation maps
            if pear == 1
            target_img = temp_img;

            file = sprintf('CorrMaps%s_pear_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_pear_1_2;
                else
                    r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_pear_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_pear_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_pear_1_2;
                else
                    z_r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_pear_1_2;
                end;            
            save_nii(target_img,target_img.fileprefix);      
            end;
            
            if spea == 1
            target_img = temp_img;
            file = sprintf('CorrMaps%s_spea_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_spea_1_2;
                else
                    r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_spea_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_spea_%d_%d.nii',str,i_run,i_run+i_sec);
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_spea_1_2;
                else
                    z_r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_spea_1_2;
                end;               
            save_nii(target_img,target_img.fileprefix);                      
            end;
            
            %compute z mean and inverse Fisher's transformation
            if pear == 1
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear%s_%d_%d',str,i_run,i_run+i_sec);
            end;
            if spea == 1
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea%s_%d_%d',str,i_run,i_run+i_sec);    
            end;
            end; 
        end;
    end;
    clear r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2 data

    % calculate average correlation for all comparisons
    disp('...creates average correlation maps ...');   
    avg_corr_4D_pear = zeros(nr_vox,count_comp);
    avg_corr_4D_spea = zeros(nr_vox,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                ind_comp = ind_comp+1;
                if pear == 1
                map = load_nii(sprintf('z_CorrMaps%s_pear_%d_%d.nii',str,i_run,i_run+i_sec));
                avg_corr_4D_pear(:,ind_comp) = map.img(:);
                clear map
                end;
                if spea == 1
                map = load_nii(sprintf('z_CorrMaps%s_spea_%d_%d.nii',str,i_run,i_run+i_sec));
                avg_corr_4D_spea(:,ind_comp) = map.img(:);
                clear map
                end;                
            end; 
        end;
    end;
    
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                

    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear.nii',str);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea.nii',str);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
    results_corr=dataset({summary(1,:),cols{:}});
    assignin('base','results_corr',results_corr);

    cd(dir_results);             
    eval(sprintf('save results_corr%s.mat results_corr',str));
    clear avg_spea avg_pear avg_corr_4D_spea avg_corr_4D_pear z_r_vec_spea_1_2    
    elseif nr_cond > 1
     % compare conditions within sessions
     for i_run = 1:runs
            if runs == 1
                i_run = single_run;
            end;
        % create data matrix with data for all voxels in cols, subjects in rows
        % and 2 conditions in third dimension
        data = zeros(nr_subj,nr_vox,2);
        for ind_split = 1:2
            for ind_subj = 1:nr_subj
                eval(sprintf('temp = FourD%s%d(:,:,:,ind_subj);',conditions{ind_split,1},i_run));
                data(ind_subj,:,ind_split)=temp(:);
            end;
        end;
        clear temp
    
    %correlations for each split
    fprintf('...creates correlation maps for two conditions in session %d ...\n',i_run);
    %load 4D images
    one=data(:,:,1);
    second=data(:,:,2);

    r1 = zeros(nr_vox,1);
    zr1 = zeros(nr_vox,1);
    r2 = zeros(nr_vox,1);
    zr2 = zeros(nr_vox,1);
    for i_vox = 1:nr_vox
        if pear == 1 
            r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
                if isnan (r(1,2))
                    r1(i_vox,1) = 0;
                    zr1(i_vox,1) = 0;
                else
                    r1(i_vox,1) = r(1,2);
                    zr1(i_vox,1) = atanh(r(1,2));
                end;
        end
        if spea == 1
                r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                r2(i_vox,1) = r;  
                zr2(i_vox,1) = atanh(r);  
        end;
    end;
    clear one second
                
    % create correlation matrices for image 
    r_vec_pear_1_2 = reshape(r1,x,y,z);
    r_vec_spea_1_2 = reshape(r2,x,y,z);
    z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
    z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
    clear r1 zr1 r2 zr2;
             
    % save correlation maps
    if pear == 1
        target_img = temp_img;
        file = sprintf('CorrMaps%s_pear_%d_between_conditions.nii',str,i_run);
        target_img.fileprefix = file;
        if roi == 0
            target_img.img = r_vec_pear_1_2;
        else
            r_vec_pear_1_2(~r_roi_ind)=0;
            target_img.img = r_vec_pear_1_2;
        end;
        save_nii(target_img,target_img.fileprefix);  

        target_img = temp_img;
        file = sprintf('z_CorrMaps%s_pear_%d_between_conditions.nii',str,i_run);
        target_img.fileprefix = file;
        if roi == 0
            target_img.img = z_r_vec_pear_1_2;
        else
            z_r_vec_pear_1_2(~r_roi_ind)=0;
            target_img.img = z_r_vec_pear_1_2;
        end;            
        save_nii(target_img,target_img.fileprefix);      
    end;

    if spea == 1
    target_img = temp_img;
    file = sprintf('CorrMaps%s_spea_%d_between_conditions.nii',str,i_run);
    target_img.fileprefix = file;
    if roi == 0
        target_img.img = r_vec_spea_1_2;
    else
        r_vec_spea_1_2(~r_roi_ind)=0;
        target_img.img = r_vec_spea_1_2;
    end;
    save_nii(target_img,target_img.fileprefix);

    target_img = temp_img;
    file = sprintf('z_CorrMaps%s_spea_%d_between_conditions.nii',str,i_run);
    target_img.fileprefix = file;
    if roi == 0
        target_img.img = z_r_vec_spea_1_2;
    else
        z_r_vec_spea_1_2(~r_roi_ind)=0;
        target_img.img = z_r_vec_spea_1_2;
    end;               
    save_nii(target_img,target_img.fileprefix);                      
    end;

    %compute z mean and inverse Fisher's transformation
    if pear == 1
        mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_pear%s_%d_between_conditions',str,i_run);
    end;
    if spea == 1
        mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_spea%s_%d_between_conditions',str,i_run);    
    end;
clear z_r_vec_spea_1_2 r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2
     end;
     
    % calculate average correlation for all comparisons
    disp('...creates average correlation maps over all between condition correlations ...');   
    avg_corr_4D_pear = zeros(nr_vox,runs);
    avg_corr_4D_spea = zeros(nr_vox,runs);
    ind_comp = 0;
    for i_run = 1:runs
        ind_comp = ind_comp+1;
        if pear == 1
        map = load_nii(sprintf('z_CorrMaps%s_pear_%d_between_conditions.nii',str,i_run));
        avg_corr_4D_pear(:,ind_comp) = map.img(:);
        clear map
        end;
        if spea == 1
        map = load_nii(sprintf('z_CorrMaps%s_spea_%d_between_conditions.nii',str,i_run));
        avg_corr_4D_spea(:,ind_comp) = map.img(:);
        clear map
        end;                
    end;
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                
    clear avg_corr_4D_pear avg_corr_4D_spea
    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear_between_conditions.nii',str);
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea_between_conditions.nii',str);
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
% compare contrasts between sessions

for i_con = 1:2
    count_comp = 0;
    % create data matrix with data for all voxels in cols, subjects in rows
    % and runs in third dimension
    data = zeros(nr_subj,nr_vox,runs);
    for ind_runs = 1:runs
        for ind_subj = 1:nr_subj
            eval(sprintf('temp = FourD%s%d(:,:,:,ind_subj);',conditions{i_con,1},ind_runs));
            data(ind_subj,:,ind_runs)=temp(:);
        end;
    end;
    clear temp;
    
    %correlations for each comparison
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                count_comp = count_comp +1;
                fprintf('...creates correlation maps for session %d and session %d...\n',i_run,i_run+i_sec);
                %load 4D images
                eval(sprintf('one=data(:,:,%d);',i_run));
                eval(sprintf('second=data(:,:,%d);',i_run+i_sec));
                
                r1 = zeros(nr_vox,1);
                zr1 = zeros(nr_vox,1);
                r2 = zeros(nr_vox,1);
                zr2 = zeros(nr_vox,1);
                for i_vox = 1:nr_vox
                    if pear == 1 
                        r = corrcoef(second(:,i_vox),one(:,i_vox), 'rows', 'pairwise');
                            if isnan (r(1,2))
                                r1(i_vox,1) = 0;
                                zr1(i_vox,1) = 0;
                            else
                                r1(i_vox,1) = r(1,2);
                                zr1(i_vox,1) = atanh(r(1,2));
                            end;
                    end
                    if spea == 1
                            r = corr(second(:,i_vox),one(:,i_vox), 'Type','Spearman', 'rows', 'pairwise');
                            r2(i_vox,1) = r;  
                            zr2(i_vox,1) = atanh(r);  
                    end;
                    
                end;
                clear one second
                
                % create correlation matrices for image 
                if pear == 1
                r_vec_pear_1_2 = reshape(r1,x,y,z);
                z_r_vec_pear_1_2 = reshape(zr1,x,y,z);
                clear r1 zr1
                end;
                
                if spea == 1
                r_vec_spea_1_2 = reshape(r2,x,y,z);
                z_r_vec_spea_1_2 = reshape(zr2,x,y,z);
                clear r2 zr2 
                end;
             

            % save correlation maps
            if pear == 1
            target_img = temp_img;

            file = sprintf('CorrMaps%s_pear_%d_%d_%s.nii',str,i_run,i_run+i_sec,conditions{i_con,1});
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_pear_1_2;
                else
                    r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_pear_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);  

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_pear_%d_%d_%s.nii',str,i_run,i_run+i_sec,conditions{i_con,1});
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_pear_1_2;
                else
                    z_r_vec_pear_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_pear_1_2;
                end;            
            save_nii(target_img,target_img.fileprefix);      
            end;
            
            if spea == 1
            target_img = temp_img;
            file = sprintf('CorrMaps%s_spea_%d_%d_%s.nii',str,i_run,i_run+i_sec,conditions{i_con,1});
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = r_vec_spea_1_2;
                else
                    r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = r_vec_spea_1_2;
                end;
            save_nii(target_img,target_img.fileprefix);

            target_img = temp_img;
            file = sprintf('z_CorrMaps%s_spea_%d_%d_%s.nii',str,i_run,i_run+i_sec,conditions{i_con,1});
            target_img.fileprefix = file;
                if roi == 0
                    target_img.img = z_r_vec_spea_1_2;
                else
                    z_r_vec_spea_1_2(~r_roi_ind)=0;
                    target_img.img = z_r_vec_spea_1_2;
                end;               
            save_nii(target_img,target_img.fileprefix);                      
            end;
            
            %compute z mean and inverse Fisher's transformation
            if pear == 1
                mean_z = mean(z_r_vec_pear_1_2(~isinf(z_r_vec_pear_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_pear%s_%d_%d_%s',str,i_run,i_run+i_sec,conditions{i_con,1});
            end;
            if spea == 1
                mean_z = mean(z_r_vec_spea_1_2(~isinf(z_r_vec_spea_1_2)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_spea%s_%d_%d_%s',str,i_run,i_run+i_sec,conditions{i_con,1});    
            end;
            end; 
        end;
    end;
    clear r_vec_spea_1_2 z_r_vec_pear_1_2 r_vec_pear_1_2

    % calculate average correlation for all comparisons
    disp('...creates average correlation maps ...');   
    avg_corr_4D_pear = zeros(nr_vox,count_comp);
    avg_corr_4D_spea = zeros(nr_vox,count_comp);
    ind_comp = 0;
    for i_run = 1:runs
        for i_sec = 1:runs-i_run
            if i_run+i_sec <= runs
                ind_comp = ind_comp+1;
                if pear == 1
                map = load_nii(sprintf('z_CorrMaps%s_pear_%d_%d_%s.nii',str,i_run,i_run+i_sec,conditions{i_con,1}));
                avg_corr_4D_pear(:,ind_comp) = map.img(:);
                end;
                if spea == 1
                map = load_nii(sprintf('z_CorrMaps%s_spea_%d_%d_%s.nii',str,i_run,i_run+i_sec,conditions{i_con,1}));
                avg_corr_4D_spea(:,ind_comp) = map.img(:);
                end;                
            end; 
        end;
    end;
    clear map
    
    avg_pear = zeros(nr_vox,1);
    avg_spea = zeros(nr_vox,1);
    if pear == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_corr_4D_pear(ind_vox,:));
            avg_pear(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if spea == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_corr_4D_spea(ind_vox,:));
            avg_spea(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                

    if pear == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_pear_%s.nii',str,conditions{i_con,1});
        target_img.img = avg_pear;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if spea == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('CorrMaps%s_spea_%s.nii',str,conditions{i_con,1});
        target_img.img = avg_spea;
        save_nii(target_img,target_img.fileprefix);        
    end;
clear avg_spea avg_pear avg_corr_4D_spea avg_corr_4D_pear z_r_vec_spea_1_2 

results_corr=dataset({summary(1,:),cols{:}});
assignin('base','results_corr',results_corr);
    
cd(dir_results);             
eval(sprintf('save results_corr%s_conditions.mat results_corr',str));         
    end;
end;
disp('...finished correlation maps')
end;
end;
%% ICCs
if cons == 1 || abs == 1
    %initiating variables
    ICC_con_ROI=zeros(nr_vox,1);
    ICC_abs_ROI=zeros(nr_vox,1);
    z_ICC_con_ROI=zeros(nr_vox,1);
    z_ICC_abs_ROI=zeros(nr_vox,1);
    summary = [];
    cols = {};
    disp('...calculating ICCs...')
    if exStats == 1
    if split == 0 && two_cons == 0
            %create data matrix
            data=zeros(nr_vox,nr_subj,runs);  
        for ind_runs = 1:runs
             for ind_subj = 1:nr_subj
                 eval(sprintf('temp = img_%d(:,:,:,ind_subj);',ind_runs));
                 data(:,ind_subj,ind_runs)=temp(:);
             end;
        end;
        clear temp
        disp('...over all sessions...')
        %calculate ICCs
        for ind_voxel = 1:nr_vox
            nsamples=nr_subj*runs;
 
            grandmean=0;
            for sub=1:nr_subj,     
                for sess=1:runs,
                    grandmean= grandmean + data(ind_voxel,sub,sess);
                end
            end;
            grandmean=grandmean./nsamples;
 
            sessionmean=zeros(runs,1);
            for sess=1:runs
                for sub=1:nr_subj,  
                    sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
                end
                sessionmean(sess)=sessionmean(sess)./nr_subj;
            end
 
            subjmean=zeros(nr_subj,1);
            for sub=1:nr_subj
                for sess=1:runs
                    subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                    WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                    EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                voxICC_con=(BMS-EMS)./(BMS+(runs-1).*EMS); 
                ICC_con_ROI(ind_voxel,1) = voxICC_con;
                z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
            end;
                     
           %absolute agreement 
           if abs==1
               voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
               runs.* (JMS-EMS)./nr_subj);
               ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
               z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
          end;
       end;

 disp('...saving ICC image over all sessions...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s.nii',str);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s.nii',str);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s.nii',str);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s.nii',str);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 disp('...computing average ICC over all sessions...')
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s',str);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s',str);    
    end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
if nr_para > 0
    for i_par = 1:nr_para
    %initiating variables
    ICC_con_ROI=zeros(numel(img_1(:,:,:,1)),1);
    ICC_abs_ROI=zeros(numel(img_1(:,:,:,1)),1);
    z_ICC_con_ROI=zeros(numel(img_1(:,:,:,1)),1);
    z_ICC_abs_ROI=zeros(numel(img_1(:,:,:,1)),1);
    disp('...calculating ICCs for parametric contrast over all sessions...')
    
            %create data matrix
            data=zeros(nr_vox,nr_subj,runs);  
        for ind_runs = 1:runs
             for ind_subj = 1:nr_subj
                 eval(sprintf('temp = img_par%d_%d(:,:,:,ind_subj);',i_par,ind_runs));
                 data(:,ind_subj,ind_runs)=temp(:);
             end;
        end;
        clear temp
        
        %calculate ICCs
        for ind_voxel = 1:nr_vox
            nsamples=nr_subj*runs;
 
            grandmean=0;
            for sub=1:nr_subj,     
                for sess=1:runs,
                    grandmean= grandmean + data(ind_voxel,sub,sess);
                end
            end;
            grandmean=grandmean./nsamples;
 
            sessionmean=zeros(runs,1);
            for sess=1:runs
                for sub=1:nr_subj,  
                    sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
                end
                sessionmean(sess)=sessionmean(sess)./nr_subj;
            end
 
            subjmean=zeros(nr_subj,1);
            for sub=1:nr_subj
                for sess=1:runs
                    subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                    WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                    EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                voxICC_con=(BMS-EMS)./(BMS+(runs-1).*EMS); 
                ICC_con_ROI(ind_voxel,1) = voxICC_con;
                z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
            end;
                     
           %absolute agreement 
           if abs==1
               voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
               runs.* (JMS-EMS)./nr_subj);
               ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
               z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
          end;
       end;

 disp('...saving ICC images for parametric contrast over all sessions...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s_par%d.nii',str,i_par);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s_par%d.nii',str,i_par);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s_par%d.nii',str,i_par);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s_par%d.nii',str,i_par);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 disp('...computing average ICC for parametric contrast over all sessions...')
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s_par%d',str,i_par);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s_par%d',str,i_par);    
    end;
clear data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
    end;
end;

%calculation ICC for each comparison
 if runs > 1
    %initiating variables
    ICC_con_ROI=zeros(nr_vox,1);
    ICC_abs_ROI=zeros(nr_vox,1);
    z_ICC_con_ROI=zeros(nr_vox,1);
    z_ICC_abs_ROI=zeros(nr_vox,1);
   
    %create data matrix
    data1=zeros(nr_vox,nr_subj,runs);  
    for ind_runs = 1:runs
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = img_%d(:,:,:,ind_subj);',ind_runs));
             data1(:,ind_subj,ind_runs)=temp(:);
         end;
    end;
    clear temp 
    
    for ind_run = 1:runs
        for ind_sec = 1:runs-1
             if ind_run + ind_sec <= runs
                 fprintf('...calculate ICC for session %d and session %d...\n',ind_run,ind_run+ind_sec);
                for ind_vox = 1:nr_vox
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        grandmean= grandmean + data1(ind_vox,sub,ind_run);
                        grandmean= grandmean + data1(ind_vox,sub,ind_run+ind_sec);
                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        if sess == 1
                            temp_ind = ind_run;
                        else
                            temp_ind = ind_run+ind_sec;
                        end;
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + data1(ind_vox,sub,temp_ind);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                            
                            subjmean(sub)=subjmean(sub) + data1(ind_vox,sub,temp_ind);
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
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                             
                            WMS = WMS + (data1(ind_vox,sub,temp_ind)-subjmean(sub)).^2;
                            EMS = EMS + (data1(ind_vox,sub,temp_ind)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                     
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
                    ICC_con_ROI(ind_vox,1) = voxICC_con;
                    z_ICC_con_ROI(ind_vox,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI(ind_vox,1) = voxICC_abs;
                    z_ICC_abs_ROI(ind_vox,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;                    
                 end;
                 
         disp('...saving ICC images...')
         ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
         ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
         z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
         z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
         if roi == 1
             ICC_con_ROI(~r_roi_ind)=0;
             z_ICC_con_ROI(~r_roi_ind)=0;
             ICC_abs_ROI(~r_roi_ind)=0;
             z_ICC_abs_ROI(~r_roi_ind)=0;
         end;     

         % save ICC maps
         target_img = temp_img;
         file1 = sprintf('ICC%s_con_%d_%d.nii',str,ind_run,ind_run+ind_sec);
         target_img.fileprefix = file1;
         target_img.img = ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
         target_img = temp_img;
         file2 = sprintf('ICC%s_abs_%d_%d.nii',str,ind_run,ind_run+ind_sec);
         target_img.fileprefix = file2;
         target_img.img = ICC_abs_ROI;
         save_nii(target_img,target_img.fileprefix); 

         target_img = temp_img;
         file3 = sprintf('z_ICC%s_con_%d_%d.nii',str,ind_run,ind_run+ind_sec);
         target_img.fileprefix = file3;
         target_img.img = z_ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
        target_img = temp_img;
        file4 = sprintf('z_ICC%s_abs_%d_%d.nii',str,ind_run,ind_run+ind_sec);
        target_img.fileprefix = file4;
        target_img.img = z_ICC_abs_ROI;
        save_nii(target_img,target_img.fileprefix);   
                             
        % computing average ICC
                 fprintf('...calculate average ICC for session %d and session %d...\n',ind_run,ind_run+ind_sec);
            %compute z mean and inverse Fisher's transformation
            if abs == 1
                mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCabs%s_%d_%d',str,ind_run,ind_run+ind_sec);
            end;
            if cons == 1
                mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCcon%s_%d_%d',str,ind_run,ind_run+ind_sec);    
            end;
            clear ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj

             end;
        end;
    end;
clear data1;
        
if nr_para > 0
     fprintf('...ICC parametric modulator for each comparison...\n')
     for ind_para = 1:nr_para

   
    %create data matrix
    data2=zeros(nr_vox,nr_subj,runs);  
    for ind_runs = 1:runs
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = img_par%d_%d(:,:,:,ind_subj);',ind_para,ind_runs));
             data2(:,ind_subj,ind_runs)=temp(:);
         end;
    end;
    clear temp 
    
    for ind_run = 1:runs
        for ind_sec = 1:runs-1
             if ind_run + ind_sec <= runs
              fprintf('...calculate ICC of parametric contrast for session %d and session %d...\n',ind_run,ind_run+ind_sec);
            %initiating variables
            ICC_con_ROI=zeros(nr_vox,1);
            ICC_abs_ROI=zeros(nr_vox,1);
            z_ICC_con_ROI=zeros(nr_vox,1);
            z_ICC_abs_ROI=zeros(nr_vox,1);
                for ind_vox = 1:nr_vox
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        grandmean= grandmean + data2(ind_vox,sub,ind_run);
                        grandmean= grandmean + data2(ind_vox,sub,ind_run+ind_sec);
                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        if sess == 1
                            temp_ind = ind_run;
                        else
                            temp_ind = ind_run+ind_sec;
                        end;
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + data2(ind_vox,sub,temp_ind);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                            
                            subjmean(sub)=subjmean(sub) + data2(ind_vox,sub,temp_ind);
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
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                             
                            WMS = WMS + (data2(ind_vox,sub,temp_ind)-subjmean(sub)).^2;
                            EMS = EMS + (data2(ind_vox,sub,temp_ind)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                     
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
                    ICC_con_ROI(ind_vox,1) = voxICC_con;
                    z_ICC_con_ROI(ind_vox,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI(ind_vox,1) = voxICC_abs;
                    z_ICC_abs_ROI(ind_vox,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;                    
                 end;
                 
         disp('...saving ICC images...')
         ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
         ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
         z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
         z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
         if roi == 1
             ICC_con_ROI(~r_roi_ind)=0;
             z_ICC_con_ROI(~r_roi_ind)=0;
             ICC_abs_ROI(~r_roi_ind)=0;
             z_ICC_abs_ROI(~r_roi_ind)=0;
         end;     

         % save ICC maps
         target_img = temp_img;
         file1 = sprintf('ICC%s_con_%d_%d_par%d.nii',str,ind_run,ind_run+ind_sec,ind_para);
         target_img.fileprefix = file1;
         target_img.img = ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
         target_img = temp_img;
         file2 = sprintf('ICC%s_abs_%d_%d_par%d.nii',str,ind_run,ind_run+ind_sec,ind_para);
         target_img.fileprefix = file2;
         target_img.img = ICC_abs_ROI;
         save_nii(target_img,target_img.fileprefix); 

         target_img = temp_img;
         file3 = sprintf('z_ICC%s_con_%d_%d_par%d.nii',str,ind_run,ind_run+ind_sec,ind_para);
         target_img.fileprefix = file3;
         target_img.img = z_ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
        target_img = temp_img;
        file4 = sprintf('z_ICC%s_abs_%d_%d_par%d.nii',str,ind_run,ind_run+ind_sec,ind_para);
        target_img.fileprefix = file4;
        target_img.img = z_ICC_abs_ROI;
        save_nii(target_img,target_img.fileprefix);   
                             
        % computing average ICC
                 fprintf('...calculate average ICC for parametric contrast in session %d and session %d...\n',ind_run,ind_run+ind_sec);
            %compute z mean and inverse Fisher's transformation
            if abs == 1
                mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCabs%s_%d_%d_par%d',str,ind_run,ind_run+ind_sec,ind_para);
            end;
            if cons == 1
                mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCcon%s_%d_%d_par%d',str,ind_run,ind_run+ind_sec,ind_para);    
            end;
            clear ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
             end;
        end;
    end;
    clearvars data2 

end;
end;
end;

    assignin('base','summary',summary);
    assignin('base','cols',cols);

    results_ICC=dataset({summary(1,:),cols{:}});
    assignin('base','results_ICC',results_ICC);

    cd(dir_results);             
    eval(sprintf('save results_ICC%s.mat results_ICC',str));

%% based on split-half
elseif split == 1
    %create data matrix
    for ind_runs = 1:runs
    data=zeros(nr_vox,nr_subj,2);  
    for ind_split = 1:2
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = img_%d_split%d(:,:,:,ind_subj);',ind_runs,ind_split));
             data(:,ind_subj,ind_split)=temp(:);
         end;
    end;
    clear temp
        
    %calculate ICCs
    fprintf('...calculate ICCs for splitted session %d...\n',ind_runs)
    for ind_voxel = 1:nr_vox
        nsamples=nr_subj*2;

        grandmean=0;
        for sub=1:nr_subj,     
            for sess=1:2,
                grandmean= grandmean + data(ind_voxel,sub,sess);
            end
        end;
        grandmean=grandmean./nsamples;

        sessionmean=zeros(2,1);
        for sess=1:2
            for sub=1:nr_subj,  
                sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
            end
            sessionmean(sess)=sessionmean(sess)./nr_subj;
        end

        subjmean=zeros(nr_subj,1);
        for sub=1:nr_subj
            for sess=1:2
                subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
            voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
            ICC_con_ROI(ind_voxel,1) = voxICC_con;
            z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
        end;

       %absolute agreement 
       if abs==1
           voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
           2.* (JMS-EMS)./nr_subj);
           ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
           z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
      end;
   end;

 disp('...saving ICC images...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s_%d_split.nii',str,ind_runs);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s_%d_split.nii',str,ind_runs);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s_%d_split.nii',str,ind_runs);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s_%d_split.nii',str,ind_runs);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 disp('...computing average ICC for each split...')
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s_%d_split',str,ind_runs);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s_%d_split',str,ind_runs);    
    end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj target_img
if nr_para > 0
    for i_par = 1:nr_para
    %initiating variables
    ICC_con_ROI=zeros(nr_vox,1);
    ICC_abs_ROI=zeros(nr_vox,1);
    z_ICC_con_ROI=zeros(nr_vox,1);
    z_ICC_abs_ROI=zeros(nr_vox,1);
    fprintf('...calculating ICCs for splitted parametric regressor in session %d...\n',ind_runs)
    
    %create data matrix
    data=zeros(nr_vox,nr_subj,2);  
    for ind_split = 1:2
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = img_%d_par%d_split%d(:,:,:,ind_subj);',ind_runs,i_par,ind_split));
             data(:,ind_subj,ind_split)=temp(:);
         end;
    end;
        clear temp
        
        %calculate ICCs
        for ind_voxel = 1:nr_vox
            nsamples=nr_subj*2;

            grandmean=0;
            for sub=1:nr_subj,     
                for sess=1:2,
                    grandmean= grandmean + data(ind_voxel,sub,sess);
                end
            end;
            grandmean=grandmean./nsamples;

            sessionmean=zeros(2,1);
            for sess=1:2
                for sub=1:nr_subj,  
                    sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
                end
                sessionmean(sess)=sessionmean(sess)./nr_subj;
            end

            subjmean=zeros(nr_subj,1);
            for sub=1:nr_subj
                for sess=1:2
                    subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                    WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                    EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
                ICC_con_ROI(ind_voxel,1) = voxICC_con;
                z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
            end;

           %absolute agreement 
           if abs==1
               voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
               2.* (JMS-EMS)./nr_subj);
               ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
               z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
          end;
       end;

 disp('...saving ICC images...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s_par%d_%d_split.nii',str,i_par,ind_runs);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s_par%d_%d_split.nii',str,i_par,ind_runs);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s_par%d_%d_split.nii',str,i_par,ind_runs);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s_par%d_%d_split.nii',str,i_par,ind_runs);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 fprintf('...computing average ICC for parametric regressor in session %d...',ind_runs)
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s_par%d_%d_split',str,i_par,ind_runs);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s_par%d_%d_split',str,i_par,ind_runs);    
    end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
    end;
    end;    
    end;
    
 % calculate voxelwise average ICC over all splits
 disp('...creates average ICC maps over all splits ...');   
    avg_ICC_4D_abs = zeros(nr_vox,runs);
    avg_ICC_4D_con = zeros(nr_vox,runs);
    for i_run = 1:runs
        if abs == 1
        map = load_nii(sprintf('z_ICC_abs%s_%d_split.nii',str,i_run));
        avg_ICC_4D_abs(:,i_run) = map.img(:);
        clear map
        end;
        if cons == 1
        map = load_nii(sprintf('z_ICC_con%s_%d_split.nii',str,i_run));
        avg_ICC_4D_con(:,i_run) = map.img(:);
        clear map
        end;                
    end;
    avg_abs = zeros(nr_vox,1);
    avg_con = zeros(nr_vox,1);
    if abs == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_ICC_4D_abs(ind_vox,:));
            avg_abs(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if cons == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_ICC_4D_con(ind_vox,:));
            avg_con(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                
    clear avg_ICC_4D_abs avg_ICC_4D_con avg_temp
    if abs == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_abs%s_split.nii',str);
        target_img.img = avg_abs;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if cons == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_con%s_split.nii',str);
        target_img.img = avg_con;
        save_nii(target_img,target_img.fileprefix);        
    end;
    clear avg_abs avg_con 
    
    assignin('base','summary',summary);
    assignin('base','cols',cols);

    results_ICC=dataset({summary(1,:),cols{:}});
    assignin('base','results_ICC',results_ICC);

    cd(dir_results);             
    eval(sprintf('save results_ICC%s_split.mat results_ICC',str));

%% two contrasts out of one statistic
elseif two_cons == 1
    %create data matrix
    for ind_runs = 1:runs
    data=zeros(nr_vox,nr_subj,2);  
    for ind_con = 1:2
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = img_%d_con%d(:,:,:,ind_subj);',ind_runs,ind_con));
             data(:,ind_subj,ind_con)=temp(:);
         end;
    end;
    clear temp
        
    %calculate ICCs
    fprintf('...calculate ICC between conditions for session %d...',ind_runs);

    for ind_voxel = 1:nr_vox
        nsamples=nr_subj*2;

        grandmean=0;
        for sub=1:nr_subj,     
            for sess=1:2,
                grandmean= grandmean + data(ind_voxel,sub,sess);
            end
        end;
        grandmean=grandmean./nsamples;

        sessionmean=zeros(2,1);
        for sess=1:2
            for sub=1:nr_subj,  
                sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
            end
            sessionmean(sess)=sessionmean(sess)./nr_subj;
        end

        subjmean=zeros(nr_subj,1);
        for sub=1:nr_subj
            for sess=1:2
                subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
            voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
            ICC_con_ROI(ind_voxel,1) = voxICC_con;
            z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
        end;

       %absolute agreement 
       if abs==1
           voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
           2.* (JMS-EMS)./nr_subj);
           ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
           z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
      end;
   end;

 disp('...saving ICC images...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s_%d_between_contrasts.nii',str,ind_runs);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s_%d_between_contrasts.nii',str,ind_runs);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s_%d_between_contrasts.nii',str,ind_runs);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s_%d_between_contrasts.nii',str,ind_runs);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 fprintf('...computing average ICC between conditions in session %d...\n',ind_runs)
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s_%d_between_contrasts',str,ind_runs);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s_%d_between_contrasts',str,ind_runs);    
    end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
% if nr_para > 0
%     for i_par = 1:nr_para
%     %initiating variables
%     ICC_con_ROI=zeros(nr_vox,1);
%     ICC_abs_ROI=zeros(nr_vox,1);
%     z_ICC_con_ROI=zeros(nr_vox,1);
%     z_ICC_abs_ROI=zeros(nr_vox,1);
%     fprintf('...calculating ICCs between parametric regressors of contrasts in session %d...',ind_runs)
%     
%     %create data matrix
%     data=zeros(nr_vox,nr_subj,2);  
%     for ind_con = 1:2
%          for ind_subj = 1:nr_subj
%              eval(sprintf('temp = img_%d_par%d_con%d(:,:,:,ind_subj);',ind_runs,i_par,ind_con));
%              data(:,ind_subj,ind_runs)=temp(:);
%          end;
%     end;
%         clear temp
%         
%         %calculate ICCs
%         for ind_voxel = 1:nr_vox
%             nsamples=nr_subj*2;
% 
%             grandmean=0;
%             for sub=1:nr_subj,     
%                 for sess=1:2,
%                     grandmean= grandmean + data(ind_voxel,sub,sess);
%                 end
%             end;
%             grandmean=grandmean./nsamples;
% 
%             sessionmean=zeros(2,1);
%             for sess=1:2
%                 for sub=1:nr_subj,  
%                     sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
%                 end
%                 sessionmean(sess)=sessionmean(sess)./nr_subj;
%             end
% 
%             subjmean=zeros(nr_subj,1);
%             for sub=1:nr_subj
%                 for sess=1:2
%                     subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
%                 end
%                 subjmean(sub)=subjmean(sub)./2;
%             end
% 
%             % mean squares
%             BMS=0; % between subject
%             WMS=0; % within subject 
%             EMS=0; % error
%             JMS=0; % session
% 
%             for sub=1:nr_subj,    
%                 BMS = BMS + (subjmean(sub)-grandmean).^2;
%                 for sess=1:2
%                     WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
%                     EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
%                 end
%             end;
% 
%             for sess=1:2
%                 JMS=  JMS + (sessionmean(sess)-grandmean).^2;
%             end;
% 
%             %define the true value of the mean square.
%             BMS= 2.*BMS./(nr_subj-1);
%             WMS= WMS./(2-1)./nr_subj;
%             JMS= nr_subj.*JMS./(2-1);
%             EMS= EMS./(2-1)./(nr_subj-1); 
% 
%             %consistency agreement  
%             if cons==1
%                 voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
%                 ICC_con_ROI(ind_voxel,1) = voxICC_con;
%                 z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
%             end;
% 
%            %absolute agreement 
%            if abs==1
%                voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
%                2.* (JMS-EMS)./nr_subj);
%                ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
%                z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
%           end;
%        end;
% 
%  disp('...saving ICC images...')
%  ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
%  ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
%  z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
%  z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
%  if roi == 1
%      ICC_con_ROI(~r_roi_ind)=0;
%      z_ICC_con_ROI(~r_roi_ind)=0;
%      ICC_abs_ROI(~r_roi_ind)=0;
%      z_ICC_abs_ROI(~r_roi_ind)=0;
%  end;
%  
%  % save ICC maps
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('ICC_con%s_par%d_%d_between_contrasts.nii',str,i_par,ind_runs);
%  target_img.img = ICC_con_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
%  clear target_img;
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('ICC_abs%s_par%d_%d_between_contrasts.nii',str,i_par,ind_runs);
%  target_img.img = ICC_abs_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('z_ICC_con%s_par%d_%d_between_contrasts.nii',str,i_par,ind_runs);
%  target_img.img = z_ICC_con_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
%  clear target_img;
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('z_ICC_abs%s_par%d_%d_between_contrasts.nii',str,i_par,ind_runs);
%  target_img.img = z_ICC_abs_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
% % computing average ICC
%  fprintf('...computing average ICC for parametric regressor between contrasts in session %d...',ind_runs)
%     %compute z mean and inverse Fisher's transformation
%     if abs == 1
%         mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
%         mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%         summary(1,end+1)=mean_r;
%         cols{1,end+1} = sprintf('mean_ICCabs%s_par%d_%d_between_contrasts',str,i_par,ind_runs);
%     end;
%     if cons == 1
%         mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
%         mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%         summary(1,end+1)=mean_r;
%         cols{1,end+1} = sprintf('mean_ICCcon%s_par%d_%d_between_contrasts',str,i_par,ind_runs);    
%     end;
% clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
%     end;
%     end;    
    end;
    
 % calculate voxelwise average ICC over all splits
 disp('...creates average ICC maps over all contrast comparisons of parametric regressors...');   
    avg_ICC_4D_abs = zeros(nr_vox,runs);
    avg_ICC_4D_con = zeros(nr_vox,runs);
    for i_run = 1:runs
        if abs == 1
        map = load_nii(sprintf('z_ICC_abs%s_%d_between_contrasts.nii',str,i_run));
        avg_ICC_4D_abs(:,i_run) = map.img(:);
        clear map
        end;
        if cons == 1
        map = load_nii(sprintf('z_ICC_con%s_%d_between_contrasts.nii',str,i_run));
        avg_ICC_4D_con(:,i_run) = map.img(:);
        clear map
        end;                
    end;
    avg_abs = zeros(nr_vox,1);
    avg_con = zeros(nr_vox,1);
    if abs == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_ICC_4D_abs(ind_vox,:));
            avg_abs(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if cons == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_ICC_4D_con(ind_vox,:));
            avg_con(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                
    clear avg_ICC_4D_abs avg_ICC_4D_con avg_temp
    if abs == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_abs%s_between_contrasts.nii',str);
        target_img.img = avg_abs;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if cons == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_con%s_between_contrasts.nii',str);
        target_img.img = avg_con;
        save_nii(target_img,target_img.fileprefix);        
    end;
    clear avg_abs avg_con 
    
    % ICCs within contrast over sessions
for ind_con = 1:2
    fprintf('...calculating ICC within contrast %d over all sessions...\n',ind_con)
    eval(sprintf('con=con%d;',ind_con));
    eval(sprintf('con_count = con%d_count;',ind_con));
    count_comp = 0;  
%initiating variables
    ICC_con_ROI=zeros(numel(img_1_con1(:,:,:,1)),1);
    ICC_abs_ROI=zeros(numel(img_1_con1(:,:,:,1)),1);
    z_ICC_con_ROI=zeros(numel(img_1_con1(:,:,:,1)),1);
    z_ICC_abs_ROI=zeros(numel(img_1_con1(:,:,:,1)),1);

    disp('...calculating ICCs...')
    
           %create data matrix
            data=zeros(nr_vox,nr_subj,runs);  
        for ind_runs = 1:runs
             for ind_subj = 1:nr_subj
                 eval(sprintf('temp = img_%d_con%d(:,:,:,ind_subj);',ind_runs,ind_con));
                 data(:,ind_subj,ind_runs)=temp(:);
             end;
        end;
        clear temp
        
        %calculate ICCs
        for ind_voxel = 1:nr_vox
            nsamples=nr_subj*runs;
 
            grandmean=0;
            for sub=1:nr_subj,     
                for sess=1:runs,
                    grandmean= grandmean + data(ind_voxel,sub,sess);
                end
            end;
            grandmean=grandmean./nsamples;
 
            sessionmean=zeros(runs,1);
            for sess=1:runs
                for sub=1:nr_subj,  
                    sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
                end
                sessionmean(sess)=sessionmean(sess)./nr_subj;
            end
 
            subjmean=zeros(nr_subj,1);
            for sub=1:nr_subj
                for sess=1:runs
                    subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                    WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                    EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                voxICC_con=(BMS-EMS)./(BMS+(runs-1).*EMS); 
                ICC_con_ROI(ind_voxel,1) = voxICC_con;
                z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
            end;
                     
           %absolute agreement 
           if abs==1
               voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
               runs.* (JMS-EMS)./nr_subj);
               ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
               z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
          end;
       end;

 disp('...saving ICC images...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s_con%d.nii',str,con_count);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s_con%d.nii',str,con_count);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s_con%d.nii',str,con_count);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s_con%d.nii',str,con_count);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 disp('...computing average ICC...')
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s_con%d',str,con_count);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s_con%d',str,con_count);    
    end;
end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
% if nr_para > 0
%     fprintf('...calculating ICC for parametric regressor within contrast %d over all sessions...',ind_con)
% 
%     for i_par = 1:nr_para
%     %initiating variables
%     ICC_con_ROI=zeros(numel(img_1(:,:,:,1)),1);
%     ICC_abs_ROI=zeros(numel(img_1(:,:,:,1)),1);
%     z_ICC_con_ROI=zeros(numel(img_1(:,:,:,1)),1);
%     z_ICC_abs_ROI=zeros(numel(img_1(:,:,:,1)),1);
%     disp('...calculating ICCs...')
%     
%             %create data matrix
%             data=zeros(nr_vox,nr_subj,runs);  
%         for ind_runs = 1:runs
%              for ind_subj = 1:nr_subj
%                  eval(sprintf('temp = img_con%d_par%d_%d(:,:,:,ind_subj);',ind_con,i_par,ind_runs));
%                  data(:,ind_subj,ind_runs)=temp(:);
%              end;
%         end;
%         clear temp
%         
%         %calculate ICCs
%         for ind_voxel = 1:nr_vox
%             nsamples=nr_subj*runs;
%  
%             grandmean=0;
%             for sub=1:nr_subj,     
%                 for sess=1:runs,
%                     grandmean= grandmean + data(ind_voxel,sub,sess);
%                 end
%             end;
%             grandmean=grandmean./nsamples;
%  
%             sessionmean=zeros(runs,1);
%             for sess=1:runs
%                 for sub=1:nr_subj,  
%                     sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
%                 end
%                 sessionmean(sess)=sessionmean(sess)./nr_subj;
%             end
%  
%             subjmean=zeros(nr_subj,1);
%             for sub=1:nr_subj
%                 for sess=1:runs
%                     subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
%                 end
%                 subjmean(sub)=subjmean(sub)./runs;
%             end
%             
%             % mean squares
%             BMS=0; % between subject
%             WMS=0; % within subject 
%             EMS=0; % error
%             JMS=0; % session
%                      
%             for sub=1:nr_subj,    
%                 BMS = BMS + (subjmean(sub)-grandmean).^2;
%                 for sess=1:runs
%                     WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
%                     EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
%                 end
%             end;
%  
%             for sess=1:runs
%                 JMS=  JMS + (sessionmean(sess)-grandmean).^2;
%             end;
%  
%             %define the true value of the mean square.
%             BMS= runs.*BMS./(nr_subj-1);
%             WMS= WMS./(runs-1)./nr_subj;
%             JMS= nr_subj.*JMS./(runs-1);
%             EMS= EMS./(runs-1)./(nr_subj-1); 
%  
%             %consistency agreement  
%             if cons==1
%                 voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
%                 ICC_con_ROI(ind_voxel,1) = voxICC_con;
%                 z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
%             end;
%                      
%            %absolute agreement 
%            if abs==1
%                voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
%                runs.* (JMS-EMS)./nr_subj);
%                ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
%                z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
%           end;
%        end;
% 
%  disp('...saving ICC images...')
%  ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
%  ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
%  z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
%  z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
%  if roi == 1
%      ICC_con_ROI(~r_roi_ind)=0;
%      z_ICC_con_ROI(~r_roi_ind)=0;
%      ICC_abs_ROI(~r_roi_ind)=0;
%      z_ICC_abs_ROI(~r_roi_ind)=0;
%  end;
%  
%  % save ICC maps
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('ICC_con%s_par%d_con%d.nii',str,i_par,con_count);
%  target_img.img = ICC_con_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
%  clear target_img;
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('ICC_abs%s_par%d_con%d.nii',str,i_par,con_count);
%  target_img.img = ICC_abs_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('z_ICC_con%s_par%d_con%d.nii',str,i_par,con_count);
%  target_img.img = z_ICC_con_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
%  clear target_img;
%  target_img = temp_img;
%  target_img.fileprefix = sprintf('z_ICC_abs%s_par%d_con%d.nii',str,i_par,con_count);
%  target_img.img = z_ICC_abs_ROI;
%  save_nii(target_img,target_img.fileprefix); 
%  
% % computing average ICC
%  disp('...computing average ICC within parametric contrast...')
%     %compute z mean and inverse Fisher's transformation
%     if abs == 1
%         mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
%         mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%         summary(1,end+1)=mean_r;
%         cols{1,end+1} = sprintf('mean_ICCabs%s_par%d_con%d',str,i_par,con_count);
%     end;
%     if cons == 1
%         mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
%         mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%         summary(1,end+1)=mean_r;
%         cols{1,end+1} = sprintf('mean_ICCcon%s_par%d_con%d',str,i_par,con_count);
%     end;
% clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
%     end;
%     end;

%calculation ICC for each comparison
 if runs > 1
    for ind_cons = 1:2
        eval(sprintf('con_count = con%d_count;',ind_cons));
    %initiating variables
    ICC_con_ROI=zeros(nr_vox,1);
    ICC_abs_ROI=zeros(nr_vox,1);
    z_ICC_con_ROI=zeros(nr_vox,1);
    z_ICC_abs_ROI=zeros(nr_vox,1);
   
    %create data matrix
    data=zeros(nr_vox,nr_subj,runs);  
    for ind_runs = 1:runs
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = img_%d_con%d(:,:,:,ind_subj);',ind_runs,ind_cons));
             data(:,ind_subj,ind_runs)=temp(:);
         end;
    end;
    clear temp 
    
    for ind_run = 1:runs
        for ind_sec = 1:runs-1
             if ind_run + ind_sec <= runs
    fprintf('...calculating ICC within contrast %d between session %d and %d...',con_count,ind_run, ind_run+ind_sec)
                 
                for ind_vox = 1:nr_vox
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        grandmean= grandmean + data(ind_vox,sub,ind_run);
                        grandmean= grandmean + data(ind_vox,sub,ind_run+ind_sec);
                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        if sess == 1
                            temp_ind = ind_run;
                        else
                            temp_ind = ind_run+ind_sec;
                        end;
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + data(ind_vox,sub,temp_ind);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                            
                            subjmean(sub)=subjmean(sub) + data(ind_vox,sub,temp_ind);
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
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                             
                            WMS = WMS + (data(ind_vox,sub,temp_ind)-subjmean(sub)).^2;
                            EMS = EMS + (data(ind_vox,sub,temp_ind)-subjmean(sub)-sessionmean(temp_ind)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                     
                        JMS=  JMS + (sessionmean(temp_ind)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
                    ICC_con_ROI(ind_vox,1) = voxICC_con;
                    z_ICC_con_ROI(ind_vox,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI(ind_vox,1) = voxICC_abs;
                    z_ICC_abs_ROI(ind_vox,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;                    
                 end;
                 
         disp('...saving ICC images...')
         ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
         ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
         z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
         z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
         if roi == 1
             ICC_con_ROI(~r_roi_ind)=0;
             z_ICC_con_ROI(~r_roi_ind)=0;
             ICC_abs_ROI(~r_roi_ind)=0;
             z_ICC_abs_ROI(~r_roi_ind)=0;
         end;     

         % save ICC maps
         target_img = temp_img;
         file1 = sprintf('ICC%s_con_%d_%d_con%d.nii',str,ind_run,ind_run+ind_sec,con_count);
         target_img.fileprefix = file1;
         target_img.img = ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
         target_img = temp_img;
         file2 = sprintf('ICC%s_abs_%d_%d_con%d.nii',str,ind_run,ind_run+ind_sec,con_count);
         target_img.fileprefix = file2;
         target_img.img = ICC_abs_ROI;
         save_nii(target_img,target_img.fileprefix); 

         target_img = temp_img;
         file3 = sprintf('z_ICC%s_con_%d_%d_con%d.nii',str,ind_run,ind_run+ind_sec,con_count);
         target_img.fileprefix = file3;
         target_img.img = z_ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
        target_img = temp_img;
        file4 = sprintf('z_ICC%s_abs_%d_%d_con%d.nii',str,ind_run,ind_run+ind_sec,con_count);
        target_img.fileprefix = file4;
        target_img.img = z_ICC_abs_ROI;
        save_nii(target_img,target_img.fileprefix);   
                             
        % computing average ICC
         disp('...computing average ICC...')
            %compute z mean and inverse Fisher's transformation
            if abs == 1
                mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCabs%s_%d_%d_con%d',str,ind_run,ind_run+ind_sec,con_count);
            end;
            if cons == 1
                mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCcon%s_%d_%d_con%d',str,ind_run,ind_run+ind_sec,con_count);
            end;
        clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
             end;
        end;
    end;
    end;
        
% if nr_para > 0
%      fprintf('...ICC parametric modulator for contrast %d between session %d and %d...\n',ind_con,ind_run, ind_run+ind_sec)
%      for ind_para = 1:nr_para
% %initiating variables
%     ICC_con_ROI=zeros(nr_vox,1);
%     ICC_abs_ROI=zeros(nr_vox,1);
%     z_ICC_con_ROI=zeros(nr_vox,1);
%     z_ICC_abs_ROI=zeros(nr_vox,1);
%    
%     %create data matrix
%     data=zeros(nr_vox,nr_subj,runs);  
%     for ind_runs = 1:runs
%          for ind_subj = 1:nr_subj
%              eval(sprintf('temp = img_con%d_par%d_%d(:,:,:,ind_subj);',ind_con,ind_para,ind_runs));
%              data(:,ind_subj,ind_runs)=temp(:);
%          end;
%     end;
%     clear temp 
%     
%     for ind_run = 1:runs
%         for ind_sec = 1:runs-1
%              if ind_run + ind_sec <= runs
%                 for ind_vox = 1:nr_vox
%                     nsamples=nr_subj*2;
% 
%                     grandmean=0;
%                     for sub=1:nr_subj,     
%                         grandmean= grandmean + data(ind_vox,sub,ind_run);
%                         grandmean= grandmean + data(ind_vox,sub,ind_run+ind_sec);
%                     end;
%                     grandmean=grandmean./nsamples;
% 
%                     sessionmean=zeros(2,1);
%                     for sess=1:2
%                         if sess == 1
%                             temp_ind = ind_run;
%                         else
%                             temp_ind = ind_run+ind_sec;
%                         end;
%                         for sub=1:nr_subj,  
%                             sessionmean(sess) = sessionmean(sess) + data(ind_vox,sub,temp_ind);
%                         end
%                         sessionmean(sess)=sessionmean(sess)./nr_subj;
%                     end
% 
%                     subjmean=zeros(nr_subj,1);
%                     for sub=1:nr_subj
%                         for sess=1:2
%                             if sess == 1
%                                 temp_ind = ind_run;
%                             else
%                                 temp_ind = ind_run+ind_sec;
%                             end;                            
%                             subjmean(sub)=subjmean(sub) + data(ind_vox,sub,temp_ind);
%                         end
%                           subjmean(sub)=subjmean(sub)./2;
%                     end
% 
%                     % mean squares
%                     BMS=0; % between subject
%                     WMS=0; % within subject 
%                     EMS=0; % error
%                     JMS=0; % session
% 
%                     for sub=1:nr_subj,    
%                         BMS = BMS + (subjmean(sub)-grandmean).^2;
%                         for sess=1:2
%                             if sess == 1
%                                 temp_ind = ind_run;
%                             else
%                                 temp_ind = ind_run+ind_sec;
%                             end;                             
%                             WMS = WMS + (data(ind_vox,sub,temp_ind)-subjmean(sub)).^2;
%                             EMS = EMS + (data(ind_vox,sub,temp_ind)-subjmean(sub)-sessionmean(temp_ind)+grandmean).^2;
%                         end
%                     end;
% 
%                     for sess=1:2
%                             if sess == 1
%                                 temp_ind = ind_run;
%                             else
%                                 temp_ind = ind_run+ind_sec;
%                             end;                     
%                         JMS=  JMS + (sessionmean(temp_ind)-grandmean).^2;
%                     end;
% 
%                     %define the true value of the mean square.
%                     BMS= 2.*BMS./(nr_subj-1);
%                     WMS= WMS./(2-1)./nr_subj;
%                     JMS= nr_subj.*JMS./(2-1);
%                     EMS= EMS./(2-1)./(nr_subj-1); 
% 
%                     %consistency agreement  
%                     if cons==1
%                     voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
%                     ICC_con_ROI(ind_vox,1) = voxICC_con;
%                     z_ICC_con_ROI(ind_vox,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
%                     end;
% 
%                     %absolute agreement 
%                     if abs==1
%                     voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
%                                                        2.* (JMS-EMS)./nr_subj);
%                     ICC_abs_ROI(ind_vox,1) = voxICC_abs;
%                     z_ICC_abs_ROI(ind_vox,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
%                     end;                    
%                  end;
%                  
%          disp('...saving ICC images...')
%          ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
%          ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
%          z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
%          z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
%          if roi == 1
%              ICC_con_ROI(~r_roi_ind)=0;
%              z_ICC_con_ROI(~r_roi_ind)=0;
%              ICC_abs_ROI(~r_roi_ind)=0;
%              z_ICC_abs_ROI(~r_roi_ind)=0;
%          end;     
% 
%          % save ICC maps
%          target_img = temp_img;
%          file1 = sprintf('ICC%s_con_%d_%d_par%d_con%d.nii',str,ind_run,ind_run+ind_sec,ind_para,con_count);
%          target_img.fileprefix = file1;
%          target_img.img = ICC_con_ROI;
%          save_nii(target_img,target_img.fileprefix); 
% 
%          clear target_img;
%          target_img = temp_img;
%          file2 = sprintf('ICC%s_abs_%d_%d_par%d_con%d.nii',str,ind_run,ind_run+ind_sec,ind_para,con_count);
%          target_img.fileprefix = file2;
%          target_img.img = ICC_abs_ROI;
%          save_nii(target_img,target_img.fileprefix); 
% 
%          target_img = temp_img;
%          file3 = sprintf('z_ICC%s_con_%d_%d_par%d_con%d.nii',str,ind_run,ind_run+ind_sec,ind_para,con_count);
%          target_img.fileprefix = file3;
%          target_img.img = z_ICC_con_ROI;
%          save_nii(target_img,target_img.fileprefix); 
% 
%          clear target_img;
%         target_img = temp_img;
%         file4 = sprintf('z_ICC%s_abs_%d_%d_par%d_con%d.nii',str,ind_run,ind_run+ind_sec,ind_para,con_count);
%         target_img.fileprefix = file4;
%         target_img.img = z_ICC_abs_ROI;
%         save_nii(target_img,target_img.fileprefix);   
%                              
%         % computing average ICC
%          disp('...computing average ICC...')
%             %compute z mean and inverse Fisher's transformation
%             if abs == 1
%                 mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
%                 mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%                 summary(1,end+1)=mean_r;
%                 cols{1,end+1} = sprintf('mean_ICCabs%s_%d_%d_par%d_con%d',str,ind_run,ind_run+ind_sec,ind_para,con_count);
%             end;
%             if cons == 1
%                 mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
%                 mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
%                 summary(1,end+1)=mean_r;
%                 cols{1,end+1} = sprintf('mean_ICCcon%s_%d_%d_par%d_con%d',str,ind_run,ind_run+ind_sec,ind_para,con_count);
%             end;
%         clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
%              end;
%         end;
%     end;
% end;
% end;
end;
    assignin('base','summary',summary);
    assignin('base','cols',cols);

    results_ICC=dataset({summary(1,:),cols{:}});
    assignin('base','results_ICC',results_ICC);

    cd(dir_results);             
    eval(sprintf('save results_ICC%s_2cons.mat results_ICC',str));    
end;
    
    


end;
if ex4D == 1
        if nr_cond == 1
        %create data matrix
        data=zeros(nr_vox,nr_subj,runs);  
        for ind_runs = 1:runs
             for ind_subj = 1:nr_subj
                 eval(sprintf('temp = FourD%d(:,:,:,ind_subj);',ind_runs));
                 data(:,ind_subj,ind_runs)=temp(:);
             end;
        end;
        clear temp
        disp('...over all sessions...')
        %calculate ICCs
        for ind_voxel = 1:nr_vox
            nsamples=nr_subj*runs;
 
            grandmean=0;
            for sub=1:nr_subj,     
                for sess=1:runs,
                    grandmean= grandmean + data(ind_voxel,sub,sess);
                end
            end;
            grandmean=grandmean./nsamples;
 
            sessionmean=zeros(runs,1);
            for sess=1:runs
                for sub=1:nr_subj,  
                    sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
                end
                sessionmean(sess)=sessionmean(sess)./nr_subj;
            end
 
            subjmean=zeros(nr_subj,1);
            for sub=1:nr_subj
                for sess=1:runs
                    subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                    WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                    EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                voxICC_con=(BMS-EMS)./(BMS+(runs-1).*EMS); 
                ICC_con_ROI(ind_voxel,1) = voxICC_con;
                z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
            end;
                     
           %absolute agreement 
           if abs==1
               voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
               runs.* (JMS-EMS)./nr_subj);
               ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
               z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
          end;
       end;

 disp('...saving ICC image over all sessions...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s.nii',str);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s.nii',str);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s.nii',str);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s.nii',str);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 disp('...computing average ICC over all sessions...')
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s',str);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s',str);    
    end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj

%calculation ICC for each comparison
 if runs > 1
    %initiating variables
    ICC_con_ROI=zeros(nr_vox,1);
    ICC_abs_ROI=zeros(nr_vox,1);
    z_ICC_con_ROI=zeros(nr_vox,1);
    z_ICC_abs_ROI=zeros(nr_vox,1);
   
    %create data matrix
    data1=zeros(nr_vox,nr_subj,runs);  
    for ind_runs = 1:runs
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = FourD%d(:,:,:,ind_subj);',ind_runs));
             data1(:,ind_subj,ind_runs)=temp(:);
         end;
    end;
    clear temp 
    
    for ind_run = 1:runs
        for ind_sec = 1:runs-1
             if ind_run + ind_sec <= runs
                 fprintf('...calculate ICC for session %d and session %d...\n',ind_run,ind_run+ind_sec);
                for ind_vox = 1:nr_vox
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        grandmean= grandmean + data1(ind_vox,sub,ind_run);
                        grandmean= grandmean + data1(ind_vox,sub,ind_run+ind_sec);
                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        if sess == 1
                            temp_ind = ind_run;
                        else
                            temp_ind = ind_run+ind_sec;
                        end;
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + data1(ind_vox,sub,temp_ind);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                            
                            subjmean(sub)=subjmean(sub) + data1(ind_vox,sub,temp_ind);
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
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                             
                            WMS = WMS + (data1(ind_vox,sub,temp_ind)-subjmean(sub)).^2;
                            EMS = EMS + (data1(ind_vox,sub,temp_ind)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                     
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
                    ICC_con_ROI(ind_vox,1) = voxICC_con;
                    z_ICC_con_ROI(ind_vox,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI(ind_vox,1) = voxICC_abs;
                    z_ICC_abs_ROI(ind_vox,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;                    
                 end;
                 
         disp('...saving ICC images...')
         ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
         ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
         z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
         z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
         if roi == 1
             ICC_con_ROI(~r_roi_ind)=0;
             z_ICC_con_ROI(~r_roi_ind)=0;
             ICC_abs_ROI(~r_roi_ind)=0;
             z_ICC_abs_ROI(~r_roi_ind)=0;
         end;     

         % save ICC maps
         target_img = temp_img;
         file1 = sprintf('ICC%s_con_%d_%d.nii',str,ind_run,ind_run+ind_sec);
         target_img.fileprefix = file1;
         target_img.img = ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
         target_img = temp_img;
         file2 = sprintf('ICC%s_abs_%d_%d.nii',str,ind_run,ind_run+ind_sec);
         target_img.fileprefix = file2;
         target_img.img = ICC_abs_ROI;
         save_nii(target_img,target_img.fileprefix); 

         target_img = temp_img;
         file3 = sprintf('z_ICC%s_con_%d_%d.nii',str,ind_run,ind_run+ind_sec);
         target_img.fileprefix = file3;
         target_img.img = z_ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
        target_img = temp_img;
        file4 = sprintf('z_ICC%s_abs_%d_%d.nii',str,ind_run,ind_run+ind_sec);
        target_img.fileprefix = file4;
        target_img.img = z_ICC_abs_ROI;
        save_nii(target_img,target_img.fileprefix);   
                             
        % computing average ICC
                 fprintf('...calculate average ICC for session %d and session %d...\n',ind_run,ind_run+ind_sec);
            %compute z mean and inverse Fisher's transformation
            if abs == 1
                mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCabs%s_%d_%d',str,ind_run,ind_run+ind_sec);
            end;
            if cons == 1
                mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCcon%s_%d_%d',str,ind_run,ind_run+ind_sec);    
            end;
            clear ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj

             end;
        end;
    end;
clear data1;
        
end;
        elseif nr_cond > 1
             %create data matrix
    for ind_runs = 1:runs
        data=zeros(nr_vox,nr_subj,2);  
        for ind_con = 1:2
             for ind_subj = 1:nr_subj
                 eval(sprintf('temp = FourD%s%d(:,:,:,ind_subj);',conditions{ind_con,1},ind_con));
                 data(:,ind_subj,ind_con)=temp(:);
             end;
        end;
    clear temp
        
    %calculate ICCs
    fprintf('...calculate ICC betweens contrasts for session %d...\n ',ind_runs);

    for ind_voxel = 1:nr_vox
        nsamples=nr_subj*2;

        grandmean=0;
        for sub=1:nr_subj,     
            for sess=1:2,
                grandmean= grandmean + data(ind_voxel,sub,sess);
            end
        end;
        grandmean=grandmean./nsamples;

        sessionmean=zeros(2,1);
        for sess=1:2
            for sub=1:nr_subj,  
                sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
            end
            sessionmean(sess)=sessionmean(sess)./nr_subj;
        end

        subjmean=zeros(nr_subj,1);
        for sub=1:nr_subj
            for sess=1:2
                subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
            voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
            ICC_con_ROI(ind_voxel,1) = voxICC_con;
            z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
        end;

       %absolute agreement 
       if abs==1
           voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
           2.* (JMS-EMS)./nr_subj);
           ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
           z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
      end;
   end;

 disp('...saving ICC images...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s_%d_between_conditions.nii',str,ind_runs);
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s_%d_between_conditions.nii',str,ind_runs);
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s_%d_between_conditions.nii',str,ind_runs);
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s_%d_between_conditions.nii',str,ind_runs);
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 fprintf('...computing average ICC between conditions in session %d...\n',ind_runs)
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s_%d_between_conditions',str,ind_runs);
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s_%d_between_conditions',str,ind_runs);    
    end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
 
    end;
    
 % calculate voxelwise average ICC over all conditions
 disp('...creates average ICC maps over all condition comparisons...');   
    avg_ICC_4D_abs = zeros(nr_vox,runs);
    avg_ICC_4D_con = zeros(nr_vox,runs);
    for i_run = 1:runs
        if abs == 1
        map = load_nii(sprintf('z_ICC_abs%s_%d_between_conditions.nii',str,i_run));
        avg_ICC_4D_abs(:,i_run) = map.img(:);
        clear map
        end;
        if cons == 1
        map = load_nii(sprintf('z_ICC_con%s_%d_between_conditions.nii',str,i_run));
        avg_ICC_4D_con(:,i_run) = map.img(:);
        clear map
        end;                
    end;
    avg_abs = zeros(nr_vox,1);
    avg_con = zeros(nr_vox,1);
    if abs == 1
        for ind_vox = 1:nr_vox
            avg_temp = mean(avg_ICC_4D_abs(ind_vox,:));
            avg_abs(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;
    if cons == 1
        for ind_vox = 1:nr_vox      
            avg_temp = mean(avg_ICC_4D_con(ind_vox,:));
            avg_con(ind_vox,1) = (exp(2*avg_temp)-1)./(exp(2*avg_temp)+1);
        end;
    end;                
    clear avg_ICC_4D_abs avg_ICC_4D_con avg_temp
    if abs == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_abs%s_between_conditions.nii',str);
        target_img.img = avg_abs;
        save_nii(target_img,target_img.fileprefix); 
    end;
    if cons == 1
        target_img = temp_img;
        target_img.fileprefix = sprintf('ICC_con%s_between_conditions.nii',str);
        target_img.img = avg_con;
        save_nii(target_img,target_img.fileprefix);        
    end;
    clear avg_abs avg_con 
    
    % ICCs within conditions over sessions
for ind_con = 1:2
    fprintf('...calculating ICC within condition %s over all sessions...\n',conditions{ind_con,1})
    count_comp = 0;  
%initiating variables
    ICC_con_ROI=zeros(nr_subj,1);
    ICC_abs_ROI=zeros(nr_subj,1);
    z_ICC_con_ROI=zeros(nr_subj,1);
    z_ICC_abs_ROI=zeros(nr_subj,1);

    disp('...calculating ICCs...')
    
           %create data matrix
            data=zeros(nr_vox,nr_subj,runs);  
        for ind_runs = 1:runs
             for ind_subj = 1:nr_subj
                 eval(sprintf('temp = FourD%s%d(:,:,:,ind_subj);',conditions{ind_con,1},ind_runs));
                 data(:,ind_subj,ind_runs)=temp(:);
             end;
        end;
        clear temp
        
        %calculate ICCs
        for ind_voxel = 1:nr_vox
            nsamples=nr_subj*runs;
 
            grandmean=0;
            for sub=1:nr_subj,     
                for sess=1:runs,
                    grandmean= grandmean + data(ind_voxel,sub,sess);
                end
            end;
            grandmean=grandmean./nsamples;
 
            sessionmean=zeros(runs,1);
            for sess=1:runs
                for sub=1:nr_subj,  
                    sessionmean(sess) = sessionmean(sess) + data(ind_voxel,sub,sess);
                end
                sessionmean(sess)=sessionmean(sess)./nr_subj;
            end
 
            subjmean=zeros(nr_subj,1);
            for sub=1:nr_subj
                for sess=1:runs
                    subjmean(sub)=subjmean(sub) + data(ind_voxel,sub,sess);
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
                    WMS = WMS + (data(ind_voxel,sub,sess)-subjmean(sub)).^2;
                    EMS = EMS + (data(ind_voxel,sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                voxICC_con=(BMS-EMS)./(BMS+(runs-1).*EMS); 
                ICC_con_ROI(ind_voxel,1) = voxICC_con;
                z_ICC_con_ROI(ind_voxel,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
            end;
                     
           %absolute agreement 
           if abs==1
               voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
               runs.* (JMS-EMS)./nr_subj);
               ICC_abs_ROI(ind_voxel,1) = voxICC_abs;
               z_ICC_abs_ROI(ind_voxel,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
          end;
       end;

 disp('...saving ICC images...')
 ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
 ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
 z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
 z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
 if roi == 1
     ICC_con_ROI(~r_roi_ind)=0;
     z_ICC_con_ROI(~r_roi_ind)=0;
     ICC_abs_ROI(~r_roi_ind)=0;
     z_ICC_abs_ROI(~r_roi_ind)=0;
 end;
 
 % save ICC maps
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_con%s_%s.nii',str,conditions{ind_con,1});
 target_img.img = ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('ICC_abs%s_%s.nii',str,conditions{ind_con,1});
 target_img.img = ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_con%s_%s.nii',str,conditions{ind_con,1});
 target_img.img = z_ICC_con_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
 clear target_img;
 target_img = temp_img;
 target_img.fileprefix = sprintf('z_ICC_abs%s_%s.nii',str,conditions{ind_con,1});
 target_img.img = z_ICC_abs_ROI;
 save_nii(target_img,target_img.fileprefix); 
 
% computing average ICC
 disp('...computing average ICC...')
    %compute z mean and inverse Fisher's transformation
    if abs == 1
        mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCabs%s_%s',str,conditions{ind_con,1});
    end;
    if cons == 1
        mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
        mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
        summary(1,end+1)=mean_r;
        cols{1,end+1} = sprintf('mean_ICCcon%s_%s',str,conditions{ind_con,1});    
    end;
end;
clearvars data ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj


%calculation ICC for each comparison
 if runs > 1
    for ind_cons = 1:2
    %initiating variables
    ICC_con_ROI=zeros(nr_vox,1);
    ICC_abs_ROI=zeros(nr_vox,1);
    z_ICC_con_ROI=zeros(nr_vox,1);
    z_ICC_abs_ROI=zeros(nr_vox,1);
   
    %create data matrix
    data=zeros(nr_vox,nr_subj,runs);  
    for ind_runs = 1:runs
         for ind_subj = 1:nr_subj
             eval(sprintf('temp = FourD%s%d(:,:,:,ind_subj);',conditions{ind_cons,1},ind_runs));
             data(:,ind_subj,ind_runs)=temp(:);
         end;
    end;
    clear temp 
    
    for ind_run = 1:runs
        for ind_sec = 1:runs-1
             if ind_run + ind_sec <= runs
    fprintf('...calculating ICC within condition %s between session %d and %d...',conditions{ind_cons,1},ind_run, ind_run+ind_sec)
                 
                for ind_vox = 1:nr_vox
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        grandmean= grandmean + data(ind_vox,sub,ind_run);
                        grandmean= grandmean + data(ind_vox,sub,ind_run+ind_sec);
                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        if sess == 1
                            temp_ind = ind_run;
                        else
                            temp_ind = ind_run+ind_sec;
                        end;
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + data(ind_vox,sub,temp_ind);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                            
                            subjmean(sub)=subjmean(sub) + data(ind_vox,sub,temp_ind);
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
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                             
                            WMS = WMS + (data(ind_vox,sub,temp_ind)-subjmean(sub)).^2;
                            EMS = EMS + (data(ind_vox,sub,temp_ind)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:2
                            if sess == 1
                                temp_ind = ind_run;
                            else
                                temp_ind = ind_run+ind_sec;
                            end;                     
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(2-1).*EMS); 
                    ICC_con_ROI(ind_vox,1) = voxICC_con;
                    z_ICC_con_ROI(ind_vox,1) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI(ind_vox,1) = voxICC_abs;
                    z_ICC_abs_ROI(ind_vox,1) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;                    
                 end;
                 
         disp('...saving ICC images...')
         ICC_con_ROI = reshape(ICC_con_ROI,x,y,z);
         ICC_abs_ROI = reshape(ICC_abs_ROI,x,y,z);
         z_ICC_con_ROI = reshape(z_ICC_con_ROI,x,y,z);
         z_ICC_abs_ROI = reshape(z_ICC_abs_ROI,x,y,z);
         if roi == 1
             ICC_con_ROI(~r_roi_ind)=0;
             z_ICC_con_ROI(~r_roi_ind)=0;
             ICC_abs_ROI(~r_roi_ind)=0;
             z_ICC_abs_ROI(~r_roi_ind)=0;
         end;     

         % save ICC maps
         target_img = temp_img;
         file1 = sprintf('ICC%s_con_%d_%d_%s.nii',str,ind_run,ind_run+ind_sec,conditions{ind_cons,1});
         target_img.fileprefix = file1;
         target_img.img = ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
         target_img = temp_img;
         file2 = sprintf('ICC%s_abs_%d_%d_%s.nii',str,ind_run,ind_run+ind_sec,conditions{ind_cons,1});
         target_img.fileprefix = file2;
         target_img.img = ICC_abs_ROI;
         save_nii(target_img,target_img.fileprefix); 

         target_img = temp_img;
         file3 = sprintf('z_ICC%s_con_%d_%d_%s.nii',str,ind_run,ind_run+ind_sec,conditions{ind_cons,1});
         target_img.fileprefix = file3;
         target_img.img = z_ICC_con_ROI;
         save_nii(target_img,target_img.fileprefix); 

         clear target_img;
        target_img = temp_img;
        file4 = sprintf('z_ICC%s_abs_%d_%d_%s.nii',str,ind_run,ind_run+ind_sec,conditions{ind_cons,1});
        target_img.fileprefix = file4;
        target_img.img = z_ICC_abs_ROI;
        save_nii(target_img,target_img.fileprefix);   
                             
        % computing average ICC
         disp('...computing average ICC...')
            %compute z mean and inverse Fisher's transformation
            if abs == 1
                mean_z = mean(z_ICC_abs_ROI(~isinf(z_ICC_abs_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCabs%s_%d_%d_%s',str,ind_run,ind_run+ind_sec,conditions{ind_cons,1});
            end;
            if cons == 1
                mean_z = mean(z_ICC_con_ROI(~isinf(z_ICC_con_ROI)),'omitnan');
                mean_r=(exp(2*mean_z)-1)./(exp(2*mean_z)+1);
                summary(1,end+1)=mean_r;
                cols{1,end+1} = sprintf('mean_ICCcon%s_%d_%d_%s',str,ind_run,ind_run+ind_sec,conditions{ind_cons,1});
            end;
        clearvars ICC_con_ROI ICC_abs_ROI z_ICC_con_ROI z_ICC_abs_ROI ROI_subj
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
    eval(sprintf('save results_ICC%s.mat results_ICC',str));
end;
end;
disp('finished ICC calculation')


cd(box_path);    
disp('DONE');
 

% % --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to axes1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: place code in OpeningFcn to populate axes1
axes(hObject);
imshow('logo.png');
