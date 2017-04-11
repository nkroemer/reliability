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

% Last Modified by GUIDE v2.5 10-Apr-2017 15:51:13

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
disp('starting similarity calculations');
study_design=evalin('base','study_design');
contrast_def = evalin('base','contrast_def');

results_dir = study_design.results_directory;
runs = str2double(study_design.number_sessions);
nr_subj = str2double(study_design.number_subjects);
load(study_design.subject_list);
stats=study_design.stats_directory;
path=study_design.stats_path;
box_path = pwd;
split = get(handles.split,'Value');

%load contrast information
two_cons = contrast_def.two_contrasts;
if two_cons == 0
    con=contrast_def.contrast;
else
    con1=contrast_def.contrast1;
    con2=contrast_def.contrast2;
    con1_count=contrast_def.number_contrast1;
    con2_count=contrast_def.number_contrast2;
end;

% load and reslice ROI
use_roi = get(handles.use_roi,'Value');
if use_roi == 1
    disp('...load and reslice ROI...')
    roi_name=get(handles.name_roi,'String');
    roi_dir=evalin('base','roi_dir');
    cd(roi_dir);
    roi_name = dir(sprintf('%s*',roi_name));
    roi_name = roi_name.name;
    roi_compl = sprintf('%s\\%s',roi_dir,roi_name);
    disp('...reslicing ROI...');
    stats_filled = sprintf(stats,1);
    if two_cons == 1
        temp = sprintf('%s\\%s\\%s\\%s,1',path,vp{1},stats_filled,con1);
    else
        temp = sprintf('%s\\%s\\%s\\%s,1',path,vp{1},stats_filled,con);
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
    r_roi = sprintf('%s\\r%s',roi_dir,roi_name);
    movefile (r_roi,results_dir,'f');

    % create index for ROI voxels
    cd(results_dir)
    r_roi = load_nii(sprintf('r%s',roi_name));
    r_roi_ind = r_roi.img==1;
end;

% load 4D images
if two_cons == 0 && split == 0
    for ind_run = 1:runs
        fprintf('...load 4D image for run %d...\n',ind_run);
        file = sprintf('%s\\4D_%d.nii',results_dir,ind_run);
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
elseif two_cons == 1
     for ind_run = 1:runs
        fprintf('...load 4D image for run %d and contrast %s...\n',ind_run,con1);
        file1 = sprintf('%s\\4D_%s_%d.nii',results_dir,con1,ind_run);
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

        fprintf('...load 4D image for run %d and contrast %s...',ind_run,con2);
        file2 = sprintf('%s\\4D_%s_%d.nii',results_dir,con2,ind_run);
        eval(sprintf('FourD%d = file2;',ind_run));
        eval(sprintf('FourD%d = load_nii(FourD%d);',ind_run,ind_run));
        eval(sprintf('FourD%d_%d = FourD%d.img;',ind_run,con2_count,ind_run));
        eval(sprintf('dims = size(FourD%d_%d(:,:,:,1));',ind_run,con2_count));
        x = dims(1);
        y = dims(2);
        z = dims(3);
        fprintf('...generate mean activation per voxel for run %d and contrast %s...',ind_run,con2); 
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
elseif split == 1
for ind_run = 1:runs
        fprintf('...load 4D image for run %d and split 1...\n',ind_run);
        file1 = sprintf('%s\\4D_split1_par1_%d.nii',results_dir,ind_run);
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
        file2 = sprintf('%s\\4D_split2_par1_%d.nii',results_dir,ind_run);
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
end;

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
            eval(sprintf('save similarity_%s_%d-%d.mat out',roi_name,ind_run, ind_run+count));           
        else
            eval(sprintf('save similarity_%d-%d.mat out',ind_run, ind_run+count));
        end;
        
        % color matrix    
        f1 = figure;
        colormap('jet');
        imagesc(out.r_mat);
        colorbar;
        hold off;
        if use_roi == 1
            name=sprintf('similarity_%s_%d-%d',roi_name,ind_run, ind_run+count);
        else
            name=sprintf('similarity_%d-%d',ind_run, ind_run+count);
        end;
        print(f1,name,'-dpng');

        % histograms
        triu_r_mat = triu(out.r_mat,1);
        vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);

        f2 = figure;
        h1 = histogram(diag(out.r_mat));
        hold on;
        h2 = histogram(vec_off_diag_r_mat);
        h1.Normalization = 'probability';
        h1.BinWidth = 0.1;
        h2.Normalization = 'probability';
        h2.BinWidth = 0.1;
        hold off;
        if use_roi == 1
            name=sprintf('histograms_%s_%d-%d',roi_name,ind_run, ind_run+count);
        else
            name=sprintf('histograms_%d-%d',ind_run, ind_run+count);
        end;
        print(f2,name,'-dpng');
                
        % ecdf - densitiy plots
        f3 = figure;
        ecdf(diag(out.r_mat),'bounds','on');
        hold on;
        ecdf(vec_off_diag_r_mat,'bounds','on');
        hold off;
        if use_roi == 1
            name=sprintf('density_%s_%d-%d',roi_name,ind_run, ind_run+count);
        else
            name=sprintf('density_%d-%d',ind_run, ind_run+count);
        end;            
        print(f3,name,'-dpng');
        
        close(f1);
        close(f2);
        close(f3);
        clearvars out f1 f2 f3;
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
            colormap('jet');
            imagesc(out.r_mat);
            colorbar;
            hold off;
            if use_roi == 1
                name=sprintf('similarity_%s_%d-%d_con%d',roi_name,ind_run, ind_run+count,con1_count);
            else
                name=sprintf('similarity_%d-%d_con%d',ind_run, ind_run+count,con1_count);
            end;
            print(f1,name,'-dpng');

            % histograms
            triu_r_mat = triu(out.r_mat,1);
            vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);

            f2 = figure;
            h1 = histogram(diag(out.r_mat));
            hold on;
            h2 = histogram(vec_off_diag_r_mat);
            h1.Normalization = 'probability';
            h1.BinWidth = 0.1;
            h2.Normalization = 'probability';
            h2.BinWidth = 0.1;
            hold off;
            if use_roi == 1
                name=sprintf('histograms_%s_%d-%d_con%d',roi_name,ind_run, ind_run+count,con1_count);
            else
                name=sprintf('histograms_%d-%d_con%d',ind_run, ind_run+count,con1_count);
            end;
            print(f2,name,'-dpng');

            % ecdf - densitiy plots
            f3 = figure;
            ecdf(diag(out.r_mat),'bounds','on');
            hold on;
            ecdf(vec_off_diag_r_mat,'bounds','on');
            hold off;
            if use_roi == 1
                name=sprintf('density_%s_%d-%d_con%d',roi_name,ind_run, ind_run+count,con1_count);
            else
                name=sprintf('density_%d-%d_con%d',ind_run, ind_run+count,con1_count);
            end;            
            print(f3,name,'-dpng');
        close(f1);
        close(f2);
        close(f3);
            clearvars out f1 f2 f3;
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
                eval(sprintf('save similarity_%s_%d-%d_con%d.mat out',roi_name,ind_run, ind_run+count,con2_count));           
            else
                eval(sprintf('save similarity_%d-%d_con%d.mat out',ind_run, ind_run+count,con2_count));
            end;

            % color matrix    
            f1 = figure;
            colormap('jet');
            imagesc(out.r_mat);
            colorbar;
            hold off;
            if use_roi == 1
                name=sprintf('similarity_%s_%d-%d_con%d',roi_name,ind_run, ind_run+count,con2_count);
            else
                name=sprintf('similarity_%d-%d_con%d',ind_run, ind_run+count,con2_count);
            end;
            print(f1,name,'-dpng');

            % histograms
            triu_r_mat = triu(out.r_mat,1);
            vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);

            f2 = figure;
            h1 = histogram(diag(out.r_mat));
            hold on;
            h2 = histogram(vec_off_diag_r_mat);
            h1.Normalization = 'probability';
            h1.BinWidth = 0.1;
            h2.Normalization = 'probability';
            h2.BinWidth = 0.1;
            hold off;
            if use_roi == 1
                name=sprintf('histograms_%s_%d-%d_con%d',roi_name,ind_run, ind_run+count,con2_count);
            else
                name=sprintf('histograms_%d-%d_con%d',ind_run, ind_run+count,con2_count);
            end;
            print(f2,name,'-dpng');

            % ecdf - densitiy plots
            f3 = figure;
            ecdf(diag(out.r_mat),'bounds','on');
            hold on;
            ecdf(vec_off_diag_r_mat,'bounds','on');
            hold off;
            if use_roi == 1
                name=sprintf('density_%s_%d-%d_con%d',roi_name,ind_run, ind_run+count,con2_count);
            else
                name=sprintf('density_%d-%d_con%d',ind_run, ind_run+count,con2_count);
            end;            
            print(f3,name,'-dpng');
         close(f1);
        close(f2);
        close(f3);           
            clearvars out f1 f2 f3;
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
                eval(sprintf('save similarity_%s_%d_con%d_con%d.mat out',roi_name,i_run,con1_count,con2_count));           
            else
                eval(sprintf('save similarity_%d_con%d_con%d.mat out',i_run,con1_count,con2_count));
            end;

            % color matrix    
            f1 = figure;
            colormap('jet');
            imagesc(out.r_mat);
            colorbar;
            hold off;
            if use_roi == 1
                name=sprintf('similarity_%s_%d_con%d_con%d',roi_name,i_run,con1_count,con2_count);
            else
                name=sprintf('similarity_%d_con%d_con%d',i_run,con1_count,con2_count);
            end;
            print(f1,name,'-dpng');

            % histograms
            triu_r_mat = triu(out.r_mat,1);
            vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);

            f2 = figure;
            h1 = histogram(diag(out.r_mat));
            hold on;
            h2 = histogram(vec_off_diag_r_mat);
            h1.Normalization = 'probability';
            h1.BinWidth = 0.1;
            h2.Normalization = 'probability';
            h2.BinWidth = 0.1;
            hold off;
            if use_roi == 1
                name=sprintf('histograms_%s_%d_con%d_con%d',roi_name,i_run,con1_count,con2_count);
            else
                name=sprintf('histograms_%d_con%d_con%d',i_run,con1_count,con2_count);
            end;
            print(f2,name,'-dpng');

            % ecdf - densitiy plots
            f3 = figure;
            ecdf(diag(out.r_mat),'bounds','on');
            hold on;
            ecdf(vec_off_diag_r_mat,'bounds','on');
            hold off;
            if use_roi == 1
                name=sprintf('density_%s_%d_con%d_con%d',roi_name,i_run,con1_count,con2_count);
            else
                name=sprintf('density_%d_con%d_con%d',i_run,con1_count,con2_count);
            end;            
            print(f3,name,'-dpng');
        close(f1);
        close(f2);
        close(f3);
            clearvars out f1 f2 f3;        
    end;
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
                eval(sprintf('mean_temp = mean_FourD%d_split1;',i_run)); 
               [out.r, out.p] = corrcoef([temp_2,mean_temp]); 
               out.r_mat(j,nr_subj+2) = out.r(1,2);
               out.p_mat(j,nr_subj+2) = out.p(1,2); 
            end;
               eval(sprintf('mean_temp = mean_FourD%d_split2;',i_run)); 
               [out.r, out.p] = corrcoef([temp_1,mean_temp]); 
               out.r_mat(i,j+1) = out.r(1,2);
               out.p_mat(i,j+1) = out.p(1,2);   
              
        end;

            cd(results_dir);
            roi_name=get(handles.name_roi,'String');

            if use_roi == 1
                eval(sprintf('save similarity_%s_%d_split.mat out',roi_name,i_run));           
            else
                eval(sprintf('save similarity_%d_split.mat out',i_run));
            end;

            % color matrix    
            f1 = figure;
            colormap('jet');
            imagesc(out.r_mat);
            colorbar;
            hold off;
            if use_roi == 1
                name=sprintf('similarity_%s_%d_split',roi_name,i_run);
            else
                name=sprintf('similarity_%d_split',i_run);
            end;
            print(f1,name,'-dpng');

            % histograms
            triu_r_mat = triu(out.r_mat,1);
            vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);

            f2 = figure;
            h1 = histogram(diag(out.r_mat));
            hold on;
            h2 = histogram(vec_off_diag_r_mat);
            h1.Normalization = 'probability';
            h1.BinWidth = 0.1;
            h2.Normalization = 'probability';
            h2.BinWidth = 0.1;
            hold off;
            if use_roi == 1
                name=sprintf('histograms_%s_%d_split',roi_name,i_run);
            else
                name=sprintf('histograms_%d_split',i_run);
            end;
            print(f2,name,'-dpng');

            % ecdf - densitiy plots
            f3 = figure;
            ecdf(diag(out.r_mat),'bounds','on');
            hold on;
            ecdf(vec_off_diag_r_mat,'bounds','on');
            hold off;
            if use_roi == 1
                name=sprintf('density_%s_%d_split',roi_name,i_run);
            else
                name=sprintf('density_%d_split',i_run);
            end;            
            print(f3,name,'-dpng');
        close(f1);
        close(f2);
        close(f3);
            clearvars out f1 f2 f3;        
    end;    
end;

cd(box_path);
disp('DONE');
