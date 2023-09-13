function varargout = SplitHalf(varargin)
% SPLITHALF MATLAB code for SplitHalf.fig
%      SPLITHALF, by itself, creates a new SPLITHALF or raises the existing
%      singleton*.
%
%      H = SPLITHALF returns the handle to a new SPLITHALF or the handle to
%      the existing singleton*.
%
%      SPLITHALF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPLITHALF.M with the given input arguments.
%
%      SPLITHALF('Property','Value',...) creates a new SPLITHALF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SplitHalf_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SplitHalf_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SplitHalf

% Last Modified by GUIDE v2.5 10-Oct-2017 15:49:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SplitHalf_OpeningFcn, ...
                   'gui_OutputFcn',  @SplitHalf_OutputFcn, ...
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


% --- Executes just before SplitHalf is made visible.
function SplitHalf_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SplitHalf (see VARARGIN)

% Choose default command line output for SplitHalf
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SplitHalf wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SplitHalf_OutputFcn(hObject, eventdata, handles) 
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
disp('starting creation of SplitHalf statistics');

name_study_design = cellstr(spm_select(1,'mat','load study design'));
load(name_study_design{1});
assignin('base','study_design',study_design);
assignin('base','name_study_design',name_study_design);

function nr_con_Callback(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr_subj as text
%        str2double(get(hObject,'String')) returns contents of nr_subj as a double

% --- Executes during object creation, after setting all properties.
function nr_con_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in par.
function par_Callback(hObject, eventdata, handles)
% hObject    handle to par (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of par
function split_name_Callback(hObject, eventdata, handles)
% hObject    handle to split_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of split_name as text
%        str2double(get(hObject,'String')) returns contents of split_name as a double


% --- Executes during object creation, after setting all properties.
function split_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to split_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% define file seperator 
f = filesep;
box_path=evalin('base','box_path');

%% set parameters

% get study design info
study_design=evalin('base','study_design');
name_study_design=evalin('base','name_study_design');
subjects = study_design.subject_list;
load(subjects);
nr_subj = str2double(study_design.number_subjects);
stats_path = study_design.stats_path;
runs = study_design.number_sessions;
stats_dir = study_design.stats_directory;

% get GUI input
con_vec = str2num(get(handles.nr_con,'String'));
split_name = get(handles.split_name,'String');

%% create list with SPM.mat files

SPM_list = cell(nr_subj,runs);
disp('...create list with all SPM.mat files...')
for j = 1:runs
    for i = 1:nr_subj
        cd(stats_path);
        cd(sprintf('%s',id{i}));
        cd(sprintf(stats_dir,j));

        if exist('SPM.mat','file')==2
            path=pwd;
            path=sprintf('%s%s%sSPM.mat',path,f,f);
            SPM_list{i,j}=path;
        end;
    end; 
end;

%% add study design info
% info split directory in study_design
split_dir = sprintf('%s_split_%s',stats_dir,split_name);
study_design.split_directory=split_dir;
study_design.nr_split_reg = length(con_vec);

% info number parametric modulators in study_design
load(SPM_list{1});
if length(con_vec) == 1
    con = con_vec(1,1);
    if strcmp(SPM.Sess(1).U(con).P(1).name,'none')
        study_design.number_parametric = 0;
        nr_para = 0;
    else
        nr_para = size(SPM.Sess(1).U(con).P,2);
        study_design.number_parametric = nr_para;
    end;
else
    for i_con = 1:length(con_vec)
        con = con_vec(1,i_con);
        if strcmp(SPM.Sess(1).U(con).P(1).name,'none')
            eval(sprintf('study_design.number_parametric_reg%d = 0;',i_con));
            eval(sprintf('nr_para_reg%d = 0',i_con));
        else
            eval(sprintf('nr_para_reg%d = size(SPM.Sess(1).U(con).P,2);',i_con));
            eval(sprintf('study_design.number_parametric_reg%d = nr_para_reg%d;',i_con,i_con));
        end;     
    end;
end;

name=name_study_design{1};
save(name,'study_design');


%% loads and modifies SPM for each participant

for m = 1:runs
    split_dir_temp = sprintf(split_dir,m);
    for count = 1:nr_subj
        fprintf('...load and modify SPM.mat for %s in session %d...',id{count},m)
        dir_spm = SPM_list{count,m};
        load(dir_spm);
        % create vector for session number (SPM.Sess(nr))
        sess_length = length(SPM.Sess);
        U_length = length(SPM.Sess(1).U);
        sess_nr = zeros(1,length(con_vec));
        if sess_length > 1
            for ind_contrast = 1:length(con_vec)
                con = con_vec(1,ind_contrast);
                for sess_ind = 1:sess_length
                        if con > (sess_ind.*U_length)-U_length && con <= sess_ind.*U_length
                            sess_nr(1,ind_contrast) = sess_ind;
                        end;
                end;
            end;
            for ind_con = 1:length(con_vec)

                sess_count = sess_nr(1,ind_con);
                if sess_count == 1
                    con = con_vec(1,ind_con);
                else
                    con = con_vec(1,ind_con)-((sess_count-1).*U_length);
                end;

                %split onsets
                onsets = SPM.Sess(sess_count).U(con).ons;
                rng('shuffle');
                pos = randperm(length(onsets))';
                onsets1 = onsets(pos(1:end/2));
                onsets2 = onsets(pos(end/2:end));

                %split parametric value
                if length(con_vec) == 1 && nr_para > 0
                    for k = 1:nr_para
                        eval(sprintf('mod%d = SPM.Sess(sess_count).U(con).P(k).P;',k));
                        eval(sprintf('mod%d_1 = mod%d(pos(1:end/2));',k,k));
                        eval(sprintf('mod%d_2 = mod%d(pos(end/2:end));',k,k));
                    end;
                end;
                if length(con_vec) > 1 && nr_para > 0
                    for k = 1:eval(sprintf('nr_para_reg%d',ind_con))
                        eval(sprintf('mod%d = SPM.Sess(sess_count).U(con).P(k).P;',k));
                        eval(sprintf('mod%d_1 = mod%d(pos(1:end/2));',k,k));
                        eval(sprintf('mod%d_2 = mod%d(pos(end/2:end));',k,k));
                    end;    
                end;

                %create new conditions 
                name = SPM.Sess(sess_count).U(con).name{1};
                U(1).name{1} = sprintf('half1_%s',name);
                U(1).ons = onsets1;
                U(1).dur = zeros(length(onsets1),1);
                if isfield(SPM.Sess(sess_count).U(con),'orth')
                    U(1).orth = SPM.Sess(sess_count).U(con).orth;
                end;
                
                if length(con_vec) == 1
                    if nr_para>0
                        for l = 1:nr_para
                            U(1).P(l).name = sprintf('par_value_1_%d',l);
                            eval(sprintf('U(1).P(l).P = mod%d_1;',l));
                            eval(evalstr);
                            U(1).P(l).h = 1;
                            U(1).P(l).i = [1,l+1];
                            U(2).P(l).name = sprintf('par_value_2_%d',l);
                            eval(sprintf('U(2).P(l).P = mod%d_2;',l));
                            U(2).P(l).h = 1;
                            U(2).P(l).i = [1,l+1];
                        end;
                    else
                        U(1).P.name = 'none';
                        U(1).P.h = 0;
                        U(1).P.i = 1;
                        U(2).P.name = 'none';
                        U(2).P.h = 0;
                        U(2).P.i = 1;
                    end;
                else
                     if eval(sprintf('nr_para_reg%d >0',ind_con))
                        for l = 1:nr_para
                            U(1).P(l).name = sprintf('par_value_1_%d',l);
                            eval(sprintf('U(1).P(l).P = mod%d_1;',l));
                            eval(evalstr);
                            U(1).P(l).h = 1;
                            U(1).P(l).i = [1,l+1];
                            U(2).P(l).name = sprintf('par_value_2_%d',l);
                            eval(sprintf('U(2).P(l).P = mod%d_2;',l));
                            U(2).P(l).h = 1;
                            U(2).P(l).i = [1,l+1];
                        end;
                    else
                        U(1).P.name = 'none';
                        U(1).P.h = 0;
                        U(1).P.i = 1;
                        U(2).P.name = 'none';
                        U(2).P.h = 0;
                        U(2).P.i = 1;
                    end;                   
                end;
                
                U(1).dt = SPM.Sess(sess_count).U(con).dt;
                U(1).u = SPM.Sess(sess_count).U(con).u;
                U(1).pst = SPM.Sess(sess_count).U(con).pst;    

                U(2).name{1} = sprintf('half2_%s',name);
                U(2).ons = onsets2;
                U(2).dur = zeros(length(onsets2),1);
                if isfield(SPM.Sess(sess_count).U(con),'orth')
                    U(2).orth = SPM.Sess(sess_count).U(con).orth;
                end;
                U(2).dt = SPM.Sess(sess_count).U(con).dt;
                U(2).u = SPM.Sess(sess_count).U(con).u;
                U(2).pst = SPM.Sess(sess_count).U(con).pst;

                Sess(sess_count).U(1:con-1) = SPM.Sess(sess_count).U(1:con-1);
                U_rest = SPM.Sess(sess_count).U(con+1:end);
                Sess(sess_count).U(con) = U(1);
                Sess(sess_count).U(con+1) = U(2);
                Sess(sess_count).U(end+1:(end+length(U_rest))) = U_rest;


                U = [];
            end;
            for sess_count = 1:sess_length
                SPM.Sess(sess_count).U = Sess(sess_count).U;
            end;
        else
            for ind_con = 1:length(con_vec)
                con = con_vec(1,ind_con);
                sess_count = 1;
                
                %split onsets
                onsets = SPM.Sess(sess_count).U(con).ons;
                rng('shuffle');
                pos = randperm(length(onsets))';
                onsets1 = onsets(pos(1:end/2));
                onsets2 = onsets(pos(end/2:end));

                %split parametric value
                if length(con_vec) == 1
                    if nr_para > 0
                        for k = 1:nr_para
                        eval(sprintf('mod%d = SPM.Sess(sess_count).U(con).P(k).P;',k));
                        eval(sprintf('mod%d_1 = mod%d(pos(1:end/2));',k,k));
                        eval(sprintf('mod%d_2 = mod%d(pos(end/2:end));',k,k));
                        end;
                    end;
                else
                    if eval(sprintf('nr_para_reg%d > 0',ind_con))
                        for k = 1:nr_para
                            eval(sprintf('mod%d = SPM.Sess(sess_count).U(con).P(k).P;',k));
                            eval(sprintf('mod%d_1 = mod%d(pos(1:end/2));',k,k));
                            eval(sprintf('mod%d_2 = mod%d(pos(end/2:end));',k,k));
                        end;
                    end;                   
                end;
                
                %create new conditions 
                name = SPM.Sess(sess_count).U(con).name{1};
                U(1).name{1} = sprintf('half1_%s',name);
                U(1).ons = onsets1;
                U(1).dur = zeros(length(onsets1),1);
                if isfield(SPM.Sess(sess_count).U(con),'orth')
                    U(1).orth = SPM.Sess(sess_count).U(con).orth;
                end; 
                
                if length(con_vec) == 1
                    if nr_para>0
                        for l = 1:nr_para
                            U(1).P(l).name = sprintf('par_value_1_%d',l);
                            eval(sprintf('U(1).P(l).P = mod%d_1;',l));
                            U(1).P(l).h = 1;
                            U(1).P(l).i = [1,l+1];
                            U(2).P(l).name = sprintf('par_value_2_%d',l);
                            eval(sprintf('U(2).P(l).P = mod%d_2;',l));
                            U(2).P(l).h = 1;
                            U(2).P(l).i = [1,l+1];
                        end;
                    else
                        U(1).P.name = 'none';
                        U(1).P.h = 0;
                        U(1).P.i = 1;
                        U(2).P.name = 'none';
                        U(2).P.h = 0;
                        U(2).P.i = 1;
                    end;
                else
                     if eval(sprintf('nr_para_reg%d>0',ind_con))
                        for l = 1:nr_para
                            U(1).P(l).name = sprintf('par_value_1_%d',l);
                            eval(sprintf('U(1).P(l).P = mod%d_1;',l));
                            eval(evalstr);
                            U(1).P(l).h = 1;
                            U(1).P(l).i = [1,l+1];
                            U(2).P(l).name = sprintf('par_value_2_%d',l);
                            eval(sprintf('U(2).P(l).P = mod%d_2;',l));
                            U(2).P(l).h = 1;
                            U(2).P(l).i = [1,l+1];
                        end;
                    else
                        U(1).P.name = 'none';
                        U(1).P.h = 0;
                        U(1).P.i = 1;
                        U(2).P.name = 'none';
                        U(2).P.h = 0;
                        U(2).P.i = 1;
                    end;                   
                end;
                
             
                U(1).dt = SPM.Sess(sess_count).U(con).dt;
                U(1).u = SPM.Sess(sess_count).U(con).u;
                U(1).pst = SPM.Sess(sess_count).U(con).pst;    

                U(2).name{1} = sprintf('half2_%s',name);
                U(2).ons = onsets2;
                U(2).dur = zeros(length(onsets2),1);
                if isfield(SPM.Sess(sess_count).U(con),'orth')
                    U(2).orth = SPM.Sess(sess_count).U(con).orth;
                end;
                U(2).dt = SPM.Sess(sess_count).U(con).dt;
                U(2).u = SPM.Sess(sess_count).U(con).u;
                U(2).pst = SPM.Sess(sess_count).U(con).pst;

                eval(sprintf('U_%d = U;',ind_con));    
                U = [];
            end;
            
            for ind_con = 1:length(con_vec)
                if ind_con == 1
                    con = con_vec(1,ind_con);
                    Sess(sess_count).U(1:con-1) = SPM.Sess(sess_count).U(1:con-1);
                    eval(sprintf('Sess(sess_count).U(con) = U_%d(1);',ind_con));
                    eval(sprintf('Sess(sess_count).U(con+1) = U_%d(2);',ind_con));
                    if ind_con ~= length(con_vec)
                        next_con = con_vec(1,ind_con+1);
                        Sess(sess_count).U(con+2) = SPM.Sess(sess_count).U(con+1:next_con-1);
                    else
                        U_rest = SPM.Sess(sess_count).U(con+1:end);
                        if length(U_rest) ~= 0
                        Sess(sess_count).U(end+1:(end+length(U_rest))) = U_rest;
                        end;
                    end;                       
                else
                    con = con_vec(1,ind_con);
                    eval(sprintf('Sess(sess_count).U(end+1) = U_%d(1);',ind_con));
                    eval(sprintf('Sess(sess_count).U(end+1) = U_%d(2);',ind_con));
                    if ind_con ~= length(con_vec)
                        next_con = con_vec(1,ind_con+1);
                        rest = SPM.Sess(sess_count).U(con+1:next_con-1);
                        Sess(sess_count).U(end+1:end+length(rest)) = rest;
                    else
                        U_rest = SPM.Sess(sess_count).U(con+1:end);
                        Sess(sess_count).U(end+1:end+length(U_rest)) = U_rest;
                    end;
                end;
            end;
            for sess_count = 1:sess_length
                SPM.Sess(sess_count).U = Sess(sess_count).U;
            end;
        end;

            
        %new stats folder
        cd(stats_path);
        cd(sprintf('%s',id{count}));
        mkdir(split_dir_temp)
        dir_results = [stats_path f id{count} f split_dir_temp];
        SPM.swd = dir_results;
        cd(dir_results);
        save SPM.mat SPM ;

        %% estimation of GLM
        %--->part of spm_fMRI_design
            fprintf('...do statistics for %s in session %d...',id{count},m);

            %-Construct Design matrix {X} %-Microtime onset and microtime resolution
            try
                fMRI_T     = SPM.xBF.T;
                fMRI_T0    = SPM.xBF.T0;
            catch
                fMRI_T     = spm_get_defaults('stats.fmri.t');
                fMRI_T0    = spm_get_defaults('stats.fmri.t0');
                SPM.xBF.T  = fMRI_T;
                SPM.xBF.T0 = fMRI_T0;
            end;
            %-Time units, dt = time bin {secs}
            SPM.xBF.dt     = SPM.xY.RT/SPM.xBF.T;
            %-Get basis functions
            SPM.xBF        = spm_get_bf(SPM.xBF);
            %-Get session specific design parameters
            Xx    = [];
            Xb    = [];
            Xname = {};
            Bname = {};
            for s = 1:length(SPM.nscan)
                %-Number of scans for this session
                k = SPM.nscan(s);
                %-Create convolved stimulus functions or inputs
                %-Get inputs, neuronal causes or stimulus functions U
                U = spm_get_ons(SPM,s);
                %-Convolve stimulus functions with basis functions
                [X,Xn,Fc] = spm_Volterra(U, SPM.xBF.bf, SPM.xBF.Volterra);
                %-Resample regressors at acquisition times (32 bin offset)
                if ~isempty(X)
                    X = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
                end
                %-Orthogonalise (within trial type)
                for i = 1:length(Fc)
                    if i<= numel(U) && ... % for Volterra kernels
                            (~isfield(U(i),'orth') || U(i).orth)
                        p = ones(size(Fc(i).i));
                    else
                        p = Fc(i).p;
                    end
                    for j = 1:max(p)
                        X(:,Fc(i).i(p==j)) = spm_orth(X(:,Fc(i).i(p==j)));
                    end
                end
                %-Get user specified regressors
                C     = SPM.Sess(s).C.C;
                Cname = SPM.Sess(s).C.name;
                %-Append mean-corrected regressors and names
                X     = [X spm_detrend(C)];
                Xn    = {Xn{:} Cname{:}};
                %-Confounds: Session effects
                B     = ones(k,1);
                Bn    = {'constant'};
                %-Session structure array
                SPM.Sess(s).U      = U;
                SPM.Sess(s).C.C    = C;
                SPM.Sess(s).C.name = Cname;
                SPM.Sess(s).row    = size(Xx,1) + (1:k);
                SPM.Sess(s).col    = size(Xx,2) + (1:size(X,2));
                SPM.Sess(s).Fc     = Fc;
                 %-Append into Xx and Xb
                Xx      = blkdiag(Xx,X);
                Xb      = blkdiag(Xb,B);
                %-Append names
                for i = 1:length(Xn)
                    Xname{end + 1} = [sprintf('Sn(%i) ',s) Xn{i}];
                end
                for i = 1:length(Bn)
                    Bname{end + 1} = [sprintf('Sn(%i) ',s) Bn{i}];
                end
            end
            %-Place design matrix structure in xX
            SPM.xX.X    = [Xx Xb];
            SPM.xX.iH   = [];
            SPM.xX.iC   = 1:size(Xx,2);
            SPM.xX.iB   = (1:size(Xb,2)) + size(Xx,2);
            SPM.xX.iG   = [];
            SPM.xX.name = {Xname{:} Bname{:}};

            save SPM.mat SPM;
        
        %--->parts of spm_spm to estimate GLM
             SVNid = '$Rev: 6015 $';
            %-Say hello
            SPMid = spm('FnBanner',mfilename,SVNid);
            spm('Pointer','Watch');
            %-Get SPM
            if ~nargin
                [P, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
                if ~sts, SPM = []; return; end
                swd      = spm_file(P,'fpath');
                load(fullfile(swd,'SPM.mat'));
                SPM.swd  = swd;
            end
            %-Change directory
            try
                cd(SPM.swd);
            catch
                SPM.swd = pwd;
            end
            %-Check input files
            try
                VY = SPM.xY.VY;
            catch
                error('Data have not been specified.');
            end

            for i = 1:numel(VY)
                if ~spm_existfile(VY(i).fname)
                    error('File not found: %s',VY(i).fname);
                end
                if ~spm_mesh_detect(VY)
                    % Backward compatibility: propagate scaling (see spm_fmri_spm_ui.m)
                    VY(i).private.dat.scl_slope = VY(i).pinfo(1);
                    VY(i).private.dat.scl_inter = VY(i).pinfo(2);
                end
            end

            spm_check_orientations(VY);
            M       = VY(1).mat;
            DIM     = VY(1).dim;
            YNaNrep = spm_type(VY(1).dt(1),'nanrep');
            if spm_mesh_detect(VY)
                file_ext = '.gii';
            else
                file_ext = spm_file_ext;
            end
            %-             A N A L Y S I S   P R E L I M I N A R I E S
            %-Get design
            xX             = SPM.xX;
            [nScan, nBeta] = size(xX.X);
            %-Get masking settings
            if isfield(SPM,'xM')
                xM         = SPM.xM;
            else
                xM         = -Inf(nScan,1);
            end
            if ~isstruct(xM)
                xM         = struct(...
                                'T',  [],...
                                'TH', xM,...
                                'I',  0,...
                                'VM', {[]},...
                                'xs', struct('Masking','analysis threshold'));
            end
            mask           = true(DIM);
            %-Check confounds (xX.K)
            if ~isfield(xX,'K')
                xX.K       = 1;
            end
            %-Get non-sphericity (xVi), otherwise assume i.i.d.
            if isfield(SPM,'xVi')
                xVi        = SPM.xVi;
            else
                xVi        = struct('form', 'i.i.d.',...
                                    'V',    speye(nScan,nScan));
            end
            %-Evoke ReML for hyperparameter estimation
            if ~isfield(xVi,'V')
                SPM.xY.VY  = VY;
                SPM.xM     = xM;
                SPM.xX.K   = xX.K;
                [xVi, am]  = spm_est_non_sphericity(SPM);
                mask       = mask & am;
                spm('FnBanner',mfilename,SVNid);
            end
            %-Get weight/whitening matrix:  W*W' = inv(V)
            if isfield(xX,'W')
                W          = xX.W;
            else
                W          = spm_sqrtm(spm_inv(xVi.V));
                W          = W.*(abs(W) > 1e-6);
                xX.W       = sparse(W);
            end
            %-Design space and projector matrix [pseudoinverse] for WLS
            xX.xKXs        = spm_sp('Set',spm_filter(xX.K,W*xX.X));    % KWX
            xX.xKXs.X      = full(xX.xKXs.X);
            xX.pKX         = spm_sp('x-',xX.xKXs);                     % Projector
            erdf           = spm_SpUtil('trRV',xX.xKXs);               % error df
            %-Use non-sphericity xVi.V to compute [effective] degrees of freedom
            xX.V           = spm_filter(xX.K,spm_filter(xX.K,W*xVi.V*W')'); % KWVW'K'
            [trRV, trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);          % trRV (for X)
            xX.trRV        = trRV;                                     % <R'*y'*y*R>
            xX.trRVRV      = trRVRV;                                   %-Satterthwaite
            xX.erdf        = trRV^2/trRVRV;                            % approximation
            xX.Bcov        = xX.pKX*xX.V*xX.pKX';                      % Cov(beta)
            %-             I N I T I A L I S E   O U T P U T   F I L E S
            %-Initialise mask file
            VM = struct(...
                'fname',   ['mask' file_ext],...
                'dim',     DIM,...
                'dt',      [spm_type('uint8') spm_platform('bigend')],...
                'mat',     M,...
                'pinfo',   [1 0 0]',...
                'descrip', 'spm_spm:resultant analysis mask');
            VM = spm_data_hdr_write(VM);
            %-Initialise beta files
            Vbeta(1:nBeta) = deal(struct(...
                'fname',   [],...
                'dim',     DIM,...
                'dt',      [spm_type('float32') spm_platform('bigend')],...
                'mat',     M,...
                'pinfo',   [1 0 0]',...
                'descrip', 'spm_spm:beta'));

            for i = 1:nBeta
                Vbeta(i).fname   = [sprintf('beta_%04d',i) file_ext];
                Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,xX.name{i});
            end
            Vbeta = spm_data_hdr_write(Vbeta);
            %-Initialise residual sum of squares file
            VResMS = struct(...
                'fname',   ['ResMS' file_ext],...
                'dim',     DIM,...
                'dt',      [spm_type('float64') spm_platform('bigend')],...
                'mat',     M,...
                'pinfo',   [1 0 0]',...
                'descrip', 'spm_spm:Residual sum-of-squares');
            VResMS = spm_data_hdr_write(VResMS);
            %-Initialise standardised residual images
            nSres = min(nScan, spm_get_defaults('stats.maxres'));
            resInMem = spm_get_defaults('stats.resmem');
            VResI(1:nSres) = deal(struct(...
                'fname',   [],...
                'dim',     DIM,...
                'dt',      [spm_type('float64') spm_platform('bigend')],...
                'mat',     M,...
                'pinfo',   [1 0 0]',...
                'descrip', 'spm_spm:StandardisedResiduals'));
            if resInMem, for i=1:nSres, VResI(i).dat = zeros(VResI(i).dim); end; end

            for i = 1:nSres
                VResI(i).fname   = [sprintf('ResI_%04d', i) file_ext];
                VResI(i).descrip = sprintf('spm_spm:ResI (%04d)', i);
            end
            VResI = spm_data_hdr_write(VResI);
            %-               G E N E R A L   L I N E A R   M O D E L
            iRes = round(linspace(1,nScan,nSres))';              % Indices for residual
            %-Get explicit mask(s)
            for i = 1:numel(xM.VM)
                %-Assume it fits entirely in memory    
                C = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
                v = true(DIM);
                [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
                for x3 = 1:DIM(3)
                    M2  = inv(M\xM.VM(i).mat);
                    y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
                    y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
                    y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
                    v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
                end
                mask = mask & v;
                clear C v x1 x2 x3 M2 y1 y2 y3
            end
            %-Split data into chunks
            chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
            nbchunks  = ceil(prod(DIM) / chunksize);
            chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

            spm_progress_bar('Init',nbchunks,'Parameter estimation','Chunks');

            for i=1:nbchunks
                chunk = chunks(i):chunks(i+1)-1;
                %-Report progress
                %======================================================================
                if i > 1, fprintf(repmat(sprintf('\b'),1,72)); end                  %-# 
                fprintf('%-40s: %30s', sprintf('Chunk %3d/%-3d',i,nbchunks),...
                                       '...processing');                            %-#
                %-Get data & construct analysis mask
                %======================================================================
                Y     = zeros(nScan,numel(chunk));
                cmask = mask(chunk);
                for j=1:nScan
                    if ~any(cmask), break, end                 %-Break if empty mask
                    Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));%-Read chunk of data
                    cmask(cmask) = Y(j,cmask) > xM.TH(j);      %-Threshold (& NaN) mask
                    if xM.I && ~YNaNrep && xM.TH(j) < 0        %-Use implicit mask
                        cmask(cmask) = abs(Y(j,cmask)) > eps;
                    end
                end
                cmask(cmask) = any(diff(Y(:,cmask),1));        %-Mask constant data
                Y            = Y(:,cmask);                     %-Data within mask

                %-Whiten/Weight data and remove filter confounds
                %======================================================================
                KWY          = spm_filter(xX.K,W*Y);

                %-Weighted Least Squares estimation
                %======================================================================
                beta         = xX.pKX*KWY;                     %-Parameter estimates
                if any(cmask)
                    res      = spm_sp('r',xX.xKXs,KWY);        %-Residuals
                else
                    res      = zeros(nScan,0);
                end
                ResSS        = sum(res.^2);                    %-Residual SSQ
                res          = res(iRes,:);

                %-Write output files
                %======================================================================
                c            = NaN(numel(chunk),1);

                %-Write mask file
                %----------------------------------------------------------------------
                mask(chunk)  = cmask;
                VM           = spm_data_write(VM, cmask', chunk);

                %-Write beta files
                %----------------------------------------------------------------------
                for j=1:nBeta
                    c(cmask) = beta(j,:);
                    Vbeta(j) = spm_data_write(Vbeta(j), c, chunk); 
                end

                %-Write ResSS into ResMS (variance) file scaled by tr(RV)
                %----------------------------------------------------------------------
                c(cmask)     = ResSS / xX.trRV;
                VResMS       = spm_data_write(VResMS, c, chunk);

                %-Write standardised residual files
                %----------------------------------------------------------------------
                for j=1:nSres
                    c(cmask) = res(j,:)./sqrt(ResSS/erdf); % or xX.erdf
                    VResI(j) = spm_data_write(VResI(j), c, chunk);
                end

                %-Report progress
                %======================================================================
                fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');             %-#
                spm_progress_bar('Set',i);
            end

            fprintf('\n');                                                          %-#
            spm_progress_bar('Clear');

            if ~any(mask(:))
                error('Please check your data: There are no inmask voxels.');
            end
            %-                 R e s M S   M O D I F I C A T I O N
            %-Modify ResMS (a form of shrinkage) to avoid problems of very low variance
            try
                if ~strcmpi(spm_get_defaults('modality'),'fmri')
                    ResMS  = spm_data_read(VResMS);
                    ResMS  = ResMS + 1e-3 * max(ResMS(isfinite(ResMS)));
                    VResMS = spm_data_write(VResMS, ResMS);
                    clear ResMS
                end
            end
            %-              S M O O T H N E S S   E S T I M A T I O N
            if ~spm_mesh_detect(VY)
                [FWHM,VRpv,R] = spm_est_smoothness(VResI,VM,[nScan erdf]);
            else
                VRpv = struct(...
                    'fname',   ['RPV' file_ext],...
                    'dim',     DIM,...
                    'dt',      [spm_type('float64') spm_platform('bigend')],...
                    'mat',     M,...
                    'pinfo',   [1 0 0]',...
                    'descrip', 'spm_spm: resels per voxel');
                VRpv = spm_data_hdr_write(VRpv);
                ResI = zeros(prod(DIM),numel(VResI));
                for i=1:numel(VResI)
                    ResI(:,i) = spm_data_read(VResI(i));
                end
                g = gifti(VY(1).fname);
                g = g.private.metadata(1).value;
                if isempty(spm_file(g,'path'))
                    g = fullfile(spm_file(VY(1).fname,'path'),g);
                end
                [R, RPV] = spm_mesh_resels(gifti(g),mask,ResI);
                RPV(~mask) = NaN;
                VRpv = spm_data_write(VRpv,RPV);
                FWHM = [1 1 1] * (1/mean(RPV(mask))).^(1/3);
            end

            %-Delete the standardised residual files
            fres = cellstr(spm_select('FPList',SPM.swd,'^ResI_.{4}\..{3}$'));
            for i=1:numel(fres)
                spm_unlink(fres{i});
            end
            %-                         S A V E   &   E X I T
            %-Compute scaled design matrix for display purposes
                xX.nKX         = spm_DesMtx('sca',xX.xKXs.X,xX.name);

            %-Compute coordinates of voxels within mask
            [x,y,z]        = ind2sub(DIM,find(mask));
            XYZ            = [x y z]';

            %-Place fields in SPM
            SPM.xVol.XYZ   = XYZ;               %-InMask XYZ coords (voxels)
            SPM.xVol.M     = M;                 %-voxels -> mm
            SPM.xVol.iM    = inv(M);            %-mm -> voxels
            SPM.xVol.DIM   = DIM';              %-Image dimensions
            SPM.xVol.FWHM  = FWHM;              %-Smoothness data
            SPM.xVol.R     = R;                 %-Resel counts
            SPM.xVol.S     = nnz(mask);         %-Volume (voxels)
            SPM.xVol.VRpv  = VRpv;              %-Filehandle - Resels per voxel
            if spm_mesh_detect(VY)
                SPM.xVol.G = g;                 %-Mesh topology
            end

            SPM.Vbeta      = Vbeta;             %-Filehandle - Beta
            SPM.VResMS     = VResMS;            %-Filehandle - Hyperparameter
            SPM.VM         = VM;                %-Filehandle - Mask

            SPM.xVi        = xVi;               %-Non-sphericity structure

            SPM.xX         = xX;                %-Design structure

            SPM.xM         = xM;                %-Mask structure

            SPM.xCon       = struct([]);        %-Contrast structure

            SPM.SPMid      = SPMid;
            SPM.swd        = pwd;
            SPM.xX.xKXs.rank = size(SPM.xX.xKXs.X,2);
            %-Save SPM.mat
            %--------------------------------------------------------------------------
            fprintf('%-40s: %30s','Saving SPM.mat','...writing');                   %-#
            save('SPM.mat','SPM', spm_get_defaults('mat.format'));
            fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#

        % building contrasts 
        disp('...creating contrast images...')
        if length(con_vec) == 1
            if  nr_para>0
                cons = cell(1,(length(con_vec).*nr_para)+length(con_vec));
                for o = 1:(length(con_vec)+nr_para)*2
                    eval(sprintf('con%d=zeros(size(SPM.xX.xKXs.X,2),1);',o));
                    eval(sprintf('con%d(%d,1)=1;',o,o));
                    eval(sprintf('cons{o}=con%d;',o));
                end;
            else
                cons = {};
                for ind_cons = 1:length(con_vec)
                    con = con_vec(1,ind_cons);
                    eval(sprintf('con%d_1=zeros(size(SPM.xX.xKXs.X,2),1);',ind_cons));
                    eval(sprintf('con%d_1(%d,1)=1;',ind_cons,con+ind_cons-1));
                    eval(sprintf('cons{end+1}=con%d_1;',ind_cons));
                    eval(sprintf('con%d_2=zeros(size(SPM.xX.xKXs.X,2),1);',ind_cons));
                    eval(sprintf('con%d_2(%d+1,1)=1;',ind_cons,con+ind_cons-1));
                    eval(sprintf('cons{end+1}=con%d_2;',ind_cons));
                end;
            end;  
        else
            cons = {};
            for ind_cons = 1:length(con_vec)
                if  eval(sprintf('nr_para_reg%d>0',ind_cons))
                    nr_para = eval(sprintf('nr_para_reg%d',ind_cons));
                    cons = cell(1,(length(con_vec).*nr_para)+length(con_vec));
                    for o = 1:(length(con_vec)+nr_para)*2
                        eval(sprintf('con%d=zeros(size(SPM.xX.xKXs.X,2),1);',o));
                        eval(sprintf('con%d(%d,1)=1;',o,o));
                        eval(sprintf('cons{o}=con%d;',o));
                    end;
                else
                    %cons = {};
                        con = con_vec(1,ind_cons);
                        eval(sprintf('con%d_1=zeros(size(SPM.xX.xKXs.X,2),1);',ind_cons));
                        eval(sprintf('con%d_1(%d,1)=1;',ind_cons,con+ind_cons-1));
                        eval(sprintf('cons{end+1}=con%d_1;',ind_cons));
                        eval(sprintf('con%d_2=zeros(size(SPM.xX.xKXs.X,2),1);',ind_cons));
                        eval(sprintf('con%d_2(%d+1,1)=1;',ind_cons,con+ind_cons-1));
                        eval(sprintf('cons{end+1}=con%d_2;',ind_cons));
                end;
            end;    
        
        end;
                    
        for ind_con = 1:length(con_vec)
            if length(con_vec) > 1
                nr_para = eval(sprintf('nr_para_reg%d',ind_cons));
            end;
            eval(sprintf('nams_%d = cell(1,2*(1+nr_para));',ind_con));
            nam1 = sprintf('half1_reg%d',con_vec(ind_con));
            eval(sprintf('nams_%d{1}=nam1;',ind_con));
            nam2 = sprintf('half2_reg%d',con_vec(ind_con));
            eval(sprintf('nams_%d{2}=nam2;',ind_con));
            if nr_para > 0
                for p = 2:nr_para+1
                    name1=sprintf('half1_con%d_value%d',con_vec(ind_con),p-1);
                    eval(sprintf('nams_%d{p}=name1;',ind_con));
                end;
                for q = (nr_para+2):(2*nr_para+2)
                    if q == nr_para+2
                        name2 = sprintf('half2_con%d',con_vec(ind_con));
                        eval(sprintf('nams_%d{q}=name2;',ind_con));
                    else
                        name3=sprintf('half2_con%d_value%d',con_vec(ind_con),q-(nr_para+1));
                        eval(sprintf('nams_%d{q} = name3;',con_vec(ind_con)));
                    end;

                end;
            end;
        end;
        nams={};
        for ind_con = 1:length(con_vec)
            eval(sprintf('nams={nams{:},nams_%d{:}};',ind_con));
        end;
                
            DxCon=struct('name',nams,'STAT', 'T', 'c', cons,'X0','','iX0','c','X1o','','eidf','1','Vcon','','Vspm','');    
            sX = SPM.xX.xKXs;
            for r = 1:length(nams)
                Fc = spm_FcUtil('Set',nams{r}, 'T', 'c', cons{r}, sX);
                DxCon(1).X0 =Fc.X0;
                DxCon(1).X1o = Fc.X1o;
                clear Fc
            end;

        %-Append to SPM.xCon
        % SPM will automatically save any contrasts that evaluate successfully
        %------------------------------------------------------------------
        if isempty(SPM.xCon)
            SPM.xCon = DxCon;
        elseif ~isempty(DxCon)
            SPM.xCon(end+1) = DxCon;
        end;
        save SPM.mat SPM
        
        if length(con_vec) == 1
            vec=1:2*(length(con_vec)*nr_para+length(con_vec));
        else
            length_vec = 0;
            for i = 1:length(con_vec)
                length_vec = length_vec + 1 + eval(sprintf('nr_para_reg%d',i));
            end;           
            %vec=1:length(vec);      
            vec=1:length_vec;
        end;
        SPM = spm_contrasts(SPM,vec);

        save SPM.mat SPM;
        clearvars -except box_path f nr_subj new_dir split_dir split_name con_vec study_design oldpointer handles run runs count id n con par path name con1 con2 con3 con4 cwd SPM_list nr_para nr_para_reg1 nr_para_reg2 stats_path stats_dir m o p q r s;
        fprintf('...new statistic for %s in session %d done...\n',id{count},m);
    end;
end;

disp('creating 4D images');
%% Create 4D images out of 3D images 
% uses SPM12
cd(box_path);
% uses template batch 'template_3D-4D.mat'
load ('template_3D-4D.mat');
% define 4D image name - number of run will be added below
name = '4D_split';
% batch name
batch_name = 'batch_3dto4d_split.mat';
% results directory
dir_results = study_design.results_directory;

% create batch
for ind_con = 1:length(con_vec)
    if length(con_vec) > 1
        nr_para = eval(sprintf('nr_para_reg%d',ind_con));
    end;
    if ind_con == 1
        cont = 1;
    else
        cont = (2*((ind_con-1)*nr_para + ind_con-1)) + 1;
    end;
    
    for i = 1:runs
        newstats_dir = sprintf(stats_dir,i); 
        new_dir = sprintf('%s_split_%s',newstats_dir,split_name);
        for j=1:length(id)
            con_img1 = sprintf('%s%s%s%s%s%s%s%s%scon_%04d.nii,1',stats_path,f,f,id{j},f,f,new_dir,f,f,cont);
            matlabbatch{i}.spm.util.cat.vols{j,1} = con_img1;
            con_img2 = sprintf('%s%s%s%s%s%s%s%s%scon_%04d.nii,1',stats_path,f,f,id{j},f,f,new_dir,f,f,cont+1+nr_para);
            matlabbatch{i+1}.spm.util.cat.vols{j,1} = con_img2;
        end;
    
        if length(con_vec) == 1
            img_name1 = sprintf('%s1_%d.nii',name,i);
            matlabbatch{i}.spm.util.cat.name = img_name1;
            matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
            img_name2 = sprintf('%s2_%d.nii',name,i);
            matlabbatch{i+1}.spm.util.cat.name = img_name2;
            matlabbatch{i+1}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
        else
            img_name1 = sprintf('%s1_reg%d_%d.nii',name,con_vec(ind_con),i);
            matlabbatch{i}.spm.util.cat.name = img_name1;
            matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
            img_name2 = sprintf('%s2_reg%d_%d.nii',name,con_vec(ind_con),i);
            matlabbatch{i+1}.spm.util.cat.name = img_name2;
            matlabbatch{i+1}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images

        end;
% save batch
save(batch_name,'matlabbatch');

% run batch
spm_jobman('serial',batch_name);
    end;

end;



% move files
% SPM saves 4D files by default in stats directory of first subject
% move to directory defined for results
for ind_con = 1:length(con_vec)
    if length(con_vec) > 1
        nr_para = eval(sprintf('nr_para_reg%d',ind_con));
    end;    
    for k = 1:runs
        newstats_dir = sprintf(stats_dir,k); 
        new_dir = sprintf('%s_split_%s',newstats_dir,split_name);

        if length(con_vec) > 1
            img_name1 = sprintf('%s1_reg%d_%d',name,con_vec(ind_con),k);
            img_name2 = sprintf('%s2_reg%d_%d',name,con_vec(ind_con),k); 
            file1 = [stats_path f id{1} f new_dir f img_name1 '.nii'];
            file2 = [stats_path f id{1} f new_dir f img_name1 '.mat'];
            file3 = [stats_path f id{1} f new_dir f img_name2 '.nii'];
            file4 = [stats_path f id{1} f new_dir f img_name2 '.mat'];
            movefile (file1,dir_results,'f');
            movefile (file2,dir_results,'f');
            movefile (file3,dir_results,'f');
            movefile (file4,dir_results,'f');
        else
            img_name1 = sprintf('%s1_%d',name,k);
            img_name2 = sprintf('%s2_%d',name,k); 
            file1 = [stats_path f id{1} f new_dir f img_name1 '.nii'];
            file2 = [stats_path f id{1} f new_dir f img_name1 '.mat'];
            file3 = [stats_path f id{1} f new_dir f img_name2 '.nii'];
            file4 = [stats_path f id{1} f new_dir f img_name2 '.mat'];
            movefile (file1,dir_results,'f');
            movefile (file2,dir_results,'f');
            movefile (file3,dir_results,'f');
            movefile (file4,dir_results,'f');
        end;
        
    end;

%create 4D and move it for parametric modulators
if nr_para>0
        for ind_par = 1:nr_para
            if ind_con == 1
                cont = 1+ind_par;
            else
                cont = (2*((ind_con-1)*nr_para + ind_con-1)) + 1 + ind_par;
            end;
            % create batch
            for i = 1:2:(runs*2)-1
                sess = round(i./2);
                newstats_dir = sprintf(stats_dir,sess); 
                new_dir = sprintf('%s_split_%s',newstats_dir,split_name);
                for j=1:length(id)
                    con_img1 = sprintf('%s%s%s%s%s%s%s%s%scon_%04d.nii,1',stats_path,f,f,id{j},f,f,new_dir,f,f,cont);
                    matlabbatch{i}.spm.util.cat.vols{j,1} = con_img1;
                    con_img2 = sprintf('%s%s%s%s%s%s%s%s%scon_%04d.nii,1',stats_path,f,f,id{j},f,f,new_dir,f,f,cont+1+nr_para);
                    matlabbatch{i+1}.spm.util.cat.vols{j,1} = con_img2;
                end;
                    img_name1 = sprintf('%s1_par%d_%d.nii',name,ind_par,sess);
                    matlabbatch{i}.spm.util.cat.name = img_name1;
                    matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
                    img_name2 = sprintf('%s2_par%d_%d.nii',name,ind_par,sess);
                    matlabbatch{i+1}.spm.util.cat.name = img_name2;
                    matlabbatch{i+1}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
            end;

            % save batch
            save(batch_name,'matlabbatch');

            % run batch
            spm_jobman('serial',batch_name);

            % move files
            % SPM saves 4D files by default in stats directory of first subject
            % move to directory defined for results
            for k = 1:runs
                    newstats_dir = sprintf(stats_dir,k); 
                    new_dir = sprintf('%s_split_%s',newstats_dir,split_name);
                    img_name1 = sprintf('%s1_par%d_%d',name,ind_par,k);
                    img_name2 = sprintf('%s2_par%d_%d',name,ind_par,k); 
                    file1 = [stats_path f id{1} f new_dir f img_name1 '.nii'];
                    file2 = [stats_path f id{1} f new_dir f img_name1 '.mat'];
                    file3 = [stats_path f id{1} f new_dir f img_name2 '.nii'];
                    file4 = [stats_path f id{1} f new_dir f img_name2 '.mat'];
                    movefile (file1,dir_results,'f');
                    movefile (file2,dir_results,'f');
                    movefile (file3,dir_results,'f');
                    movefile (file4,dir_results,'f');
            end;
        end;
end;
end;


cd(box_path);
disp('...DONE');


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
axes(hObject);
imshow('logo.png');
