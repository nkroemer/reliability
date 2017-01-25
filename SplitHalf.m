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

% Last Modified by GUIDE v2.5 17-Jan-2017 13:34:34

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

function nr_subj_Callback(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr_subj as text
%        str2double(get(hObject,'String')) returns contents of nr_subj as a double

% --- Executes during object creation, after setting all properties.
function nr_subj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in list_spm.
function list_spm_Callback(hObject, eventdata, handles)
% hObject    handle to list_spm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n =  str2double(get(handles.nr_subj,'String'));
list_spm = cellstr(spm_select(n,'mat','add list with SPM paths'));
assignin('base','list_spm',list_spm);

function subj_list_Callback(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr_subj as text
%        str2double(get(hObject,'String')) returns contents of nr_subj as a double

% --- Executes during object creation, after setting all properties.
function subj_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nr_con_Callback(hObject, eventdata, handles)
% hObject    handle to nr_con (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr_con as text
%        str2double(get(hObject,'String')) returns contents of nr_con as a double


% --- Executes during object creation, after setting all properties.
function nr_con_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr_con (see GCBO)
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



function path_Callback(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path as text
%        str2double(get(hObject,'String')) returns contents of path as a double


% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stats_Callback(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stats as text
%        str2double(get(hObject,'String')) returns contents of stats as a double


% --- Executes during object creation, after setting all properties.
function stats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con1_Callback(hObject, eventdata, handles)
% hObject    handle to con1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con1 as text
%        str2double(get(hObject,'String')) returns contents of con1 as a double


% --- Executes during object creation, after setting all properties.
function con1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con2_Callback(hObject, eventdata, handles)
% hObject    handle to con2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con2 as text
%        str2double(get(hObject,'String')) returns contents of con2 as a double


% --- Executes during object creation, after setting all properties.
function con2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con3_Callback(hObject, eventdata, handles)
% hObject    handle to con3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con3 as text
%        str2double(get(hObject,'String')) returns contents of con3 as a double


% --- Executes during object creation, after setting all properties.
function con3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con4_Callback(hObject, eventdata, handles)
% hObject    handle to con4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con4 as text
%        str2double(get(hObject,'String')) returns contents of con4 as a double


% --- Executes during object creation, after setting all properties.
function con4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con4 (see GCBO)
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
cwd=pwd;
subjects = get(handles.subj_list,'String');
load(subjects);
n = str2double(get(handles.nr_subj,'String'));
con = str2double(get(handles.nr_con,'String'));
par = get(handles.par,'value');
path = get(handles.path,'String');
name = get(handles.stats,'String');

for count = 1:n
%%---loads and modifies SPM for each participant
    list_spm=evalin('base','list_spm');
    dir_spm = list_spm{count};
    load(dir_spm);
    
    %split onsets
    onsets = SPM.Sess.U(con).ons;
    pos = randperm(length(onsets))';
    onsets1 = onsets(pos(1:end/2));
    onsets2 = onsets(pos(end/2:end));
    
    %split parametric value
    if par==1
    sv = SPM.Sess.U(con).P.P;
    sv1 = sv(pos(1:end/2));
    sv2 = sv(pos(end/2:end));
    end;
    
    %create new conditions 
    U(1).name{1} = 'half1';
    U(1).ons = onsets1;
    U(1).dur = zeros(length(onsets1),1);
    U(1).orth = 1;
    if par==1
        U(1).P.name = 'value';
        U(1).P.P = sv1;
        U(1).P.h = 1;
        U(1).P.i = [1,2];
    else
        U(1).P.name = 'none';
        U(1).P.h = 0;
        U(1).P.i = 1;
    end;
    U(1).dt = SPM.Sess.U(con).dt;
    U(1).u = SPM.Sess.U(con).u;
    U(1).pst = SPM.Sess.U(con).pst;

    U(2).name{1} = 'half2';
    U(2).ons = onsets2;
    U(2).dur = zeros(length(onsets2),1);
    U(2).orth = 1;
    if par==1
        U(2).P.name = 'value';
        U(2).P.P = sv2;
        U(2).P.h = 1;
        U(2).P.i = [1,2];
    else
        U(2).P.name = 'none';
        U(2).P.h = 0;
        U(2).P.i = 1;
    end;
    U(2).dt = SPM.Sess.U(con).dt;
    U(2).u = SPM.Sess.U(con).u;
    U(2).pst = SPM.Sess.U(con).pst;
    
    SPM.Sess.U = U;

    
    dir_subj = sprintf('%s\\%d',path,vp(count));
    cd(dir_subj);
    
    %new stats folder
    dir_stats = sprintf('%s\\%s',dir_subj,name);
    mkdir(dir_stats);
    SPM.swd = dir_stats;
    evalstr = sprintf('save %s\\SPM.mat SPM',dir_stats);
    eval(evalstr);
    
    %estimation of GLM
    cd(dir_stats);
    %% part of spm_fMRI_design
        %-Construct Design matrix {X}
    %-Microtime onset and microtime resolution
    try
        fMRI_T     = SPM.xBF.T;
        fMRI_T0    = SPM.xBF.T0;
    catch
        fMRI_T     = spm_get_defaults('stats.fmri.t');
        fMRI_T0    = spm_get_defaults('stats.fmri.t0');
        SPM.xBF.T  = fMRI_T;
        SPM.xBF.T0 = fMRI_T0;
    end
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
    %% parts of spm_spm to estimate GLM
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
    %-            C H E C K   F I L E S   A N D   F O L D E R S
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
    if  par==1
        con1=zeros(size(SPM.xX.xKXs.X,2),1);
        con1(1,1)=1;
        con2=zeros(size(SPM.xX.xKXs.X,2),1);
        con2(2,1)=1;
        con3=zeros(size(SPM.xX.xKXs.X,2),1);
        con3(3,1)=1;
        con4=zeros(size(SPM.xX.xKXs.X,2),1);
        con4(4,1)=1;
    else
        con1=zeros(size(SPM.xX.xKXs.X,2),1);
        con1(1,1)=1;
        con2=zeros(size(SPM.xX.xKXs.X,2),1);
        con2(2,1)=1;
    end;    
    
        if par == 1
            nam1 = 'half1';
            nam2 = 'half1_value';
            nam3 = 'half2';
            nam4 = 'half2_value';
            DxCon=struct('name',{nam1,nam2,nam3,nam4},'STAT', 'T', 'c', {con1,con2,con3,con4},'X0','','iX0','c','X1o','','eidf','1','Vcon','','Vspm','');    
            sX = SPM.xX.xKXs;
            for m = 1:4
                if m == 1
                    Fc = spm_FcUtil('Set',nam1, 'T', 'c', con1, sX);
                    DxCon(1).X0 =Fc.X0;
                    DxCon(1).X1o = Fc.X1o;
                    clear Fc
                elseif m == 2 
                    Fc = spm_FcUtil('Set',nam2, 'T', 'c', con2, sX);
                    DxCon(2).X0 =Fc.X0;
                    DxCon(2).X1o = Fc.X1o;
                    clear Fc
                elseif m == 3
                    Fc = spm_FcUtil('Set',nam3, 'T', 'c', con3, sX);
                    DxCon(3).X0 =Fc.X0;
                    DxCon(3).X1o = Fc.X1o;
                    clear Fc
                elseif m == 4
                    Fc = spm_FcUtil('Set',nam4, 'T', 'c', con4, sX);
                    DxCon(4).X0 =Fc.X0;
                    DxCon(4).X1o = Fc.X1o;
                    clear Fc
                end;
            end;
        else
            nam1 = 'half1';
            nam2 = 'half2';
            DxCon=struct('name',{nam1,nam2},'STAT', 'T', 'c', {con1,con2},'X0','','iX0','c','X1o','','eidf','1','Vcon','','Vspm','');    
            sX = SPM.xX.xKXs;
            for m = 1:2
                if m == 1
                    Fc = spm_FcUtil('Set',nam1, 'T', 'c', con1, sX);
                    DxCon(1).X0 =Fc.X0;
                    DxCon(1).X1o = Fc.X1o;
                    clear Fc
                elseif m == 2 
                    Fc = spm_FcUtil('Set',nam2, 'T', 'c', con2, sX);
                    DxCon(2).X0 =Fc.X0;
                    DxCon(2).X1o = Fc.X1o;
                    clear Fc
                end;
            end;
        end;
                        
        %-Append to SPM.xCon
        % SPM will automatically save any contrasts that evaluate successfully
        %------------------------------------------------------------------
        if isempty(SPM.xCon)
            SPM.xCon = DxCon;
        elseif ~isempty(DxCon)
            SPM.xCon(end+1) = DxCon;
        end
        SPM = spm_contrasts(SPM,[1,2,3,4]);
    
        i1=SPM.Sess.Fc(1).i;
        i2=SPM.Sess.Fc(2).i;

        p1=SPM.Sess.Fc(1).p;
        p2=SPM.Sess.Fc(2).p;
        
        if par==1
        Fc=struct('i',{i1,i2},'name',{nam1,nam3},'p',{p1,p2});
        else
        Fc=struct('i',{i1,i2},'name',{nam1,nam2},'p',{p1,p2});
        end;
        
        SPM.Sess.Fc = Fc;

        save SPM.mat SPM;
 
    clearvars -except run count vp n con par path name con1 con2 con3 con4 cwd;
    fprintf('...new statistic for %d done...',vp(count));
end;
disp('... done');
cd(cwd);