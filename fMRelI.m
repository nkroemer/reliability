function varargout = fMRelI(varargin)
% FMRELI MATLAB code for fMRelI.fig
%      FMRELI, by itself, creates a new FMRELI or raises the existing
%      singleton*.
%
%      H = FMRELI returns the handle to a new FMRELI or the handle to
%      the existing singleton*.
%
%      FMRELI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FMRELI.M with the given input arguments.
%
%      FMRELI('Property','Value',...) creates a new FMRELI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fMRelI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fMRelI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fMRelI

% Last Modified by GUIDE v2.5 26-Jul-2017 09:14:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fMRelI_OpeningFcn, ...
                   'gui_OutputFcn',  @fMRelI_OutputFcn, ...
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


% --- Executes just before fMRelI is made visible.
function fMRelI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fMRelI (see VARARGIN)

% Choose default command line output for fMRelI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fMRelI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = fMRelI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
disp('Starting fMRelI...');
set(handles.specify,'TooltipString','Specification of sample, data folders etc.');
set(handles.define,'TooltipString','Specification of contrast names and regressor numbers');
set(handles.split,'TooltipString','Use split-half procedure to look at within-session reliability');
set(handles.CorrMaps,'TooltipString','Calculation of correlations and ICCs voxel-/ROI-wise');
set(handles.dice_ROI,'TooltipString','Calculation of Dice-/Jaccard-coefficients');
set(handles.similarity,'TooltipString','Calculation of global similarity measures');
set(handles.summary,'TooltipString','Creates reliability summary for all ROIs in atlas');
set(handles.threshold,'TooltipString','Creates ROIs with certain reliability threshold');






% --- Executes on button press in specify.
function specify_Callback(hObject, eventdata, handles)
% hObject    handle to specify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
design;

% --- Executes on button press in CorrMaps.
function CorrMaps_Callback(hObject, eventdata, handles)
% hObject    handle to CorrMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CorrICC;


% --- Executes on button press in split.
function split_Callback(hObject, eventdata, handles)
% hObject    handle to split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SplitHalf;


% --- Executes on button press in dice_ROI.
function dice_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to dice_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Overlap;

% --- Executes on button press in define.
function define_Callback(hObject, eventdata, handles)
% hObject    handle to define (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Contrast_Def;

% --- Executes on button press in similarity.
function similarity_Callback(hObject, eventdata, handles)
% hObject    handle to similarity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Similarity;

% --- Executes on button press in summary.
function summary_Callback(hObject, eventdata, handles)
% hObject    handle to summary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Summary;

% --- Executes on button press in threshold.
function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Corr2ROI;

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
axes(hObject);
imshow('logo.png');
