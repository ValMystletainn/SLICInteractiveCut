function varargout = ui(varargin)
% UI MATLAB code for ui.fig
%      UI, by itself, creates a new UI or raises the existing
%      singleton*.
%
%      H = UI returns the handle to a new UI or the handle to
%      the existing singleton*.
%
%      UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UI.M with the given input arguments.
%
%      UI('Property','Value',...) creates a new UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ui

% Last Modified by GUIDE v2.5 12-Dec-2019 02:19:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ui_OpeningFcn, ...
                   'gui_OutputFcn',  @ui_OutputFcn, ...
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


% --- Executes just before ui is made visible.
function ui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ui (see VARARGIN)

% Choose default command line output for ui
handles.output = hObject;

% Update handles structure
handles.isDrawing = false;
handles.forePoly = []; handles.backPoly = [];
handles.ImgOriginal = [];
handles.ISLIC = [];
handles.maskedI = [];
handles.mannulSeedsPoly = [];
guidata(hObject, handles);

% UIWAIT makes ui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in openBtn.
function openBtn_Callback(hObject, eventdata, handles)
% hObject    handle to openBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, path] = uigetfile('*.*');
handles.OriginalFullFile = fullfile(path, filename);
handles.ImgOriginal = im2double(imread(handles.OriginalFullFile));
guidata(hObject, handles);
axes(handles.axeso);
imshow(handles.ImgOriginal)
clear handles.L handles.num
for i = 1 : length(handles.forePoly)
    delete(handles.forePoly(i));
end
for i = 1 : length(handles.backPoly)
    delete(handles.backPoly(i));
end
for i = 1 : length(handles.mannulSeedsPoly)
    delete(handles.mannulSeedsPoly(i));
end
handles.forePoly = []; handles.backPoly = []; handles.mannulSeedsPoly = [];
guidata(hObject, handles);

% --- Executes on button press in saveBtn.
function saveBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, path] = uiputfile('*.jpg');
fullPath = fullfile(path, filename);
imwrite(handles.maskedI, fullPath);
guidata(hObject, handles);

 % --- support function ---
 % use to combine color features like rgb, lab, ... in a m by n by x matrix
function Ifeature = CombineFeatures(handles)
Ifeature = [];
labI = rgb2lab(handles.ImgOriginal);

if handles.LCbx.Value == 1
    Ifeature = cat(3, Ifeature, labI(:,:,1));
end

if handles.ACbx.Value == 1
    Ifeature = cat(3, Ifeature, labI(:,:,2));
end

if handles.labBCbx.Value == 1
    Ifeature = cat(3, Ifeature, labI(:,:,3));
end
sumI = sum(handles.ImgOriginal, 3);
if handles.RCbx.Value == 1
    Ifeature = cat(3, Ifeature, handles.ImgOriginal(:,:,1) ./ sumI * 100);
end
if handles.GCbx.Value == 1
    Ifeature = cat(3, Ifeature, handles.ImgOriginal(:,:,2) ./ sumI * 100);
end
if handles.BCbx.Value == 1
    Ifeature = cat(3, Ifeature, handles.ImgOriginal(:,:,3) ./ sumI * 100);
end
hsv = rgb2hsv(handles.ImgOriginal);
if handles.SCbx.Value == 1
    Ifeature = cat(3, Ifeature, hsv * 100);
end
    

% handles.Ifeature = Ifeature;
% guidata(hObject, handles);

% --- Executes on button press in SLICBtn.
function SLICBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SLICBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = CombineFeatures(handles);
k = str2double(handles.kEdit.String);
m = str2double(handles.mEdit.String);
info = ['k:', num2str(k), '  m: ', num2str(m), '  feature layers: ', num2str(size(I, 3))];
disp(info)
% free the memory of the deleted poly
deletedMask = false(1, length(handles.mannulSeedsPoly));
for i = 1 : length(handles.mannulSeedsPoly)
    if (~isvalid(handles.mannulSeedsPoly(i)))
        deletedMask(i) = true;
    end
end
handles.mannulSeedsPoly(deletedMask) = [];

MannulSeeds = [];
disp('SLIC..')
if (isempty(handles.mannulSeedsPoly))
    [L, num, centerFeatures] = mySLIC(I, k, m);
else    
    for i = 1 : length(handles.mannulSeedsPoly)
        polyobj = handles.mannulSeedsPoly(i);
        for j = 1 : size(polyobj.Position, 1)
           row = round(polyobj.Position(j, 2));
           col = round(polyobj.Position(j, 1));
           MannulSeeds = [MannulSeeds;row col];
        end
    end
    [L, num, centerFeatures] = mySLICMannul(I, k, m, MannulSeeds);
end
isModify = str2double(handles.modifyEdit.String);
if (isModify ~= 0)
    disp('post process..')
    tic
    L = postProcess(L); 
    toc
end
handles.L = L;
handles.num = num;
handles.centerFeatures = centerFeatures;
guidata(hObject, handles);
edgeMask = boundarymask(L);
edgeInd = find(edgeMask);
newI = handles.ImgOriginal;
[M, N,~] = size(I);
newI(edgeInd) = 1;
newI(edgeInd + M * N) = 1;
newI(edgeInd + 2 * M * N) = 1;
axes(handles.axess);
imshow(newI)
disp('done');
handles.ISLIC = newI;
guidata(hObject, handles);


handles.mannulSeedsPoly = [];
for i = 1 : length(MannulSeeds)
    h = drawpolyline(gca, 'Color', 'red', 'Position', [MannulSeeds(i, 2), MannulSeeds(i, 1)]);
    handles.mannulSeedsPoly = [handles.mannulSeedsPoly, h];
end
guidata(hObject, handles);

% --- Executes on button press in foreBtn.
function foreBtn_Callback(hObject, eventdata, handles)
% hObject    handle to foreBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.isDrawing)
    return
end
axes(handles.axeso);
handles.isDrawing = true;
guidata(hObject, handles);
h = drawpolyline(gca, 'Color', 'blue');
handles.isDrawing = false;
guidata(hObject, handles);
if ~isfield(handles,'forePoly')
    handles.forePoly = [];
end
handles.forePoly = [handles.forePoly, h];
guidata(hObject, handles);

% --- Executes on button press in backBtn.
function backBtn_Callback(hObject, eventdata, handles)
% hObject    handle to backBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.isDrawing)
    return
end
axes(handles.axeso);
handles.isDrawing = true;
guidata(hObject, handles);
h = drawpolyline(gca, 'Color', 'black');
handles.isDrawing = false;
guidata(hObject, handles);
if ~isfield(handles,'backPoly')
    handles.backPoly = [];
end
handles.backPoly = [handles.backPoly, h];
guidata(hObject, handles);

% --- Executes on button press in lazySnapBtn.
function lazySnapBtn_Callback(hObject, eventdata, handles)
% hObject    handle to lazySnapBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('lazy sanpping...')
tic
foreMask = false(size(handles.ImgOriginal, 1), size(handles.ImgOriginal, 2));
if isfield(handles,'forePoly')
   for i = 1 : length(handles.forePoly)
       polyobj = handles.forePoly(i);
       if (~isvalid(polyobj))
            continue
       end
       % foreMask = [foreMask; zeros(size(polyobj.Position, 1), 1)];
       for j = 1 : size(polyobj.Position, 1)
           m = round(polyobj.Position(j, 2));
           n = round(polyobj.Position(j, 1));
           foreMask(m, n) = true;
       end
   end
end

backMask = false(size(handles.ImgOriginal, 1), size(handles.ImgOriginal, 2));
if isfield(handles,'backPoly')
   for i = 1 : length(handles.backPoly)
       polyobj = handles.backPoly(i);
       if (~isvalid(polyobj))
            continue
       end
       for j = 1 : size(polyobj.Position, 1)
           m = round(polyobj.Position(j, 2));
           n = round(polyobj.Position(j, 1));
           backMask(m, n) = true;
       end
   end
end

mask = lazysnapping(handles.ImgOriginal, handles.L, foreMask, backMask);
maskedI = handles.ImgOriginal;
maskedI(repmat(~mask,[1 1 3])) = 0;
axes(handles.axesm);
imshow(maskedI)
toc
handles.maskedI = maskedI;
guidata(hObject, handles);


% --- Executes on button press in deleteBtn.
function deleteBtn_Callback(hObject, eventdata, handles)
% hObject    handle to deleteBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear handles.L handles.num
for i = 1 : length(handles.forePoly)
    delete(handles.forePoly(i));
end
for i = 1 : length(handles.backPoly)
    delete(handles.backPoly(i));
end
for i = 1 : length(handles.mannulSeedsPoly)
    delete(handles.mannulSeedsPoly(i));
end
handles.forePoly = []; handles.backPoly = []; handles.mannulSeedsPoly = [];
guidata(hObject, handles);

% --- Executes on button press in showProcessBtn.
function showProcessBtn_Callback(hObject, eventdata, handles)
% hObject    handle to showProcessBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = CombineFeatures(handles);
k = str2double(handles.kEdit.String);
m = str2double(handles.mEdit.String);
tempf = figure;
dynamicSLIC(I, k, m, handles.ImgOriginal);
% close(tempf)

% --- Executes on button press in addSeedBtn.
function addSeedBtn_Callback(hObject, eventdata, handles)
% hObject    handle to addSeedBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.isDrawing)
    return
end
axes(handles.axess);
handles.isDrawing = true;
guidata(hObject, handles);
h = drawpolyline(gca, 'Color', 'red');
handles.isDrawing = false;
guidata(hObject, handles);
if ~isfield(handles,'mannulSeedsPoly')
    handles.mannulSeedsPoly = [];
end
handles.mannulSeedsPoly = [handles.mannulSeedsPoly, h];
guidata(hObject, handles);

%% useless callback



function mEdit_Callback(hObject, eventdata, handles)
% hObject    handle to mEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mEdit as text
%        str2double(get(hObject,'String')) returns contents of mEdit as a double


% --- Executes during object creation, after setting all properties.
function mEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kEdit_Callback(hObject, eventdata, handles)
% hObject    handle to kEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kEdit as text
%        str2double(get(hObject,'String')) returns contents of kEdit as a double


% --- Executes during object creation, after setting all properties.
function kEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function modifyEdit_Callback(hObject, eventdata, handles)
% hObject    handle to modifyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modifyEdit as text
%        str2double(get(hObject,'String')) returns contents of modifyEdit as a double


% --- Executes during object creation, after setting all properties.
function modifyEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modifyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rgbCbx.
function rgbCbx_Callback(hObject, eventdata, handles)
% hObject    handle to rgbCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rgbCbx


% --- Executes on button press in labCbx.
function labCbx_Callback(hObject, eventdata, handles)
% hObject    handle to labCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of labCbx


% --- Executes on button press in RCbx.
function RCbx_Callback(hObject, eventdata, handles)
% hObject    handle to RCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RCbx


% --- Executes on button press in GCbx.
function GCbx_Callback(hObject, eventdata, handles)
% hObject    handle to GCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GCbx


% --- Executes on button press in BCbx.
function BCbx_Callback(hObject, eventdata, handles)
% hObject    handle to BCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BCbx


% --- Executes on button press in LCbx.
function LCbx_Callback(hObject, eventdata, handles)
% hObject    handle to LCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LCbx


% --- Executes on button press in ACbx.
function ACbx_Callback(hObject, eventdata, handles)
% hObject    handle to ACbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ACbx


% --- Executes on button press in labBCbx.
function labBCbx_Callback(hObject, eventdata, handles)
% hObject    handle to labBCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of labBCbx


% --- Executes on button press in SCbx.
function SCbx_Callback(hObject, eventdata, handles)
% hObject    handle to SCbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SCbx
