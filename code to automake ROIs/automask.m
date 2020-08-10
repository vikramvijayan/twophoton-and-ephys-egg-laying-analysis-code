function varargout = automask(varargin)
% AUTOMASK MATLAB code for automask.fig
%      AUTOMASK, by itself, creates a new AUTOMASK or raises the existing
%      singleton*.
%
%      H = AUTOMASK returns the handle to a new AUTOMASK or the handle to
%      the existing singleton*.
%
%      AUTOMASK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOMASK.M with the given input arguments.
%
%      AUTOMASK('Property','Value',...) creates a new AUTOMASK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before automask_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to automask_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help automask

% Last Modified by GUIDE v2.5 08-May-2019 13:36:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @automask_OpeningFcn, ...
                   'gui_OutputFcn',  @automask_OutputFcn, ...
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
end

% --- Executes just before automask is made visible.
function automask_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to automask (see VARARGIN)

% Choose default command line output for automask
handles.output = hObject;
handles.Ym = varargin{1};
handles.filename = varargin{2};
handles.size_mat = varargin{3};
handles.minsize = 1;
handles.threshold = 375;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes automask wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end

% --- Outputs from this function are returned to the command line.
function varargout = automask_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.threshold = str2double(get(hObject,'String'));

Ym_bw = handles.Ym;
Ym = handles.Ym; 
    axes(handles.axes1); hold on;

imagesc(flipud(Ym),[0,512]); colormap('jet'); axis tight; colorbar;

Ym_bw(Ym<handles.threshold) = 0;
Ym_bw(Ym>=handles.threshold) = 1;
    axes(handles.axes2); hold on;

    imshow(Ym_bw); axis tight;

BW2 = bwareaopen(Ym_bw, handles.minsize);
cc = bwconncomp(BW2);
L = labelmatrix(cc);
RGB = label2rgb(L);
axes(handles.axes3); hold on;
imshow(RGB); axis tight;
ROI_masks_from_GUI = [];

for i =1:1:length(cc.PixelIdxList)
    tmp_v = zeros(cc.ImageSize(1)*cc.ImageSize(2),1);
    tmp_v(cc.PixelIdxList{i}) = 1;
    for j = 1:1:handles.size_mat(5)
            ROI_masks_from_GUI(:,:,j,i) = reshape(tmp_v, cc.ImageSize);
    end
end
handles.ROI_masks = ROI_masks_from_GUI;

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.minsize = str2double(get(hObject,'String'));

Ym_bw = handles.Ym;
Ym = handles.Ym; 
    axes(handles.axes1); hold on;

imagesc(flipud(Ym),[0,512]); colormap('jet'); axis tight; colorbar;

Ym_bw(Ym<handles.threshold) = 0;
Ym_bw(Ym>=handles.threshold) = 1;
    axes(handles.axes2); hold on;

    imshow(Ym_bw); axis tight;

BW2 = bwareaopen(Ym_bw, handles.minsize);
cc = bwconncomp(BW2);
L = labelmatrix(cc);
RGB = label2rgb(L);
axes(handles.axes3); hold on;
imshow(RGB); axis tight;
ROI_masks_from_GUI = [];

for i =1:1:length(cc.PixelIdxList)
    tmp_v = zeros(cc.ImageSize(1)*cc.ImageSize(2),1);
    tmp_v(cc.PixelIdxList{i}) = 1;
    for j = 1:1:handles.size_mat(5)
            ROI_masks_from_GUI(:,:,j,i) = reshape(tmp_v, cc.ImageSize);
    end
end
handles.ROI_masks = ROI_masks_from_GUI;
guidata(hObject, handles);

end
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



end


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.ROI_masks;
ROI_masks_from_GUI = handles.ROI_masks;

save([handles.filename '.mat'], 'ROI_masks_from_GUI');

end