function varargout = LoadTrimapGUI(varargin)
% LOADTRIMAPGUI MATLAB code for LoadTrimapGUI.fig
%      LOADTRIMAPGUI, by itself, creates a new LOADTRIMAPGUI or raises the existing
%      singleton*.
%
%      H = LOADTRIMAPGUI returns the handle to a new LOADTRIMAPGUI or the handle to
%      the existing singleton*.
%
%      LOADTRIMAPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADTRIMAPGUI.M with the given input arguments.
%
%      LOADTRIMAPGUI('Property','Value',...) creates a new LOADTRIMAPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LoadTrimapGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LoadTrimapGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LoadTrimapGUI

% Last Modified by GUIDE v2.5 25-Feb-2011 16:26:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoadTrimapGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LoadTrimapGUI_OutputFcn, ...
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


% --- Executes just before LoadTrimapGUI is made visible.
function LoadTrimapGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoadTrimapGUI (see VARARGIN)

% Choose default command line output for LoadTrimapGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LoadTrimapGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LoadTrimapGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadSource.
function LoadSource_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[SourceFileName,SourcePathName] = uigetfile('*.jpg;*.bmp;*.png','Select the source image');
filePath = strcat(SourcePathName, SourceFileName);
handles.SrcImage = imread(filePath);
guidata(hObject, handles);


% --- Executes on button press in LoadTrimap.
function LoadTrimap_Callback(hObject, eventdata, handles)
% hObject    handle to LoadTrimap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[TrimapFileName,TrimapPathName] = uigetfile('*.jpg;*.bmp;*.png','Select the source image');
filePath = strcat(TrimapPathName, TrimapFileName);
handles.TriImage = imread(filePath);
imshow(handles.TriImage,'Parent',handles.axes2);
guidata(hObject, handles);

I = handles.SrcImage;
IT = handles.TriImage;

% Convert color image to gray image.
% See http://en.wikipedia.org/wiki/Grayscale for "HowTo".
GrayI = 0.3 * I(:,:,1) + 0.59 * I(:,:,2) + 0.11 * I(:,:,3);
imshow(GrayI,'Parent',handles.axes1);

iWidth = size(GrayI,1);
iHeight = size(GrayI,2);
% Both F-B and I are measured in grayscale channel.
RawFore = zeros(iWidth, iHeight);
Fore = zeros(iWidth, iHeight);
RawBack = zeros(iWidth, iHeight);
Back = zeros(iWidth, iHeight);
RawAlpha = zeros(iWidth, iHeight);
InnerBW = zeros(iWidth, iHeight);
OuterBW = zeros(iWidth, iHeight);
UR = zeros(iWidth, iHeight);

% initialization
for i=1:size(IT,1)
    for j=1:size(IT,2)
        if IT(i,j)==0
            RawAlpha(i,j)=0;
            Back(i,j)=GrayI(i,j);
            InnerBW(i,j)=0;
            OuterBW(i,j)=0;
            UR(i,j)=0;
        end;
        if IT(i,j)==128
            RawAlpha(i,j)=0.5;
            InnerBW(i,j)=0;
            OuterBW(i,j)=1;
            UR(i,j)=1;
        end;
        if IT(i,j)==255
            RawAlpha(i,j)=1;
            Fore(i,j)=GrayI(i,j);
            InnerBW(i,j)=1;
            OuterBW(i,j)=1;
            UR(i,j)=0;
        end;
    end;
end;

% imshow(uint8(Fore),'parent',handles.axes3);
% imshow(uint8(Back),'parent',handles.axes3);

for i = 1:iWidth
    for j = 1:iHeight
        if RawAlpha(i,j) == 1
            RawFore(i,j) = Fore(i,j);
        end;
        if RawAlpha(i,j) == 0
            RawBack(i,j) = Back(i,j);
        end;
        % if the pixel is in the Unknown Region...
        if RawAlpha(i,j) == 0.5
            % find its nearest foreground pixel through a r * r mask,
            % FindFore and FindBack.
            r = 1;
            while 1
                FindFore = InnerBW(max(i-r,1):min(i+r,iWidth),max(j-r,1):min(j+r,iHeight));
                if size(find(FindFore),1)
                    % GrayI1 and j1 are the row and column index respectively
                    % that the entris of the mask are nonzeros
                    [i1, j1] = find(FindFore);
                    
                    % find the nearest pixel in the mask
                    big=10000;
                    nearesti=i;
                    nearestj=j;
                    totalpixel=size(i1,1);
                    p=1;q=1;
                    while totalpixel>0
                        if abs(i1(p)-r-1)+abs(j1(q)-r-1)<big
                            nearesti=i1(p);
                            nearestj=j1(q);
                            big=abs(i1(p)-r-1)+abs(j1(q)-r-1);
                        end;
                        totalpixel=totalpixel-1;
                        p=p+1;
                        q=q+1;
                    end;
                    nearesti=nearesti+max(i-r,1)-1;
                    nearestj=nearestj+max(j-r,1)-1;
                    
                    RawFore(i,j)=Fore(nearesti,nearestj);
                    break;
                else
                    r=r+1;
                end;
            end;
            % find its nearest background pixel through a r * r mask
            r=1;
            while 1
                FindBack = ~OuterBW(max(i-r,1):min(i+r,iWidth),max(j-r,1):min(j+r,iHeight));
                if size(find(FindBack),1)
                    [i1, j1] = find(FindBack);
                    
                    big=10000;
                    nearesti=i;
                    nearestj=j;
                    totalpixel=size(i1,1);
                    p=1;q=1;
                    while totalpixel>0
                        if abs(i1(p)-r-1)+abs(j1(q)-r-1)<big
                            nearesti=i1(p);
                            nearestj=j1(q);
                            big=abs(i1(p)-r-1)+abs(j1(q)-r-1);
                        end;
                        totalpixel=totalpixel-1;
                        p=p+1;
                        q=q+1;
                    end;
                    nearesti=nearesti+max(i-r,1)-1;
                    nearestj=nearestj+max(j-r,1)-1;
                    
                    RawBack(i,j)=Back(nearesti,nearestj);
                    break;
                else
                    r=r+1;
                end;
            end;
        end;
    end;
end;
imshow(uint8(RawFore),'parent',handles.axes3);
imshow(uint8(RawBack),'parent',handles.axes3);

% when to stop iteration in Poisson Equation solving
Threshold = 0.1;

% --- Calculating Alpha ---
again = 1;
OldAlpha = RawAlpha;
NewAlpha = RawAlpha;
LastAlpha = RawAlpha;

% iteration times
count = 1;

while again == 1
    again = 0;
    FMinusB = RawFore - RawBack;
    GrayI1=double(GrayI);

    % Gauss Filter to smooth F-B.
    for i = 2:(size(FMinusB,1)-1)
        for j = 2:(size(FMinusB,2)-1)
            FMinusB(i,j) = (FMinusB(i-1,j-1)+FMinusB(i-1,j+1)+FMinusB(i+1,j-1)+FMinusB(i+1,j+1)+2*(FMinusB(i,j-1)+FMinusB(i-1,j)+FMinusB(i,j+1)+FMinusB(i+1,j))+4*FMinusB(i,j))/16;
        end;
    end;

    % In case divided by 0.
    for i = 1:size(FMinusB,1)
        for j = 1: size(FMinusB,2)
            if FMinusB(i,j) == 0
                FMinusB(i,j) = 1;
            end;
        end;
    end;

    iter = 0;
    while 1
        for i=2:(size(OldAlpha,1)-1)
            for j=2:(size(OldAlpha,2)-1)
                if UR(i,j)==1
%                     NewAlpha(i,j) = OldAlpha(i,j);
%                 if OldAlpha(i,j)>0 && OldAlpha(i,j)<1
                    % Solve Poisson Equation using Gauss-Sidel Method. See
                    % "Numerical Solution of Poisson Equation.pdf" for detail.
                    DivI = ((GrayI1(i+1,j) + GrayI1(i-1,j) - 2 * GrayI1(i,j)) * FMinusB(i,j) - (GrayI1(i+1,j) - GrayI1(i,j)) * (FMinusB(i+1,j) - FMinusB(i,j)))/(FMinusB(i,j) * FMinusB(i,j));
                    DivJ = ((GrayI1(i,j+1) + GrayI1(i,j-1) - 2 * GrayI1(i,j)) * FMinusB(i,j) - (GrayI1(i,j+1) - GrayI1(i,j)) * (FMinusB(i,j+1) - FMinusB(i,j)))/(FMinusB(i,j) * FMinusB(i,j));
                    Div = DivI + DivJ;
                    NewAlpha(i,j) = (OldAlpha(i+1,j) + NewAlpha(i-1,j) + OldAlpha(i,j+1) + NewAlpha(i,j-1) - Div) / 4;
%                     if NewAlpha(i,j)<0
%                         NewAlpha(i,j)=0;
%                     end;
%                     if NewAlpha(i,j)>1
%                         NewAlpha(i,j)=1;
%                     end;
                end;
            end;
        end;
        DifferenceAlpha = abs(NewAlpha - OldAlpha);
        residual=sum(sum(DifferenceAlpha));
        OldAlpha = NewAlpha;
        if  residual < Threshold
            break;
        end;
        iter = iter + 1;
        if iter > 50
            break;
        end;
    end;

    % Foreground and Background refinement
    % After 1 iteration, if there're elements in OmegaF or OmegaB (means Foreground or Background can be updated),
    % then do the next iteration. Otherwise don't.
    for i = 1:size(OldAlpha,1)
        for j = 1:size(OldAlpha,2)
            if OldAlpha(i,j) < 0.05 && UR(i,j)==1
                OldAlpha(i,j) = 0;
                RawBack(i,j) = Back(i,j);
                UR(i,j)=0;
                OuterBW(i,j)=0;
                InnerBW(i,j)=0;
                again=1;
            end;
            if OldAlpha(i,j) > 0.95 && UR(i,j)==1
                OldAlpha(i,j) = 1;
                RawFore(i,j) = Fore(i,j);
                UR(i,j)=0;
                OuterBW(i,j)=0;
                InnerBW(i,j)=1;
                again=1;
            end;
        end;
    end;
    
    imshow(uint8(OldAlpha*255),'parent',handles.axes3);
    
    if again == 1
        for i=1:size(OldAlpha,1)
            for j=1:size(OldAlpha,2)
                if UR(i,j)==1
                    r = 1;
                    while 1
                        FindFore = InnerBW(max(i-r,1):min(i+r,iWidth),max(j-r,1):min(j+r,iHeight));
                        if size(find(FindFore),1)
                            % i1 and j1 are the row and column index respectively
                            % that the entris of the mask are nonzeros
                            [i1, j1] = find(FindFore);

                            % find the nearest pixel in the mask
                            big=10000;
                            nearesti=i;
                            nearestj=j;
                            totalpixel=size(i1,1);
                            p=1;q=1;
                            while totalpixel>0
                                if abs(i1(p)-r-1)+abs(j1(q)-r-1)<big
                                    nearesti=i1(p);
                                    nearestj=j1(q);
                                    big=abs(i1(p)-r-1)+abs(j1(q)-r-1);
                                end;
                                totalpixel=totalpixel-1;
                                p=p+1;
                                q=q+1;
                            end;
                            nearesti=nearesti+max(i-r,1)-1;
                            nearestj=nearestj+max(j-r,1)-1;

                            RawFore(i,j)=Fore(nearesti,nearestj);
                            break;
                        else
                            r=r+1;
                        end;
                    end;
                    % find its nearest background pixel through a r * r mask
                    r=1;
                    while 1
                        FindBack = ~OuterBW(max(i-r,1):min(i+r,iWidth),max(j-r,1):min(j+r,iHeight));
                        if size(find(FindBack),1)
                            [i1, j1] = find(FindBack);

                            % find the nearest pixel in the mask
                            big=10000;
                            nearesti=i;
                            nearestj=j;
                            totalpixel=size(i1,1);
                            p=1;q=1;
                            while totalpixel>0
                                if abs(i1(p)-r-1)+abs(j1(q)-r-1)<big
                                    nearesti=i1(p);
                                    nearestj=j1(q);
                                    big=abs(i1(p)-r-1)+abs(j1(q)-r-1);
                                end;
                                totalpixel=totalpixel-1;
                                p=p+1;
                                q=q+1;
                            end;
                            nearesti=nearesti+max(i-r,1)-1;
                            nearestj=nearestj+max(j-r,1)-1;

                            RawBack(i,j)=Back(nearesti,nearestj);
                            break;
                        else
                            r=r+1;
                        end;
                    end;
                end;
            end;
        end;
    end;

    % If change in alpha is suffciently small, don't do the next iteration.
    AlphaChange = abs(OldAlpha - LastAlpha);
    if sum(sum(AlphaChange)) < 0.01
        break;
    end;
    LastAlpha = OldAlpha;
    
    if count > 50
        break;
    end;
    count = count + 1;
end;

fprintf('%i\n',count);

imshow(uint8(OldAlpha*255),'parent',handles.axes3);
imwrite(OldAlpha, 'Alpha.bmp', 'bmp');

handles.NewAlpha = OldAlpha;
guidata(hObject, handles);