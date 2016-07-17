function varargout = NewGUI(varargin)
% NEWGUI MATLAB code for NewGUI.fig
%      NEWGUI, by itself, creates a new NEWGUI or raises the existing
%      singleton*.
%
%      H = NEWGUI returns the handle to a new NEWGUI or the handle to
%      the existing singleton*.
%
%      NEWGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEWGUI.M with the given input arguments.
%
%      NEWGUI('Property','Value',...) creates a new NEWGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NewGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NewGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NewGUI

% Last Modified by GUIDE v2.5 27-Apr-2011 20:52:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NewGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NewGUI_OutputFcn, ...
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


% --- Executes just before NewGUI is made visible.
function NewGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NewGUI (see VARARGIN)

% Choose default command line output for NewGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NewGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NewGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function [foreColor, backColor] = findNearestPixel(row, col, srcImg, Trimap)
[rows, cols] = size(srcImg);

r = 1;
while 1
    mask = Trimap(max(row-r,1):min(row+r,rows), max(col-r,1):min(col+r,cols));
    [rows_m, cols_m] = size(mask);
    xIndex = rows_m / 2 + 1;
    yIndex = cols_m / 2 + 1;
    dist = 10000;
    x = 0;
    y = 0;
    for col_m = 1 : cols_m
        for row_m = 1 : rows_m
            if mask(row_m, col_m) == 255
                temp = abs(row_m-xIndex) + abs(col_m-yIndex);
                if  temp < dist && temp > 0
                    dist = temp;
                    x = row_m;
                    y = col_m;
                end;
            end;
        end;
    end;
    if x > 0 || y > 0
        break;
    else
        r = r + 1;
    end;
end;
foreColor = srcImg(x, y, :);

r = 1;
while 1
    mask = Trimap(max(row-r,1):min(row+r,rows), max(col-r,1):min(col+r,cols));
    [rows_m, cols_m] = size(mask);
    xIndex = rows_m / 2 + 1;
    yIndex = cols_m / 2 + 1;
    dist = 10000;
    x = 0;
    y = 0;
    for col_m = 1 : cols_m
        for row_m = 1 : rows_m
            if mask(row_m, col_m) == 0
                temp = abs(row_m-xIndex) + abs(col_m-yIndex);
                if  temp < dist && temp > 0
                    dist = temp;
                    x = row_m;
                    y = col_m;
                end;
            end;
        end;
    end;
    if x > 0 || y > 0
        break;
    else
        r = r + 1;
    end;
end;
backColor = srcImg(x, y, :);

% --- Executes on button press in LoadSource.
function LoadSource_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [SourceFileName,SourcePathName] = uigetfile('*.jpg;*.bmp;*.png','Select the source image');
% filePath = strcat(SourcePathName, SourceFileName);

% TODO:
filePath = 'C:\Users\luofucong\Documents\MATLAB\MyMatting\NewGUI\cat.jpg';

handles.MyImage = imread(filePath);
imshow(handles.MyImage,'InitialMagnification','fit');
guidata(hObject, handles);

% --- Global Poisson Matting Method begins.
A = handles.MyImage;

% TODO:
filePath = 'C:\Users\luofucong\Documents\MATLAB\MyMatting\NewGUI\Trimap.bmp';
Trimap = imread(filePath);

[rows, cols] = size(Trimap);
srcImg = A;
RawAlpha = zeros(rows, cols);
RawFore = zeros(rows, cols, 3);
RawBack = zeros(rows, cols, 3);

for col = 1 : cols
    for row = 1 : rows
        if Trimap(row, col) == 0
            RawAlpha(row, col) = 0;
            RawFore(row, col, :) = 0;
            RawBack(row, col, :) = srcImg(row, col, :);
        else if Trimap(row, col) == 128
                RawAlpha(row, col) = 0.5;
                [RawFore(row, col, :), RawBack(row, col, :)] = findNearestPixel(row, col, srcImg, Trimap);
            else
                RawAlpha(row, col) = 1;
                RawFore(row, col, :) = srcImg(row, col, :);
                RawBack(row, col, :) = 0;
            end;
        end;
    end;
end;

% % --- Trimap entering begins ---
% hold on;
% 
% OuterBW = roipoly(A);
% B(:,:,1) = immultiply(A(:,:,1),OuterBW);% B is the Foreground and Unknown Region.
% B(:,:,2) = immultiply(A(:,:,2),OuterBW);
% B(:,:,3) = immultiply(A(:,:,3),OuterBW);
% Back(:,:,1) = immultiply(A(:,:,1),~OuterBW);% Back is the Background.
% Back(:,:,2) = immultiply(A(:,:,2),~OuterBW);
% Back(:,:,3) = immultiply(A(:,:,3),~OuterBW);
% imshow(B);
% 
% InnerBW = roipoly(B);
% C(:,:,1) = immultiply(B(:,:,1),~InnerBW);% C is the Unknown Region.
% C(:,:,2) = immultiply(B(:,:,2),~InnerBW);
% C(:,:,3) = immultiply(B(:,:,3),~InnerBW);
% Fore(:,:,1) = immultiply(A(:,:,1),InnerBW);% Fore is the Foreground.
% Fore(:,:,2) = immultiply(A(:,:,2),InnerBW);
% Fore(:,:,3) = immultiply(A(:,:,3),InnerBW);
% imshow(C);
% 
% hold off;
% % --- Trimap entering ends ---
% 
% RawAlpha = (double(OuterBW)+double(InnerBW))/2;
% I = double(A);
% 
% % saving trimap
% Trimap=zeros(size(RawAlpha,1), size(RawAlpha,2));
% for i=1:size(RawAlpha,1)
%     for j=1:size(RawAlpha,2)
%         if RawAlpha(i,j)==0
%             Trimap(i,j)=0;
%         end;
%         if RawAlpha(i,j)==0.5
%             Trimap(i,j)=128;
%         end;
%         if RawAlpha(i,j)==1
%             Trimap(i,j)=255;
%         end;
%     end;
% end;
% imwrite(uint8(Trimap), 'C:\Users\luofucong\Downloads\Trimap.bmp', 'bmp');

% Get Foreground and Background colors.
% RawFore = zeros(size(RawAlpha,1), size(RawAlpha,2), 3);
% RawBack = zeros(size(RawAlpha,1), size(RawAlpha,2), 3);

% for i = 1:size(RawAlpha,1)
%     for j = 1:size(RawAlpha,2)
%         if RawAlpha(i,j) == 1
%             RawFore(i,j,:) = Fore(i,j,:);
%         end;
%         if RawAlpha(i,j) == 0
%             RawBack(i,j,:) = Back(i,j,:);
%         end;
%         % if the pixel is in the Unknown Region...
%         if RawAlpha(i,j) == 0.5
%             % find its nearest foreground pixel through a r * r mask,
%             % FindFore and FindBack.
%             r = 1;
%             while 1
%                 FindFore = InnerBW(max(i-r,1):min(i+r,size(A,1)),max(j-r,1):min(j+r,size(A,2)));
%                 if size(find(FindFore),1)
%                     % i1 and j1 are the row and column index respectively
%                     % that the entris of the mask are nonzeros
%                     [i1, j1] = find(FindFore);
%                     % i2 and j2 are the row and column index respectively
%                     % of the real image
%                     i2 = i1 + max(i-r,1) - 1;
%                     j2 = j1 + max(j-r,1) - 1;
%                     
%                     AllNonzeros = double(Fore(i2,j2,:));
%                     % the entries' color are in the diagonal of AllNonzeros
%                     Color(:,1) = diag(AllNonzeros(:,:,1));
%                     Color(:,2) = diag(AllNonzeros(:,:,2));
%                     Color(:,3) = diag(AllNonzeros(:,:,3));
%                     
%                     if size(Color,1) == 1
%                         % if grayscale image (the row of Color is 1)
%                         RawFore(i,j,:) = Color;
%                     else
%                         % if color image
%                         % average of R,G,B intensity
%                         RawFore(i,j,:) = sum(Color) / size(Color,1);
%                     end;
%                     clear Color;
%                     break;
%                 else
%                     r=r+1;
%                 end;
%             end;
%             % find its nearest background pixel through a r * r mask
%             r=1;
%             while 1
%                 FindBack = ~OuterBW(max(i-r,1):min(i+r,size(A,1)),max(j-r,1):min(j+r,size(A,2)));
%                 if size(find(FindBack),1)
%                     [i1, j1] = find(FindBack);
%                     i2 = i1 + max(i-r,1) - 1;
%                     j2 = j1 + max(j-r,1) - 1;
%                     AllNonzeros = double(Back(i2,j2,:));
%                     Color(:,1) = diag(AllNonzeros(:,:,1));
%                     Color(:,2) = diag(AllNonzeros(:,:,2));
%                     Color(:,3) = diag(AllNonzeros(:,:,3));
%                     if size(Color,1) == 1
%                         RawBack(i,j,:) = Color;
%                     else
%                         RawBack(i,j,:) = sum(Color) / size(Color,1);
%                     end;
%                     clear Color;
%                     break;
%                 else
%                     r=r+1;
%                 end;
%             end;
%         end;
%     end;
% end;

% when to stop iteration in Poisson Equation solving
Threshold = 1;

% --- Calculating Alpha ---
again = 1;
OldAlpha = RawAlpha;
NewAlpha = RawAlpha;

% iteration times
count = 1;

while again == 1
    again = 0;
    FMinusB = RawFore - RawBack;

    % Calculate Grayscale Channel.
%     GrayI = 0.3 * I(:,:,1) + 0.59 * I(:,:,2) + 0.11 * I(:,:,3);
%     FMinusB1 = 0.3*FMinusB(:,:,1) + 0.59*FMinusB(:,:,2) + 0.11*FMinusB(:,:,3);
    GrayI = srcImg(:, :, 2);
    FMinusB1 = FMinusB(:, :, 2);

    % Gauss Filter to smooth F-B.
    for i = 2:(size(FMinusB1,1)-1)
        for j = 2:(size(FMinusB1,2)-1)
            FMinusB1(i,j) = (FMinusB1(i-1,j-1)+FMinusB1(i-1,j+1)+FMinusB1(i+1,j-1)+FMinusB1(i+1,j+1)+2*(FMinusB1(i,j-1)+FMinusB1(i-1,j)+FMinusB1(i,j+1)+FMinusB1(i+1,j))+4*FMinusB1(i,j))/16;
        end;
    end;

    % In case divided by 0.
    for i = 1:size(FMinusB1,1)
        for j = 1: size(FMinusB1,2)
            if FMinusB1(i,j) == 0
                FMinusB1(i,j) = 1;
            end;
        end;
    end;

    while 1
        for i=1:size(OldAlpha,1)
            for j=1:size(OldAlpha,2)
                NewAlpha(i,j) = OldAlpha(i,j);
                if OldAlpha(i,j)>0 && OldAlpha(i,j)<1
                    % Solve Poisson Equation using Gauss-Sidel Method. See
                    % "Numerical Solution of Poisson Equation.pdf" for detail.
                    DivI = ((GrayI(i+1,j) + GrayI(i-1,j) - 2 * GrayI(i,j)) * FMinusB1(i,j) - (GrayI(i+1,j) - GrayI(i,j)) * (FMinusB1(i+1,j) - FMinusB1(i,j)))/(FMinusB1(i,j) * FMinusB1(i,j));
                    DivJ = ((GrayI(i,j+1) + GrayI(i,j-1) - 2 * GrayI(i,j)) * FMinusB1(i,j) - (GrayI(i,j+1) - GrayI(i,j)) * (FMinusB1(i,j+1) - FMinusB1(i,j)))/(FMinusB1(i,j) * FMinusB1(i,j));
                    Div = DivI + DivJ;
                    NewAlpha(i,j) = (OldAlpha(i+1,j) + NewAlpha(i-1,j) + OldAlpha(i,j+1) + NewAlpha(i,j-1) - Div) / 4;
%                     fprintf('%f\t %f\n',NewAlpha(i,j),Div);

                    if NewAlpha(i,j)<0
                        NewAlpha(i,j)=0;
                    end;
                    if NewAlpha(i,j)>1
                        NewAlpha(i,j)=1;
                    end;
                end;
            end;
        end;
        DifferenceAlpha = abs(NewAlpha - OldAlpha);
        OldAlpha = NewAlpha;
        if sum(sum(DifferenceAlpha)) < Threshold
            break;
        end;                
    end;

    % Foreground and Background refinement
    for i = 1:size(OldAlpha,1)
        for j = 1:size(OldAlpha,2)
            if OldAlpha(i,j) < 0.05 && OldAlpha(i,j) > 0
                OldAlpha(i,j) = 0;
                RawBack(i,j,:) = Back(i,j,:);
                again = 1;
            end;
            if OldAlpha(i,j) > 0.95 && OldAlpha(i,j) < 1
                OldAlpha(i,j) = 1;
                RawFore(i,j,:) = Fore(i,j,:);
                again = 1;
            end;
        end;
    end;
    count = count + 1;
end;

fprintf('%i',count);

OldAlpha = uint8(OldAlpha*255);
imshow(OldAlpha);
imwrite(OldAlpha, 'C:\Users\luofucong\Downloads\Alpha.bmp', 'bmp');

handles.NewAlpha = OldAlpha;
guidata(hObject, handles);


% --- Executes on button press in Composite.
function Composite_Callback(hObject, eventdata, handles)
% hObject    handle to Composite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[SourceFileName,SourcePathName] = uigetfile('*.jpg;*.bmp;*.png','Select the source image');
filePath = strcat(SourcePathName, SourceFileName);
handles.backImg = imread(filePath);
handles.alphaImg = imread('C:\Users\luofucong\Downloads\Alpha.bmp');
handles.sourceImg = imread('C:\Users\luofucong\Downloads\source.bmp');
guidata(hObject, handles);

backImg = handles.backImg;
alphaImg = handles.alphaImg;
mask = double(alphaImg)/255;
sourceImg = handles.sourceImg;
compositeImg = backImg;

for i=1:size(mask,1)
    for j=1:size(mask,2)
        compositeImg(i,j,:) = backImg(i,j,:) * (1-mask(i,j)) + sourceImg(i,j,:) * mask(i,j);
    end;
end;
imshow(compositeImg,'InitialMagnification','fit');
imwrite(compositeImg, 'C:\Users\luofucong\Downloads\composite.bmp', 'bmp');
guidata(hObject, handles);
