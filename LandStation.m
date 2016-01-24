function varargout = LandStation(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LandStation_OpeningFcn, ...
                   'gui_OutputFcn',  @LandStation_OutputFcn, ...
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

% --- Executes just before LandStation is made visible.
function LandStation_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

clc
imaqreset
imaqmem(1e11);
axes(handles.axes1);
axis off;

% --- Outputs from this function are returned to the command line.
function varargout = LandStation_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

clc

imaqreset
global video s
s = 1;
%video = videoinput('winvideo', 1, 'UYVY_720x480');
video = videoinput('winvideo', 2, 'MJPG_320x240');
set(video, 'FramesPerTrigger', Inf);
set(video, 'ReturnedColorspace', 'rgb');
video.FrameGrabInterval = 5;
start(video);
while( s == 1)
    imshow(getsnapshot(video))
    axes(handles.axes1);
    hold on
    if(s ~= 1)
        stop(video)
        flushdata(video)
    end
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

global s video
s = 0;
stop(video);
clc;
clear all

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

global s y

%s=0;
close(gcbf);
close
clear all
clc
      
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)

global video s

%s=0;
stop(video)
msgbox({'Ing. Armando René Narvaez Contreras';'Ing. Hector Manuel Zamora García de León';'Ing. Hilario Salomón Galván Guzmán';'Ing. Rodrigo Torres Arrazate'},'Students')
flushdata(video)

% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)

global ShapeF video suma ColorF

flushdata(video)

if (hObject == handles.Circle)
    ShapeF = 10;
    flushdata(video)
end

if (hObject == handles.Rectangle)
    ShapeF = 20;
    flushdata(video)
end

if (hObject == handles.Square)
    ShapeF = 30;
    flushdata(video)
end


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)

global video ShapeF ColorF suma s

flushdata(video)

if (hObject == handles.Red)
    ColorF = 1;
    flushdata(video);
end

if (hObject == handles.Green)
    ColorF = 2;
    flushdata(video);
end

if (hObject == handles.Blue)
    ColorF = 3;
    flushdata(video);
end

suma = ShapeF+ColorF;

%Rojo con circulo

while(suma == 11)
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,1), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    %diff_im = im2bw(diff_im,0.2); % Hue Webcam
    diff_im = im2bw(diff_im,0.072); % Hue
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 1
                plot(centroid(1),centroid(2),'wO');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end
end

%Rojo con cuadrado
        
while(suma == 21)
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,1), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    %diff_im = im2bw(diff_im,0.2); % Hue Webcam
    diff_im = im2bw(diff_im,0.072); % Hue
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 3
                plot(centroid(1),centroid(2),'wS');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end
end
    
%Rojo con rectangulo

while(suma == 31)
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,1), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    %diff_im = im2bw(diff_im,0.2); % Hue Webcam
    diff_im = im2bw(diff_im,0.072); % Hue
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 2
                plot(centroid(1),centroid(2),'wX');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end    
end

%Verde con circulo

while(suma == 12)
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,2), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    %diff_im = im2bw(diff_im,0.1);
    diff_im = im2bw(diff_im,0.032);    
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 1
                plot(centroid(1),centroid(2),'wO');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end    
end

%Verde con cuadrado
        
while(suma == 22)
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,2), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);    
    %diff_im = im2bw(diff_im,0.1);
    diff_im = im2bw(diff_im,0.032);    
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 3
                plot(centroid(1),centroid(2),'wS');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end    
end

%Verde con rectangulo

while(suma == 32)  
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,2), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    %diff_im = im2bw(diff_im,0.1);
    diff_im = im2bw(diff_im,0.032);    
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 2
                plot(centroid(1),centroid(2),'wX');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end    
end

%Azul con circulo

while(suma == 13)
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,3), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    diff_im = im2bw(diff_im,0.087);
    %diff_im = im2bw(diff_im,0.13);    
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 1
                    plot(centroid(1),centroid(2),'wO');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end    
end

%Azul con cuadrado
        
while(suma == 23)
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,3), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    diff_im = im2bw(diff_im,0.087);
    %diff_im = im2bw(diff_im,0.13);    
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 3
                plot(centroid(1),centroid(2),'wS');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end    
end

%Azulk mit Rectangulart

while(suma == 33)  
    dataR = getsnapshot(video); % Snapshot with camera
    diff_im = imsubtract(dataR(:,:,3), rgb2gray(dataR));
    diff_im = medfilt2(diff_im, [3 3]);
    diff_im = im2bw(diff_im,0.087);
    %diff_im = im2bw(diff_im,0.13);    
    diff_im = bwareaopen(diff_im,300); 
    bw = bwlabel(diff_im, 8);
    stats = regionprops(bw, 'BoundingBox', 'Centroid');       
    imshow(dataR) 
    axes(handles.axes1);
    hold on    
    for object = 1:length(stats)
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;
        plot(bc(1),bc(2), '-m+')
        a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
        set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    end
    hold off
    imwrite(getsnapshot(video),'Fig.bmp');
    RGB = imread('Fig.bmp');
    GRAY = rgb2gray(RGB);
    threshold = graythresh(GRAY);
    BW = im2bw(GRAY, threshold);
    BW = ~ BW;
    [B,L] = bwboundaries(BW, 'noholes');
    STATS = regionprops(L, 'all');
    imshow(RGB);
    axes(handles.axes1);
    hold on
    for i = 1 : length(STATS)
        W(i) = uint8(abs(STATS(i).BoundingBox(3)-STATS(i).BoundingBox(4)) < 0.1);
        W(i) = W(i) + 2 * uint8((STATS(i).Extent - 1) == 0 );
        centroid = STATS(i).Centroid;
        switch W(i)
            case 2
                plot(centroid(1),centroid(2),'wX');
        end
    end
    suma=ShapeF+ColorF;
    if(s ~= 1)
        stop(video)
        flushdata(video)
        suma=0;
    end    
end
