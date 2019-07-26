function varargout = Browser(varargin)
% browser MATLAB code for browser.fig
%      BROWSER, by itself, creates a new BROWSER or raises the existing
%      singleton*.
%
%      H = BROWSER returns the handle to a new BROWSER or the handle to
%      the existing singleton*.
%
%      BROWSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BROWSER.M with the given input arguments.
%
%      BROWSER('Property','Value',...) creates a new BROWSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Browser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Browser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Browser

% Last Modified by GUIDE v2.5 07-Sep-2013 14:54:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Browser_OpeningFcn, ...
                   'gui_OutputFcn',  @Browser_OutputFcn, ...
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

% --- Executes just before Browser is made visible.
function Browser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Browser (see VARARGIN)

% Choose default command line output for Browser
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.opt_Combine, 'tooltipstring', sprintf('If there were multiple standards used in the experiment and\n you wish to treat all the standards the same, use this option'));
set(handles.opt_Separate, 'tooltipstring', sprintf('If there were multiple standards used in the experiment and\n you wish to treat all the standards separately, use this option.\nThis will create separate graphs for each standard.'));
% UIWAIT makes browser wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Executes on button press in opt_freq.
function opt_freq_Callback(hObject, eventdata, handles)
% hObject    handle to opt_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of opt_freq
end

% --- Executes on button press in opt_int.
function opt_int_Callback(hObject, eventdata, handles)
% hObject    handle to opt_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of opt_int
end

% --- Outputs from this function are returned to the command line.
function varargout = Browser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in btn_Browse.
function btn_Browse_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rawdata = [];

if get(handles.opt_LeftRight, 'value') == 1;
    Opts.resptype = 'LR';
else
    Opts.resptype = '12';
end
if get(handles.opt_Combine, 'value') == 1;
    Opts.Stands = 'combine';
else
    Opts.Stands = 'separate';
end
if get(handles.opt_freq, 'value') == 1;
    Opts.suff = 'Hz';
else
    Opts.suff = 'um';
end
paramcol = str2double(get(handles.txt_param, 'string'));


[filename, pathname, ~] = uigetfile('.csv','',['D:\Dynamic Playwave\Dynamic Load Playwave\Experiments\' datestr(now,29)], 'multiselect', 'on');
if iscell(filename)  %have multiple files been selected
    if get(handles.btn_new, 'value') == 1;
        ln = csvread([pathname filename{1}],1,0,[1,0,1,0]);  %read how many waveforms there are
        try
            for i = 1:numel(filename)   %concatenate files
                temp = csvread([pathname filename{i}],ln+6,0);
                rawdata = [rawdata; temp];
            end
        catch
            msgbox('error with file format. Try change File Format option.')
        end
    else
        ln = csvread([pathname filename{1}],1,0,[1,0,1,0]);  %read how many waveforms there are
        try
            for i = 1:numel(filename)   %concatenate files
                temp = csvread([pathname filename{i}],1,0);
                rawdata = [rawdata; temp];
            end
        catch
            msgbox('error with file format. Try change File Format option.')
        end
    end
    filename = [filename{1}(1:end-6) '_all' filename{1}(end-3:end)]; %just select first one for csv check below (and add 'all' for combined) - check probably not entirely necessary 

else
    if get(handles.btn_new, 'value') == 1;
        ln = csvread([pathname filename],1,0,[1,0,1,0]);   %read how many waveforms there are
        try
            rawdata = csvread([pathname filename],ln+6,0);
        catch
            msgbox('error with file format. Try change File Format option.')
        end
    else
        ln = csvread([pathname filename],1,0,[1,0,1,0]);   %read how many waveforms there are
        try
            rawdata = csvread([pathname filename],1,0);
        catch
            msgbox('error with file format. Try change File Format option.')
        end
    end
end
if filename ~= 0 
    if strcmp(filename(end-2:end),'csv') 
        id = strfind(pathname,filesep);
        stimParams = csvread(strcat(pathname(1:id(end-1)),'StimParams.csv'),4,0);
        %checktype(stimParams);  %maybe future use for checking if selected options
        %are suitable
        [DATA,bandNo,mydata] = FitPsychFn9(pathname, rawdata, Opts, stimParams, paramcol);  
    end
    plotting(bandNo,DATA,mydata,filename,Opts)
end
end    

function plotting(bandNo,DATA,mydata,fname,Opts)
%create axes for plotting
for i = 1:bandNo
    name = fname(21:end-4);
%     ax(i) = axes('parent',handles.pan_plot, 'position', [);
    h = figure('name', [name,'_',num2str(DATA.stndf(i)),'Hz'], 'position', [500 300 800 500]);
    subplot(1,2,1)
    plot(DATA.compfreq{i},DATA.pCompFaster_adj{i},'k.','markersize',10);
%     set(handles.axes1,'YLim',[])
    hold on;
    %plot fitted function
    plot(DATA.freqFine{i},DATA.curveFine{i},'r-','linewidth',1);
    hold off;
    xlabel(['Comparison Parameter (' Opts.suff ')']);
    ylabel('Proportion Comparison Greater');
    axis([min(DATA.compfreq{i})-std(DATA.compfreq{i})/3 max(DATA.compfreq{i})+std(DATA.compfreq{i})/3 -0.1 1.1])
    subplot(1,2,2)
    axis off
    text(-0.1,1,['Std = ' num2str(DATA.stndf(i)) ' ' Opts.suff])
    text(-0.1,0.9,['Band = ' num2str(i)])
    text(-0.1,0.8,['PSE = ' num2str(DATA.PSE(i)) ' ' Opts.suff])
    text(-0.1,0.7,['WeberLower = ' num2str(DATA.WeberLower(i))])
    text(-0.1,0.6,['WeberHigher = ' num2str(DATA.WeberHigher(i))])
    text(-0.2,0.5,'comp nCompGreater nTrials pCompGreater pCompGreater-adj logit')
    text(-0.1,0.3,num2str(mydata{i}))
    print(h,'-dtiffnocompression', strcat(DATA.figdir,name,'_band',num2str(i),'_',num2str(DATA.stndf(i)),'Hz.tif'));
end          
end

function checktype(stimParams)    



end


% --- Executes when selected object is changed in uipanel5.
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel5 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.opt_freq, 'value') == 1;
    set(handles.txt_param,'string','1')
else
    set(handles.txt_param,'string','3')
end
end

function txt_param_Callback(hObject, eventdata, handles)
% hObject    handle to txt_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_param as text
%        str2double(get(hObject,'String')) returns contents of txt_param as a double
end

% --- Executes during object creation, after setting all properties.
function txt_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
