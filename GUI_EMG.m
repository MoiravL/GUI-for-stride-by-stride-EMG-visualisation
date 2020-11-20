function varargout = GUI_EMG(varargin)
% Gui for visualizing EMG data and filters
% Inputs:
%       1) Detrended EMG data (time x channels)
%       2) Vector of events (e.g. samples of left heelstrikes)
%       3) Filename (string for when you want to save NaN indices and new
%       EMG)
%       4) Sample frequency
%       
%   If "save EMG" is pressed, one saves the indices of the removed (made NaN) strides
%   as well as the edited EMG in a mat-file.


%   Authors: 
%   Moira van Leeuwen, Vrije Universiteit Amsterdam (september 2019)
        


% GUI_EMG MATLAB code for GUI_EMG.fig
%      GUI_EMG, by itself, creates a new GUI_EMG or raises the existing
%      singleton*.
%
%      H = GUI_EMG returns the handle to a new GUI_EMG or the handle to
%      the existing singleton*.
%
%      GUI_EMG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EMG.M with the given input arguments.
%
%      GUI_EMG('Property','Value',...) creates a new GUI_EMG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_EMG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_EMG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_EMG

% Last Modified by GUIDE v2.5 09-Aug-2019 15:54:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_EMG_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_EMG_OutputFcn, ...
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


% --- Executes just before GUI_EMG is made visible.
function GUI_EMG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_EMG (see VARARGIN)

emg_data = varargin{1};
handles.events   = varargin{2};

handles.filename = varargin{3};

handles.fs       = varargin{4};

handles.n_strides = 0;

handles.emg = emg_data;
handles.emg_current = emg_data;

handles.avstride_click = 3;

handles.emg_output = emg_data;
handles.emg_channel = 1;
handles.emg_window = handles.events(1):handles.events(end);%[1:length(emg_data)];

plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))


% Choose default command line output for GUI_EMG
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes GUI_EMG wait for user response (see UIRESUME)
% uiwait(handles.twotag);


% --- Outputs from this function are returned to the command line.
function [varargout] = GUI_EMG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.emg_output;


% --- Executes on button press in nxt.
function nxt_Callback(hObject, eventdata, handles)
% hObject    handle to nxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.n_strides == 0;
disp('error: only one window')
elseif handles.n_strides >0  && handles.n_window + handles.n_strides -1 <= length(handles.events)-handles.n_strides;
    if handles.n_strides ==1 && handles.n_window +2 <= length(handles.events)
        handles.n_window = handles.n_window+1;
    elseif handles.n_strides~= 1
    handles.n_window = handles.n_window + handles.n_strides -1; % 1 stride overlap for each window
    end
    handles.emg_window = [];
    handles.emg_window = handles.events(handles.n_window):handles.events(handles.n_window+handles.n_strides);
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
end
% Update handles structure
guidata(hObject, handles);

 


% --- Executes on button press in one.
function one_Callback(hObject, eventdata, handles)
% hObject    handle to one (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.n_window   = 1;
handles.emg_window = handles.events(1):handles.events(2);
handles.n_strides  = 1;
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.n_strides == 0;
disp('error: only one window')
elseif handles.n_strides > 0 && handles.n_window +1 > handles.n_strides;
    if handles.n_strides == 1 && handles.n_window > 1
      handles.n_window = handles.n_window-1;
    elseif handles.n_strides ~= 1
    handles.n_window = handles.n_window -handles.n_strides+1; % same here one step overlap
    end
    handles.emg_window = [];
    handles.emg_window = handles.events(handles.n_window):handles.events(handles.n_window+handles.n_strides);
end
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in five.
function five_Callback(hObject, eventdata, handles)
% hObject    handle to five (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.n_window   = 1;
handles.n_strides  = 5;

handles.emg_window = handles.events(1):handles.events(handles.n_strides)+1;

plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in twfive.
function twfive_Callback(hObject, eventdata, handles)
% hObject    handle to twfive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.n_window   = 1;
handles.n_strides  = 25;

handles.emg_window = handles.events(1):handles.events(handles.n_strides)+1;

plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in fifty.
function fifty_Callback(hObject, eventdata, handles)
% hObject    handle to fifty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.n_window   = 1;
handles.n_strides  = 50;

handles.emg_window = handles.events(1):handles.events(handles.n_strides)+1;

plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in hndrd.
function hndrd_Callback(hObject, eventdata, handles)
% hObject    handle to hndrd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.n_window   = 1;
handles.n_strides  = 100;

handles.emg_window = handles.events(1):handles.events(handles.n_strides)+1;

plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in all.
function all_Callback(hObject, eventdata, handles)
% hObject    handle to all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.n_strides = 0;
handles.emg_window = [handles.events(1):handles.events(end)];
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in nanwdw.
function nanwdw_Callback(hObject, eventdata, handles)
% hObject    handle to nanwdw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.emg_output(handles.emg_window,handles.emg_channel) = NaN;
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
guidata(hObject, handles);


% --- Executes on button press in Undo_nan.
function Undo_nan_Callback(hObject, eventdata, handles)
% hObject    handle to Undo_nan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.emg_output(handles.emg_window,handles.emg_channel) = handles.emg_current(handles.emg_window,handles.emg_channel)
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
guidata(hObject, handles);

% --- Executes on button press in ch_down.
function ch_down_Callback(hObject, eventdata, handles)
% hObject    handle to ch_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.emg_channel > 1

handles.emg_channel = handles.emg_channel -1;
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

set(handles.text_channel,'string',['EMG channel : ' num2str(handles.emg_channel)])

%disp(['EMG channel ' num2str(handles.emg_channel)])

else
    disp('This is EMG channel 1')
end

guidata(hObject, handles);


% --- Executes on button press in ch_up.
function ch_up_Callback(hObject, eventdata, handles)
% hObject    handle to ch_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.emg_channel < size(handles.emg,2)
    handles.emg_channel = handles.emg_channel + 1;
    
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

set(handles.text_channel,'string',['EMG channel : ' num2str(handles.emg_channel)])

 %   disp(['EMG channel ' num2str(handles.emg_channel)])
    
end
guidata(hObject, handles)
    

function min1stride_Callback(hObject, eventdata, handles)
if handles.n_strides ~= 0

if handles.n_window > 1
    handles.n_window = handles.n_window -1;
    handles.emg_window = handles.events(handles.n_window):handles.events(handles.n_window+handles.n_strides);
    
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

end

end

guidata(hObject, handles)


function plus1stride_Callback(hObject, eventdata, handles)
if handles.n_strides ~= 0

if handles.n_window +1 + handles.n_strides <= length(handles.events)
    handles.n_window = handles.n_window +1;
    handles.emg_window = handles.events(handles.n_window):handles.events(handles.n_window+handles.n_strides);
    
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

end
end

guidata(hObject, handles)


function min10stride_Callback(hObject, eventdata, handles)
if handles.n_strides ~= 0

if handles.n_window > 10
    handles.n_window = handles.n_window -10;
    handles.emg_window = handles.events(handles.n_window):handles.events(handles.n_window+handles.n_strides);
    
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

end
end

guidata(hObject, handles)


function plus10stride_Callback(hObject, eventdata, handles)

if handles.n_strides ~= 0

if handles.n_window +10 + handles.n_strides <= length(handles.events)
    handles.n_window = handles.n_window +10;
    handles.emg_window = handles.events(handles.n_window):handles.events(handles.n_window+handles.n_strides);
    
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

end
end

guidata(hObject, handles)


function nan_all_Callback(hObject,eventdata,handles)
handles.emg_output(handles.emg_window,:) = NaN;
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

guidata(hObject, handles)


function rest_all_Callback(hObject,eventdata,handles)
handles.emg_output(handles.emg_window,:) = handles.emg_current(handles.emg_window,:); 
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

guidata(hObject, handles)


function zoom2stride_Callback(hObject,eventdata,handles)
handles.n_strides = 1;

[x_stride,y_stride] = ginput(1);


distance =  abs(x_stride - handles.events);

[B, closest_2] = mink(distance,2);
handles.n_window = min(closest_2);
handles.emg_window = handles.events(handles.n_window):handles.events(handles.n_window+handles.n_strides);
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))


guidata(hObject, handles)

function raw_Callback(hObject,eventdata,handles)
current_nans = find(isnan(handles.emg_output));
handles.emg_output = handles.emg;
handles.emg_current = handles.emg_output;
handles.emg_output(current_nans) = NaN;
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
guidata(hObject, handles)

function highpass_Callback(hObject,eventdata,handles)
    
current_nans = find(isnan(handles.emg_output));

fc_high = str2num(get(handles.edit1,'string'))
   
Wn = fc_high/(handles.fs/2); % conform Rankin 20
[b,a] = butter(2,Wn,'high');

handles.emg_output = filtfilt(b,a,double(handles.emg_current)); % input is the current data without the NaNs
handles.emg_current   = handles.emg_output;
handles.emg_output(current_nans) = NaN;
plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
guidata(hObject, handles)


function rectify_Callback(hObject,eventdata,handles)

handles.emg_output = abs(handles.emg_output);
handles.emg_current = abs(handles.emg_current);

plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
guidata(hObject, handles)


function env_Callback(hObject,eventdata,handles)
    
current_nans = find(isnan(handles.emg_output));

fc_low = str2num(get(handles.edit2,'string'))
   
Wn = fc_low/(handles.fs/2); % conform Rankin 20
[c,d] = butter(2,Wn,'low');

handles.emg_output = filtfilt(c,d,double(handles.emg_current)); % input is the current output without the NaNs
handles.emg_current = handles.emg_output;
handles.emg_output(current_nans) = NaN;

plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
guidata(hObject, handles)


function saveEMG_Callback(hObject,eventdata,handles)
emg_save = handles.emg_output;
nan_values = find(isnan(handles.emg_output));
save(handles.filename,'emg_save','nan_values')


function power_Callback(hObject,eventdata,handles)
     [Pxx,F] = pwelch(handles.emg_current(handles.emg_window,handles.emg_channel)-mean(handles.emg_current(handles.emg_window,handles.emg_channel)),[],[],[],handles.fs); 
     plot(F,Pxx.^2)
     xlabel('frequency (Hz)')
     xlim([1 200])   
     
    function time_Callback(hObject,eventdata,handles)
         plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))

         
    function avstride_Callback(hObject,eventdata,handles)
        
        switch handles.avstride_click
            
            case 1
                handles.avstride_click = 2;
            case 2
                handles.avstride_click = 3;
            case 3
                handles.avstride_click = 1;  
        end
        
        
        if handles.n_strides == 0
        i_stride = 1:length(handles.events);
        nan_strides = zeros(1,length(handles.events)-1);
   
        else 
        i_stride = handles.n_window:handles.n_window+handles.n_strides;
        nan_strides = zeros(1,handles.n_strides);
        end
        
        strides = handles.events(i_stride);

       % newlength = round(mean(diff(strides)));

        newlength = 1000;
        for i = i_stride(1:end-1)
            
            temp = find(isnan(handles.emg_output(i:i+1,handles.emg_channel)));
            
            if length(temp)>0
                nan_strides(i) = 1;
            else
                nan_strides(i) = 0;
            end
            
            temp = [];
        end
      
        norm_strides =normalize_gait_events(handles.emg_current(:,handles.emg_channel),handles.events(i_stride(1:end-1)),handles.events(i_stride(2:end)),newlength);
        norm_data    = norm_strides.data;
        
        norm_strides =squeeze(norm_data);
        
        strides4nan = find(nan_strides==1);
        
        norm_strides(strides4nan,:) = NaN;
        av_stride = nanmean(norm_strides,1);
        plot(av_stride/max(av_stride))
        
        %xlabel('AVERAGE')
        
        if handles.avstride_click ==2
            plot(norm_strides'/max(av_stride))
            hold on
            plot(av_stride/max(av_stride),'k','linewidth',3)
            hold off
        end
        
        if handles.avstride_click == 3
            
             std_stride = nanstd(norm_strides,[],1);
            
            plotError(av_stride/max(av_stride),std_stride/max(av_stride))
            
            
        end
        
%         if handles.avstride_click == 4
%             plot(nanmedian(norm_strides,1)/max(nanmedian(norm_strides,1)))
%             xlabel('MEDIAN')
%         end
        
        guidata(hObject, handles)  
        
        function notch_Callback(hObject,eventdata,handles)
            
            current_nans = find(isnan(handles.emg_output));

           
            notch_freq = str2num(get(handles.valnotch,'string'));
     
            %[b,a]=butter(2,[(notch_freq-0.5)/(handles.fs/2),(notch_freq+0.5)/(handles.fs/2)],'stop')
            
            BW = 1/(handles.fs/2); %bandwidth of 1
            
            [b,a] = iircomb(round(handles.fs/notch_freq),BW,'notch');  

            
            handles.emg_output = filtfilt(b,a,double(handles.emg_current)); % input is the current data without the NaNs
            handles.emg_current   = handles.emg_output;
            handles.emg_output(current_nans) = NaN;
            plot(handles.emg_window,handles.emg_output(handles.emg_window,handles.emg_channel))
            guidata(hObject, handles)

            
              
         

 



    
