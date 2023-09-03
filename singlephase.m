function varargout = singlephase(varargin)
% SINGLEPHASE MATLAB code for singlephase.fig
%      SINGLEPHASE, by itself, creates a new SINGLEPHASE or raises the existing
%      singleton*.
%
%      H = SINGLEPHASE returns the handle to a new SINGLEPHASE or the handle to
%      the existing singleton*.
%
%      SINGLEPHASE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLEPHASE.M with the given input arguments.
%
%      SINGLEPHASE('Property','Value',...) creates a new SINGLEPHASE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before singlephase_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to singlephase_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help singlephase

% Last Modified by GUIDE v2.5 23-Dec-2022 16:28:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @singlephase_OpeningFcn, ...
                   'gui_OutputFcn',  @singlephase_OutputFcn, ...
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


% --- Executes just before singlephase is made visible.
function singlephase_OpeningFcn(hObject, eventdata, handles, varargin)


% Choose default command line output for singlephase
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


function varargout = singlephase_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

%INPUTS
v = str2double( get(handles.v, 'string'));
p = str2double( get(handles.p, 'string'));
f = str2double( get(handles.f, 'string'));
pole = str2double( get(handles.pole, 'string'));
n = str2double( get(handles.n, 'string'));
vdc_M = str2double( get(handles.vdc_M, 'string'));
idc_M = str2double( get(handles.idc_M, 'string'));
vdc_A = str2double( get(handles.vdc_A, 'string'));
idc_A = str2double( get(handles.idc_A, 'string'));
vlr_M = str2double( get(handles.vlr_M, 'string'));
ilr_M = str2double( get(handles.ilr_M, 'string'));
plr_M = str2double( get(handles.plr_M, 'string'));
vlr_A = str2double( get(handles.vlr_A, 'string'));
ilr_A = str2double( get(handles.ilr_A, 'string'));
plr_A = str2double( get(handles.plr_A, 'string'));
Vnl = str2double( get(handles.Vnl, 'string'));
Inl = str2double( get(handles.Inl, 'string'));
Pnl = str2double( get(handles.Pnl, 'string'));
ns = 2*60*f/pole;
s = (ns - n) / ns;
%DC test
R1m = vdc_M/idc_M;                                          %R1 main winding
R1a = vdc_A/idc_A;                                          %R1 aux winding
%LRT
Rlr_m = plr_M/ilr_M^2;                                        %R main
R2m = Rlr_m - R1m;                                            %R2' main
xlr_m = sqrt((vlr_M/ilr_M)^2 - Rlr_m^2);                      %X main
x1m = xlr_m/2;
x2m = xlr_m/2;
Rlr_a = plr_A/ilr_A^2;                                        %R aux
R2a = Rlr_a - R1a;                                            %R2' aux
xlr_a = sqrt((vlr_A/ilr_A)^2 - Rlr_a^2);                      %X aux
x1a = xlr_a/2;
x2a = xlr_a/2;

Zm = Rlr_m + xlr_m*1i;                            %standstill main impedance
Za = Rlr_a + xlr_a*1i;                            %standstill aux impedance

Effective_turns_ratio = sqrt(R2a/R2m);
ETR = Effective_turns_ratio;

%NLT
Rnl = R1m + R2m/4;
Prot = Pnl - (Inl^2)*Rnl;                      %Rotational losses
Xnl = sqrt((Vnl/Inl)^2 - Rnl^2);
Xm = 2*(Xnl - x1m - x2m/2);                    %leakage
%Resistance starting
R_add = 2*Rlr_a;
Im = v / Zm;                            %Main winding current
Ia_R = v/(Za + R_add);
Ist_R = Im + Ia_R;                      %Starting current w/ added resistance
a1 = abs((angle(Im) - angle(Ia_R)))*(180/pi);    %ALPHA R
Tst_R = ((2*ETR*R2m)/(ns*pi/30))*abs(Im)*abs(Ia_R)*sind(a1);

%Capacitor starting
XC_add = 2i*xlr_a;                     %added capacitance
Ia_C = v/(Za - XC_add);
Ist_C = Im + Ia_C;                     %Starting current w/ added capacitor
a2 = abs((angle(Im) - angle(Ia_C)))*(180/pi);        %ALPHA C
Tst_C = ((2*ETR*R2m)/(ns*pi/30))*abs(Im)*abs(Ia_C)*sind(a2);

%Running operation
Zf = 0.5i * (R2m/s + x2m*1i) * Xm / (R2m/s + x2m*1i +Xm*1i);       %Z forward
Zb = 0.5i * (R2m/(2-s) + x2m*1i)*Xm / (R2m/(2-s) + x2m*1i +Xm*1i); %Z backward
I_motor = v / (R1m + x1m*1i + Zf + Zb);      %Motor current
Pcu1 = (abs(I_motor)^2) * R1m;
Pgf = (abs(I_motor)^2) * real(Zf);
Pgb = (abs(I_motor)^2) * real(Zb);
Pcu2_f = Pgf * s;
Pcu2_b = Pgb * (2-s);
Pcu2 = Pcu2_f + Pcu2_b;
P_losses = Pcu1 + Pcu2 + Prot;              %Total power losses
Pd = (1-s) * (Pgf - Pgb);                 %Developed power at desired speed
P_output = Pd - Prot;
P_input = P_output + P_losses;
Power_factor = cosd(atand(imag(I_motor)/real(I_motor)));
Motor_Efficiency = (P_output / P_input) * 100;

set(handles.R1m, 'string',num2str(R1m));
set(handles.R2m, 'string',num2str(R2m));
set(handles.x1m, 'string',num2str(x1m));
set(handles.x2m, 'string',num2str(x2m));
set(handles.R1a, 'string',num2str(R1a));
set(handles.R2a, 'string',num2str(R2a));
set(handles.x1a, 'string',num2str(x1a));
set(handles.x2a, 'string',num2str(x2a));
set(handles.Zm, 'string',num2str(Zm));
set(handles.Za, 'string',num2str(Za));
set(handles.Effective_turns_ratio, 'string',num2str(Effective_turns_ratio));
set(handles.Ist_R, 'string',num2str(Ist_R));
set(handles.Tst_R, 'string',num2str(Tst_R));
set(handles.Ist_C, 'string',num2str(Ist_C));
set(handles.Tst_C, 'string',num2str(Tst_C));
set(handles.I_motor, 'string',num2str(I_motor));
set(handles.P_losses, 'string',num2str(P_losses));
set(handles.P_output, 'string',num2str(P_output));
set(handles.P_input, 'string',num2str(P_input));
set(handles.Power_factor, 'string',num2str(Power_factor));
set(handles.Motor_Efficiency, 'string',num2str(Motor_Efficiency));



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

% Clear inputs
q=0;
w= num2str(q);
set(handles.v, 'string',w);
set(handles.p, 'string',w);
set(handles.f, 'string',w);
set(handles.pole, 'string',w);
set(handles.n, 'string',w);
set(handles.vdc_M, 'string',w);
set(handles.idc_M, 'string',w);
set(handles.vdc_A, 'string',w);
set(handles.idc_A, 'string',w);
set(handles.vlr_M, 'string',w);
set(handles.ilr_M, 'string',w);
set(handles.plr_M, 'string',w);
set(handles.vlr_A, 'string',w);
set(handles.ilr_A, 'string',w);
set(handles.plr_A, 'string',w);
set(handles.Vnl, 'string',w);
set(handles.Inl, 'string',w);
set(handles.Pnl, 'string',w);
% Clear output
set(handles.R1m, 'string','');
set(handles.R2m, 'string','');
set(handles.x1m, 'string','');
set(handles.x2m, 'string','');
set(handles.R1a, 'string','');
set(handles.R2a, 'string','');
set(handles.x1a, 'string','');
set(handles.x2a, 'string','');
set(handles.Zm, 'string','');
set(handles.Za, 'string','');
set(handles.Effective_turns_ratio, 'string','');
set(handles.Ist_R, 'string','');
set(handles.Tst_R, 'string','');
set(handles.Ist_C, 'string','');
set(handles.Tst_C, 'string','');
set(handles.I_motor, 'string','');
set(handles.P_losses, 'string','');
set(handles.P_output, 'string','');
set(handles.P_input, 'string','');
set(handles.Power_factor, 'string','');
set(handles.Motor_Efficiency, 'string','');



function Vnl_Callback(hObject, eventdata, handles)

function Vnl_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Inl_Callback(hObject, eventdata, handles)

function Inl_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pnl_Callback(hObject, eventdata, handles)

function Pnl_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vlr_A_Callback(hObject, eventdata, handles)

function vlr_A_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ilr_A_Callback(hObject, eventdata, handles)

function ilr_A_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plr_A_Callback(hObject, eventdata, handles)

function plr_A_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vlr_M_Callback(hObject, eventdata, handles)

function vlr_M_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ilr_M_Callback(hObject, eventdata, handles)

function ilr_M_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plr_M_Callback(hObject, eventdata, handles)

function plr_M_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vdc_A_Callback(hObject, eventdata, handles)

function vdc_A_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function idc_A_Callback(hObject, eventdata, handles)

function idc_A_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vdc_M_Callback(hObject, eventdata, handles)

function vdc_M_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function idc_M_Callback(hObject, eventdata, handles)

function idc_M_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_Callback(hObject, eventdata, handles)

function f_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pole_Callback(hObject, eventdata, handles)

function pole_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_Callback(hObject, eventdata, handles)

function n_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_Callback(hObject, eventdata, handles)

function v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_Callback(hObject, eventdata, handles)

function p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
