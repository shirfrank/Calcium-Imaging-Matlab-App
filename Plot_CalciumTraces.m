%% Plot calcium traces


%% Zoom in - plot peaks of an ROI

%insert signal parameters here:
data_mat=load('test.mat');
data_struct = data_mat.myStruct; % Access the struct inside the MAT file
time_vec=data_struct.Time_vec;
cell_num=74
sig = data_struct(cell_num).Signal;


% Extract the signal for the specified cell
cell_num = 74; % Specify the cell number


max_sig=max(sig); %max signal amplitude
normsig = (sig-min(sig))/(max(sig)-min(sig));
O_duration_min = 21; %sec
O_duration_max = 180; %sec
O_prominence_min = 0.3*max_sig; 
O_dis_min = 60; %sec
O_amp_min = max(max_sig)*0.5;

figure()
findpeaks(sig,time_vec,'MinPeakWidth',O_duration_min,'MaxPeakWidth',O_duration_max,'MinPeakProminence',O_prominence_min,'MinPeakHeight',O_amp_min,'MinPeakDistance',O_dis_min,'Annotate','extents','WidthReference','halfprom');


%% Sanity check - plot all calcium traces from a record

%load data
data_mat=load('test.mat');
data_struct = data_mat.myStruct % Access the struct inside the MAT file
n_cells=size(data_struct,2);
time_vec=data_struct.Time_vec;

% Glutamate administration time
dt=time_vec(2)-time_vec(1); %frame rate
trim=10;%10 sec trim
glut_injection_time = 120; %DEAFULT:110, REAL; wt1-106 wt8-110 wt9-60 / mut1-60 mut1-2-120 mut3-120
glut_injection_acc=glut_injection_time+dt-trim;
glut_inj_idx = round(glut_injection_acc/dt);

for i=1:n_cells
    f=figure('visible','off');
    sig = data_struct(i).Signal;

    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = analys_oscillations(sig,time_vec,glut_inj_idx,1,i);

    %change file name here
    saveas(f,sprintf('file_name_ROI%d.png',i))
end