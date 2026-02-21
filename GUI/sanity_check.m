%% Sanity check - plot all calcium traces in one recording
% The script plots all ROIs (calcium traces) from one record
% Load file after preprocessing and classifiction
% This script should be in the same folder with the 'analys_oscillation' latest function

% Load data
% change file name here
data_mat=load('file_name');
n_cells=size(data_mat.ocillatory.Oscillatory_sig,2);
time_vec=data_mat.time_vec;

% Parameters
dt=time_vec(2)-time_vec(1); %frame rate
trim=10; %DEAFULT
glut_injection_time = 120; %DEAFULT
glut_injection_acc=glut_injection_time+dt-trim;
glut_inj_idx = round(glut_injection_acc/dt);

% Save PNG
% this section saves plots of each ROI in different PNG file in the path
for i=1:n_cells
    f=figure('visible','off');
    sig=data_mat.ocillatory.Oscillatory_sig(:,i);
    % make sure to update to the latest 'analys_oscillation' function
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = analys_oscillations_030724(sig,time_vec,glut_inj_idx,1,i);
    % change file name here
    saveas(f,sprintf('file_name_ROI%d.png',i))
end
close all

