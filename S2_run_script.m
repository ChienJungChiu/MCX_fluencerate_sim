clc;clear;close all;

%% clear the stop flag
to_save=0;
save('stop_flag.txt','to_save','-ascii','-tabs');

%% simulate 

% each line for one source type
S2_fun_main_simulation_literature_param_anglepattern(1,6,1:64);
S2_fun_main_simulation_literature_param_anglepattern(13,6,1:64);

% can run for part of the subject or part of the OPs
%{
S2_fun_main_simulation_literature_param_anglepattern(3,1:4,1:64);
S2_fun_main_simulation_literature_param_anglepattern(3,1:9,30:60);
%}