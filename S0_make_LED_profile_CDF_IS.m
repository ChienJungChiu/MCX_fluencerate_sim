%{
Make the CDF of LED profile for MCX anglepattern simulation

Benjamin Kao
Last update: 2020/12/02
%}

clc;clear;close all;

%% param
input_patten='LED_mouse_power_profile.txt'; % in (normalized power/sr) sr=steradian,the angle pattern of the source
num_angle=100000; % how many angle interval in the 0~90 degree range
num_CDF_IS_interval=10000; % how many interval in the 0~1 CDF range
do_plot=1;

%% main
if do_plot
    figure('Units','pixels','position',[0 0 1920 1080]);
    ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
end

% load the angle pattern
led_profile=load(input_patten);

% plot the angle pattern
if do_plot
    nexttile();
    plot(led_profile(:,1),led_profile(:,2));
    xlabel('angle');
    ylabel('power/unit solid angle');
    title('origianl angle profile');
end

% turn the origninal angle pattern (power/sr) to the probability of power comes out from certain angle (PDF)
output_angle=transpose(linspace(0,90,num_angle+1));
led_profile=interp1(led_profile(:,1),led_profile,output_angle,'pchip');
% multiply the angle by the circle perimeter
led_profile(:,2)=led_profile(:,2).*abs(sind(led_profile(:,1)+180/num_angle/2));

if do_plot
    nexttile();
    plot(led_profile(:,1),led_profile(:,2));
    xlabel('angle');
    ylabel('total power/each angle');
    title('power distribution PDF');
end

% calculate the CDF and normalize to 1
led_profile(:,2)=cumsum(led_profile(:,2));
led_profile(:,2)=led_profile(:,2)./max(led_profile(:,2));

if do_plot
    nexttile();
    plot(led_profile(:,1),led_profile(:,2));
    xlabel('angle');
    ylabel('CDF');
    title('power distribution CDF');
end

% turn the angle into radian
led_profile(:,1)=led_profile(:,1)./180*pi;
% make the CDF uniform distributed
output_CDF=transpose(linspace(0,1,num_CDF_IS_interval+1));
led_profile2=interp1(led_profile(:,2),led_profile(:,1),output_CDF);
led_profile2(isnan(led_profile2))=0;

% check
assert(issorted(led_profile2),'ERROR: CDF corresponding angle not increasing!');

save([strtok(input_patten,'.') '_CDF.txt'],'led_profile','-ascii','-tabs');
save([strtok(input_patten,'.') '_CDF_IS.txt'],'led_profile2','-ascii','-tabs');

if do_plot
    nexttile();
    plot(led_profile(:,1),led_profile(:,2),'-.',led_profile2,output_CDF,'--');
    xlabel('angle');
    ylabel('CDF');
    title('power distribution CDF');
    legend({'orig CDF','interp CDF'},'Location','best');
    title(ti,strrep(strtok(input_patten,'.'),'_',' '));
    saveas(gcf,[strtok(input_patten,'.') '_CDF.png']);
end

disp('Done!');