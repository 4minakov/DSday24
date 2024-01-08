clear all
close all

% Now let's listen to an Earth's son
y=load('earthquake_data.txt');
fs = 1;
dt=1/fs;
t = (0:length(y)-1)*dt;
[wt,f] = cwt(y,fs);

figure,
subplot(2,1,1)
plot(t,y), axis tight
title('Kobe Earthquake and Scalogram'),xlabel('Time (s)'),ylabel('Acceleration (nm/s^2)')
subplot(2,1,2)
surface(t,f,abs(wt))
axis tight,shading flat,xlabel('Time (s)'),ylabel('Frequency (Hz)')
set(gca,'yscale','log'), 

%% 
% let's map data to audio domain
fs_ac = 44.1e3;
t_scaled = t/1000; %about 3 sec
[y_ac, t_ac] = resample(y,t_scaled,fs_ac);
y_ac = y_ac/max(abs(y_ac));
figure, plot(t_ac,y_ac,'-'), title('signal scaled in acoustic range')
xlabel('Time (s)'), ylabel('Acceleration (nm/s^2)'), axis tight


[wt_ac,f_ac] = cwt(y_ac,fs_ac);

sound(y_ac, fs_ac)
% save the scaled earthquake signal to a file
audiowrite('kobe_earthquake.wav', y_ac, fs_ac);

figure,
subplot(2,1,1)
plot(t_ac,y_ac), axis tight
title('Kobe Earthquake and Scalogram'),xlabel('Scaled Time (s)'),ylabel('Acceleration (nm/s^2)')
subplot(2,1,2)
surface(t_ac,f_ac,abs(wt_ac))
axis tight,shading flat,xlabel('Scaled time (s)'),ylabel('Scaled frequency (Hz)')
set(gca,'yscale','log'), ylim([1 1000])
