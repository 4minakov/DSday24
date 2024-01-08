clear all
close all

% Now let's listen to an Earth's son
y=load('earthquake_data.txt');
fs = 1;
dt=1/fs;
t = (0:length(y)-1)*dt;
[wt,f] = cwt(y,fs);

P = abs(wt);

figure,
subplot(2,1,1)
plot(t,y), axis tight
title('Kobe Earthquake and Scalogram'),xlabel('Time (s)'),ylabel('Acceleration (nm/s^2)')
subplot(2,1,2)
surface(t,f,P)
axis tight,shading flat,xlabel('Time (s)'),ylabel('Frequency (Hz)')
set(gca,'yscale','log'), 

%% extract and play frequencies
P1 = P/max(abs(P(:))); % loudness
% P1(P1<0.1) = NaN;
figure, imagesc(P1), colorbar
%%
[nfrq,nt] = size(P1);
% initialize sound array
y1 = zeros(nt,1);
%
fs_ac = 44.1e3;
dur = 0.5; % duration of note
dt_ac = 1/fs_ac;
nt_ac = dur/dt_ac;
t_ac = 0:dt_ac:dt_ac*(nt_ac-1);
sc_coeff = 1000;
yout = 0;

for i = 1:10:nt
    i
    p = P1(:,i);
    y = zeros(1,nt_ac);
    for j = 1:nfrq
            f_sc = f(j)*sc_coeff;
            % synthetize note and add
            y = y + p(j)*sin(pi*t_ac/dur).*sin(2*pi*t_ac*f_sc);
    end
    y = y/max(abs(y)) * max(abs(p));
    % add some zeros at the end to separate notes and append notes together
    %yout = [yout, [y, zeros(1, floor(0.3*dur*fs_ac))]];
    yout = [yout, y];
end

% sound(yout,fs_ac)
% save the song to a file
audiowrite('kobe_earthquake_2.wav', yout, fs_ac);
tout = [0:(length(yout)-1)]*dt_ac;
figure,
subplot(2,1,1)
plot(tout,yout), axis tight
title('Sonified  Earthquake and Scalogram'),xlabel('Time (s)'),ylabel('Amplitude')
subplot(2,1,2)
imagesc(tout,f*sc_coeff,P1)
axis tight,shading flat,xlabel('Time (s)'),ylabel('Frequency (Hz)')
set(gca,'yscale','log','Ydir','normal'),

