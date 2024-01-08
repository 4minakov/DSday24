clear all
close all

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
figure, imagesc(P1), colorbar
%%
[nfrq,nt] = size(P1);
% initialize sound array
y1 = zeros(nt,1);
%
fs_ac = 44.1e3;
dur = 0.3; % duration of note
dt_ac = 1/fs_ac;
nt_ac = dur/dt_ac;
t_ac = 0:dt_ac:dt_ac*(nt_ac-1);
sc_coeff = 1000;
yout = 0;

type = 'saw' %choose synthesis type

for i = 1:5:nt
    i
    p = P1(:,i);
    y = 0;
    for j = 1:nfrq
        
        f_sc = f(j)*sc_coeff;
        
        if (strcmp(type,'sine'))
            
            y1 = p(j)*sin(pi*t_ac/dur).*sin(2*pi*t_ac*f_sc);
            
        elseif (strcmp(type,'saw'))
            
            T = (1/f_sc)*fs_ac;     % period in fractional samples
            ramp = (0:(nt_ac-1))/T;
            y1 = ramp-fix(ramp);
            y1 = p(j).*y1;
            y1 = y1 - mean(y1);
            
        elseif (strcmp(type,'fm'))
            envel = interp1([0 dur/6 dur/3 dur/5 dur], [0 1 .75 .6 0], t_ac);
            I_env = 5.*envel;
            y1 = envel.*sin(2.*pi.*f_sc.*t_ac + I_env.*sin(2.*pi.*f_sc.*t_ac));
            
        else
            error('Unknown synthesis type');
        end
        y = y + y1;
    end
    y = y/max(abs(y)) * max(abs(p));
    % add some zeros at the end to separate notes and append notes together
    yout = [yout, [y, zeros(1, floor(0.01*dur*fs_ac))]];
end

% sound(yout,fs_ac)
% save the song to a file
audiowrite('kobe_earthquake_2.wav', yout, fs_ac);
%%
tout = [0:(length(yout)-1)]*dt_ac;
figure,
subplot(2,1,1)
plot(tout,yout), axis tight
title('Sonified  Earthquake and Scalogram'),xlabel('Time (s)'),ylabel('Amplitude')
subplot(2,1,2)
imagesc(tout,f*sc_coeff,P1)
axis tight,shading flat,xlabel('Time (s)'),ylabel('Frequency (Hz)')
set(gca,'yscale','log','Ydir','normal'),


