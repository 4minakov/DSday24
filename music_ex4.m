clear all
close all


%piano notes in sequence starting from A4 (440 Hz)
notename = {'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'}; 

% duration of each note (seconds)
dur = 0.3; 

%sampling frequency (hertz)
fs = 8.192e3; 

% time vector
dt = 1/fs;
nt = dur/dt;

t = 0:dt:dt*(nt-1);

% let's play a simple song
song = {'A' 'A' 'E' 'E' 'F#' 'F#' 'E' 'E' 'D' 'D' 'C#' 'C#' 'B' 'B' 'A' 'A'};

yout = []; % initial the song file
for i=1:length(song)
    
    n = find(strcmp(song(i), notename));

    %convert sequential note number to frequency
    f = 440*2.^((n-1)/12);
    
    % harmonic oscillation + amplitude is damped at the begining and at the
    % end using a sine window to make smooth transitions
    y = sin(pi*t/dur).*sin(2*pi*t*f);
    % add some zeros at the end to separate notes and append notes together 
    yout = [yout, [y, zeros(1, floor(0.3*dur*fs))]];
    
    % plot a numerical representation of a musical note
    figure(1),clf,
    plot(t,y,'o-'), title(notename{n}), xlim([0.1 0.11])
    xlabel('time (s)'), ylabel('amplitude')
    drawnow
    % play note
    soundsc(y,fs)
    pause(dur)
end

% save the song to a file
audiowrite('little_song_smooth.wav', yout, fs);

tout = (0:length(yout)-1)*dt;
%% Time-frequency analysis using CWT
[wt,frq] = cwt(yout,fs);
figure,
subplot(2,1,1)
plot(tout,yout), axis tight
title('Tune and CWT'),ylabel('Amplitude')
subplot(2,1,2)
surface(tout,frq,abs(wt))
axis tight,shading flat,xlabel('Time (s)'),ylabel('Frequency (Hz)')
set(gca,'yscale','log'), 
ylim([200 2000])
