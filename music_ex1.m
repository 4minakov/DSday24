clear all
close all


%piano notes in sequence starting from A4 (440 Hz)
notename = {'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'}; 

% duration of each note (seconds)
dur = 0.3; 

%sampling frequency (hertz)
fs = 44.1e3; 

% time vector
dt = 1/fs;
nt = dur/dt;

t = 0:dt:dt*(nt-1);

for n=1:length(notename)

    %convert sequential note number to frequency
    f = 440*2.^((n-1)/12);
    
    % simple harmonic oscillation
    y = sin(2*pi*t*f);
    
    % plot a numerical representation of a musical note
    figure(1),clf,
    plot(t,y,'o-'), title(notename{n}), xlim([0.1 0.11])
    xlabel('time (s)'), ylabel('amplitude')
    drawnow
    % play note
    soundsc(y,fs)
    pause(dur)
end
%%
