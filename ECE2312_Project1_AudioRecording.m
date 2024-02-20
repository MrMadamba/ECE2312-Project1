%  ECE 2312 - Project 1
%  Matthew Madamba

%********************************************
%           Record your Speech 
% *******************************************
%      Getting the device info IDs

% deviceInfo = audiodevinfo;
% 
% disp('Available audio input devices:');       %Display the different inputs and their ID's
% for i = 1:numel(deviceInfo.input)
%     disp(['Input Device ', num2str(i), ':']);
%     disp(deviceInfo.input(i));
%     disp(' ');
% end


%Creating recorder object
Fs = 48000;             % Sampling frequency (Hz)
BitsPerSample = 16;     % Bits per sample
NumberOfChannels = 1;   % Number of audio channels (mono)

% Input ID for each phrase (using 0 as built in mic)
InputDeviceIDs = [0, 0, 0]; 

% 5 seconds for each phrase
recordDuration = 5;  

% Store the recorded audio data for each phrase
recordedData = cell(1, numel(InputDeviceIDs));

% Create an audiorecorder object for the current phrase 
MyAudio = audiorecorder(Fs, BitsPerSample, NumberOfChannels, 0);

% Loop through each phrase
for i = 1:numel(InputDeviceIDs)      % For loop iterates 3 for each plot
    if i == 1                        % Display correct title each recording
        PhraseTitle = 'Quick Brown Fox';
    elseif i == 2
        PhraseTitle = 'We Promptly Judged';
    elseif i == 3
        PhraseTitle = 'Crazy Fredrick';
    end    
    disp(['Recording Phrase ', PhraseTitle]);
    % Record audio
    record(MyAudio);
    
    pause(recordDuration);  % Pause for the duration (5s)
    stop(MyAudio);

    % Store the recorded audio data
    recordedData{i} = getaudiodata(MyAudio);
end

% Plot the recorded waveforms for each phrase
figure;
for i = 1:numel(recordedData)   % For loop iterates 3 for each *subplot*
    if i == 1
        PhraseTitle = 'Quick Brown Fox';
    elseif i == 2
        PhraseTitle = 'We Promptly Judged';
    elseif i == 3
        PhraseTitle = 'Crazy Fredrick';
    end
    subplot(numel(recordedData), 1, i);
    timeInSeconds = (0:length(recordedData{i})-1) / Fs;
    plot(timeInSeconds, recordedData{i});
    title(['Recorded Audio - ', PhraseTitle])
    xlabel('Time (seconds)');
    ylabel('Amplitude');
end


%*****************************************************************
%                Visualization via Spectrogram
%*****************************************************************

%  Define parameters for spectrogram

windowLength = 1024;                 % Length of the window for computing the spectrogram 
overlap = round(windowLength * 0.9); % Overlap 
nfft = 2048;                         % Number of FFT points
frequencyRange = [0 8000];           % Given 8000Hz limit
figure;
for i = 1:numel(recordedData) % For loop iterates 3 for each plot 
    if i == 1
        PhraseTitle = 'Quick Brown Fox';
    elseif i == 2
        PhraseTitle = 'We Promptly Judged';
    elseif i == 3
        PhraseTitle = 'Crazy Fredrick';
    end

    % Generate spectrogram
    [s, f, t] = spectrogram(recordedData{i}, windowLength, overlap, nfft, Fs, 'yaxis');

    % Plot the spectrogram
    subplot(numel(recordedData), 1, i);
    imagesc(t, f, 10*log10(abs(s)));
    title(['Spectrogram of Recorded Audio -  ', PhraseTitle]);
    xlabel('Time (seconds)');
    ylabel('Frequency (Hz)');
    colorbar;
    axis xy;
    clim([-1 20]);
    % Limit the frequency range
    ylim(frequencyRange);
end


%*****************************************************************
%                Saving and Loading WAV Files
%*****************************************************************

% Define file names for the phrases
fileNames = {'the_quick_brown_fox.wav', 'we_promptly_judged.wav', 'crazy_fredrick.wav'};

% Save recorded speech signals as WAV files
for i = 1:numel(recordedData)
    fileName = fileNames{i};
    audiowrite(fileName, recordedData{i}, Fs);
end

% Input WAV files in a for loop to generate spectrograms
figure;
for i = 1:numel(fileNames)
    % Load the WAV file
    [wavData{i}, Fs] = audioread(fileNames{i});

    % Generate spectrogram
    [s, f, t] = spectrogram(wavData{i}, windowLength, overlap, nfft, Fs, 'yaxis');

    % Plot the spectrogram in subplots
    subplot(numel(fileNames), 1, i);
    imagesc(t, f, 10*log10(abs(s)));
    title(['Spectrogram of WAV File - ', strrep(fileNames{i}, '_', ' ')]); 
    xlabel('Time (seconds)');
    ylabel('Frequency (Hz)');
    colorbar;
    axis xy;
    clim([-1 20]);
    ylim(frequencyRange);
end


%*****************************************************************
%               Fun with stereo speech files
%*****************************************************************

%--------------------- Add Sample Delay --------------------------

% Load a previously recorded speech signal (fox)
[y, Fs] = audioread('the_quick_brown_fox.wav');

% Duplicate the column data
OriginalStereo = [y, y];

% Save the duplicated double column stereo signal as a WAV file
Original_0ms_wav = 'teamMatthewMadamba-stereosoundfile-0ms.wav'; % Original non-delayed signal 
audiowrite(Original_0ms_wav, OriginalStereo, Fs);


% Implement delays

% Ear Delay
% Delay the second column by inserting 26 zeros (calculated # of samples) at the beginning
Delayed_avghead_Stereo = [OriginalStereo(:,1), [zeros(26, 1); OriginalStereo(1:end-26, 2)]]; % Second columns get 26 zero's at the start and last 26 samples removed 
% Save the delayed signal as a WAV file
Delayed_avghead_wav = 'teamMatthewMadamba-stereosoundfile-avghead.wav'; 
audiowrite(Delayed_avghead_wav, Delayed_avghead_Stereo, Fs);


% 1ms, 10ms, and 100ms delays
delaySamples_1ms = round(1e-3 * Fs);  % Convert 1ms to samples = 48 samples 
delaySamples_10ms = round(10e-3 * Fs);  % Convert 10ms to samples = 480 samples
delaySamples_100ms = round(100e-3 * Fs);  % Convert 100ms to samples = 4800 samples

% Apply delays to the stereo signal - Only shifting 2nd column
% by the corresponding # of samples
Delayed_1ms_Stereo = [OriginalStereo(:,1), [zeros(delaySamples_1ms, 1); OriginalStereo(1:end-delaySamples_1ms, 2)]]; % Delay 2nd Column by 48 samples
Delayed_10ms_Stereo = [OriginalStereo(:,1), [zeros(delaySamples_10ms, 1); OriginalStereo(1:end-delaySamples_10ms, 2)]]; % Delay 2nd Column by 480 samples
Delayed_100ms_Stereo = [OriginalStereo(:,1), [zeros(delaySamples_100ms, 1); OriginalStereo(1:end-delaySamples_100ms, 2)]]; % Delay 2nd Column by 4800 samples

% Save the delayed signals as WAV files
Delayed_1ms_wav = 'teamMatthewMadamba-stereosoundfile-1msDelay.wav';
Delayed_10ms_wav = 'teamMatthewMadamba-stereosoundfile-10msDelay.wav';
Delayed_100ms_wav = 'teamMatthewMadamba-stereosoundfile-100msDelay.wav';
audiowrite(Delayed_1ms_wav, Delayed_1ms_Stereo, Fs);
audiowrite(Delayed_10ms_wav, Delayed_10ms_Stereo, Fs);
audiowrite(Delayed_100ms_wav, Delayed_100ms_Stereo, Fs);



% ----------------------- Attenuation ---------------------------

% Get 0ms delay file information columns 
[y0ms, Fs] = audioread('teamMatthewMadamba-stereosoundfile-0ms.wav');

% Get average head delay file information columns
[yAvgHead, ~] = audioread('teamMatthewMadamba-stereosoundfile-EarDelay.wav');

% -1.5dB, -3dB, and -6dB
attenuations = [-1.5, -3, -6]; 

% For loop to make decibel drop signals 
for i = 1:numel(attenuations)
    attenuationFactor = 10^(attenuations(i) / 20);      % Calculate the attenuation factor with equation
    
    yAttenuated = y0ms;                                 % Make new channel with identical columns to 0ms delay
    yAttenuated(:, 2) = y0ms(:, 2) * attenuationFactor; % Multiply the second column of the new channel by factor
    
    % Save the attenuated signal as a WAV file
    filename = sprintf('teamMatthewMadamba-stereosoundfile-0ms-%ddB.wav', abs(attenuations(i))); % Absolute value just to remove the negative sign 
    audiowrite(filename, yAttenuated, Fs);
end


% AvgHead delay for loop
for i = 1:numel(attenuations)
    attenuationFactor = 10^(attenuations(i) / 20);      % Calculate the attenuation factor with equation

    yAttenuated = yAvgHead;                             % Make new channel with identical columns to avghead
    yAttenuated(:, 2) = yAvgHead(:, 2) * attenuationFactor;  % Apply the attenuation to the second channel
    
    % Save the attenuated signal as a WAV file
    filename = sprintf('teamMatthewMadamba-avghead-%ddB.wav', abs(attenuations(i)));
    audiowrite(filename, yAttenuated, Fs);
end



