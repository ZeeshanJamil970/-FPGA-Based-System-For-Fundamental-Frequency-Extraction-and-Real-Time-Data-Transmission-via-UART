% Parameters
bitsPerSample = 14;
sampleRate = 100000; % High sample rate to accurately capture the 5 kHz signal
signalFreq = 8000; % Frequency of the square wave
voltageRange = 3; % 0 to 3V
duration = 0.01; % Duration of the signal in seconds
numSamples = sampleRate * duration; % Number of samples to generate

% Generate a synthetic square wave (0 to 3V) at 8 kHz
t = (0:numSamples-1) / sampleRate;
originalSignal = 1.5 * (square(2 * pi * signalFreq * t) + 1); % 0 to 3V square wave

% Step 1: ADC - Convert voltage signal to 14-bit samples
maxVal = 2^bitsPerSample - 1;
data14Bit = round((originalSignal / voltageRange) * maxVal);

% Step 2: Serialize the data with start and stop bits
startBit = 1;
stopBit = 0;
frameSize = bitsPerSample + 2; % 14 bits data + 1 start bit + 1 stop bit
numFrames = length(data14Bit);
dataStream = zeros(1, numFrames * frameSize);

for i = 1:numFrames
    sample = data14Bit(i);
    % Convert sample to binary string
    binSample = dec2bin(sample, bitsPerSample) - '0';
    % Concatenate start bit, sample, and stop bit
    frame = [startBit, binSample, stopBit];
    % Insert the frame into the data stream
    startIdx = (i - 1) * frameSize + 1;
    endIdx = startIdx + frameSize - 1;
    dataStream(startIdx:endIdx) = frame;
end

% Step 3: Deserialize the data and regenerate the signal
receivedData = dataStream;

% Decode received data (strip start and stop bits)
decodedData = zeros(1, numFrames);

for i = 1:numFrames
    startIdx = (i - 1) * frameSize + 1;
    endIdx = startIdx + frameSize - 1;
    frame = receivedData(startIdx:endIdx);
    
    if frame(1) == startBit && frame(end) == stopBit
        % Extract the 14-bit sample
        sampleBits = frame(2:end-1);
        decodedData(i) = bin2dec(char(sampleBits + '0'));
    end
end

% Step 4: DAC - Convert 14-bit samples back to voltage
regeneratedSignal = (decodedData / maxVal) * voltageRange;

% Step 5: Design and apply Kaiser window FIR filter to filter out first harmonic
Fs = sampleRate;  % Sampling Frequency
N = 60;           % Order
Fc = 8000;        % Cutoff Frequency
Beta = 0.5;       % Window Parameter
flag = 'scale';   % Sampling Flag

% Create the window vector for the design algorithm
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function
b = fir1(N, Fc/(Fs/2), 'low', win, flag);

% Filter the regenerated signal to extract the first harmonic
filteredSignal = filter(b, 1, regeneratedSignal);

% Step 6: Compute FFT of the original, regenerated, and filtered signals
n = length(originalSignal);
f = (0:n-1)*(sampleRate/n); % Frequency range

% Compute FFT
fftOriginal = fft(originalSignal);
fftRegenerated = fft(regeneratedSignal);
fftFiltered = fft(filteredSignal);

% Plot the original, regenerated, and filtered signals
figure;
subplot(3,1,1);
plot(t, originalSignal, 'b-', 'DisplayName', 'Original Signal');
title('Original Signal');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend show;
axis([0 duration 0 voltageRange]); % Adjust axis to show the signal clearly

subplot(3,1,2);
plot(t, regeneratedSignal, 'r--', 'DisplayName', 'Regenerated Signal');
title('Regenerated Signal');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend show;
axis([0 duration 0 voltageRange]); % Adjust axis to show the signal clearly

subplot(3,1,3);
plot(t, filteredSignal, 'g-', 'DisplayName', 'Filtered Signal (First Harmonic)');
title('Filtered Signal (First Harmonic)');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend show;
axis([0 duration 0 voltageRange]); % Adjust axis to show the signal clearly

% Plot FFT of the original, regenerated, and filtered signals
figure;
subplot(3,1,1);
plot(f, abs(fftOriginal), 'b-', 'DisplayName', 'Original Signal FFT');
title('FFT of Original Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend show;
axis([0 sampleRate/2 0 max(abs(fftOriginal))]); % Show up to Nyquist frequency

subplot(3,1,2);
plot(f, abs(fftRegenerated), 'r--', 'DisplayName', 'Regenerated Signal FFT');
title('FFT of Regenerated Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend show;
axis([0 sampleRate/2 0 max(abs(fftRegenerated))]); % Show up to Nyquist frequency

subplot(3,1,3);
plot(f, abs(fftFiltered), 'g-', 'DisplayName', 'Filtered Signal FFT');
title('FFT of Filtered Signal (First Harmonic)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend show;
axis([0 sampleRate/2 0 max(abs(fftFiltered))]); % Show up to Nyquist frequency