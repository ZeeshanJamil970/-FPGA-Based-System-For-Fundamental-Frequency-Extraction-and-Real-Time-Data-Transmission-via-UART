% List available serial ports
availablePorts = serialportlist("available");
disp('Available ports:');
disp(availablePorts);

% Define the serial port and baud rate
port = 'COM5'; % Adjust to your specific COM port
baudRate = 115200;

% Create and configure the serial port object
s = serial(port, 'BaudRate', baudRate, 'InputBufferSize', 1024);

% Open the serial port connection
fopen(s);

% Check if the serial port is open
if strcmp(s.Status, 'open')
    disp('Serial port is open');
else
    disp('Failed to open serial port');
end

% Number of 16-bit samples to read
numOfSamples = 1000; % Adjust based on your requirements
numOfBytesToRead = numOfSamples * 2; % Assuming 16-bit samples

% Read 8-bit binary data from the serial port
rawData = fread(s, numOfBytesToRead, 'uint8');

% Initialize an array to hold the 16-bit data
data = zeros(1, numOfSamples, 'uint16');

% Combine pairs of 8-bit values into 16-bit values
for i = 1:numOfSamples
    data(i) = bitor(bitshift(uint16(rawData(2*i-1)), 8), uint16(rawData(2*i)));
end

% Display the received 16-bit data
disp('Data received:');
disp(data);

% Close the serial port connection
fclose(s);
delete(s);
clear s;

% Define the sampling frequency
fs = 64000;  % 64 kHz

% Define the filter specifications
f0 = 8000;   % Center frequency (8 kHz)
bw = 2000;   % Bandwidth (2 kHz)
Fpass1 = 7000;  % Lower passband edge
Fpass2 = 9000;  % Upper passband edge
Fstop1 = 6000;  % Lower stopband edge
Fstop2 = 10000; % Upper stopband edge

% Normalize the frequencies by the Nyquist frequency
Fpass1 = Fpass1 / (fs / 2);
Fpass2 = Fpass2 / (fs / 2);
Fstop1 = Fstop1 / (fs / 2);
Fstop2 = Fstop2 / (fs / 2);

% Design the FIR bandpass filter using designfilt
bpFilt = designfilt('bandpassfir', ...
         'StopbandFrequency1', Fstop1, ...
         'PassbandFrequency1', Fpass1, ...
         'PassbandFrequency2', Fpass2, ...
         'StopbandFrequency2', Fstop2, ...
         'FilterOrder', 50, ...
         'SampleRate', fs);

% Visualize the filter response
fvtool(bpFilt);

% Export the filter coefficients for use in FPGA (optional)
b = bpFilt.Coefficients;

% Apply the filter to the received data
filteredData = filter(bpFilt, double(data));

% Convert the received data to voltage if needed
% Assuming 8-bit ADC with a reference voltage of 3.3V
% Adjust based on your ADC resolution and Vref
receivedVoltage = data * (3.3 / 255);

% Plot the received signal
t = (0:numOfSamples-1) / fs; % Assuming a 64 kHz sampling rate
figure;
subplot(2, 1, 1);
plot(t, receivedVoltage);
title('Received Signal');
xlabel('Time (s)');
ylabel('Amplitude (V)');

% Plot the filtered signal
subplot(2, 1, 2);
plot(t, filteredData);
title('Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude (V)');
