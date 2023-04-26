% Transmitter

% Parameters
centerFrequency = 2.2e9;
bandwidth = 20e6;
gain = 10; % Adjust according to your requirements
rollOffFactor = 0.5;
% Set the seed for the random number generator
seed = 12345;
rng(seed);
% Generate a random payload
payload = randi([0, 1], 1000, 1);

% Create a non-HT (802.11g) format configuration object
cfg = wlanNonHTConfig('MCS', 0); % Default MCS (BPSK, rate 1/2)

% Generate the 802.11g frame
txPSDU = wlanMACFrame(payload, 'Data', cfg);
wifiSignal = wlanWaveformGenerator(txPSDU, cfg);

% Raised Cosine Transmit Filter
txFilter = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', rollOffFactor, ...
    'FilterSpanInSymbols', 8, ... % Adjust according to your requirements
    'OutputSamplesPerSymbol', 1); % Adjust according to your requirements

% Apply pulse shaping
shapedSignal = txFilter(wifiSignal);

% Create SDR transmitter object
rx = sdrrx('ZC706 and FMCOMMS2/3/4', ...
           'CenterFrequency', centerFrequency, ...
           'BasebandSampleRate', bandwidth, ...
           'Gain', gain);

% Transmit the shaped 802.11g signal
transmitRepeat(tx, shapedSignal);

% Release the transmitter
release(tx);
