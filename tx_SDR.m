% Transmitter

% Parameters
centerFrequency = 2.2e9;
bandwidth = 20e6;
gain = -10; % Adjust according to your requirements
rollOffFactor = 0.5;
% Set the seed for the random number generator
seed = 12345;
rng(seed);
% Generate a random payload
payload = randi([0, 1], 1000, 1);

% Create a non-HT (802.11g) format configuration object
cfg = wlanNonHTConfig('MCS', 0); % Default MCS (BPSK, rate 1/2)

% Generate the 802.11g frame
macConfig = wlanMACFrameConfig(FrameType="Data");
% macConfig.Address1 = [255,255,255,255,255,255]; % Broadcast
Address2 = dec2hex(randi([0, 255],1,6))';
Address3 = dec2hex(randi([0,255],1,6))';
macConfig.Address2 = Address2(:)';
macConfig.Address3 = Address3(:)';

[txPSDU,PSDUlen] = wlanMACFrame(payload, macConfig,"OutputFormat","bits");
% txPSDU = reshape(hex2dec(reshape(txPSDU,2,[])),1,[]);
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
tx = sdrtx('AD936x', ...
           'IPAddress', '192.168.3.2', ...
           'CenterFrequency', centerFrequency, ...
           'BasebandSampleRate', bandwidth, ...
           'Gain', gain, ...
           'ChannelMapping', [1, 2]); % Enable both transmit channels

% Transmit the shaped 802.11g signal
% Duplicate the shaped signal for both antennas
shapedSignal2Antennas = repmat(shapedSignal, 1, 2);
transmitRepeat(tx, shapedSignal2Antennas);

% Release the transmitter
release(tx);