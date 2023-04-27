% Transmitter
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
rollOffFactor = 0.5;
txFilter = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', rollOffFactor, ...
    'FilterSpanInSymbols', 8, ...
    'OutputSamplesPerSymbol', 1);

% Apply pulse shaping
shapedSignal = txFilter(wifiSignal);

% Simulate a more realistic channel
KFactor = 4; % Rician K-factor
pathDelays = [0, 1.5e-8, 3e-8]; % Path delays
avgPathGains = [0, -3, -6]; % Average path gains in dB

channel = comm.RicianChannel(...
    'SampleRate', bandwidth, ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', avgPathGains, ...
    'KFactor', KFactor, ...
    'MaximumDopplerShift', 0, ...
    'DirectPathInitialPhase', 0, ...
    'DirectPathDopplerShift', 0, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', seed, ...
    'PathGainsOutputPort', true);

[directPathSignal, pathGains] = channel(shapedSignal);

% Receiver

% Raised Cosine Receive Filter
rxFilter = comm.RaisedCosineReceiveFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', rollOffFactor, ...
    'FilterSpanInSymbols', 8, ...
    'InputSamplesPerSymbol', 1, ...
    'DecimationFactor', 1);

% Apply matched filtering
filteredSignal = rxFilter(directPathSignal);

% Coarse frequency offset estimation and correction
coarseEst = wlanCoarseCFOEstimate(filteredSignal, cfg);
correctedSignal = wlanFrequencyCorrection(filteredSignal, coarseEst);

% Fine frequency offset estimation and correction
fineEst = wlanFineCFOEstimate(correctedSignal, cfg);
correctedSignal = wlanFrequencyCorrection(correctedSignal, fineEst);

% Symbol timing synchronization
[ind, symOffset] = wlanSymbolTimingEstimate(correctedSignal, cfg);
offsetCorrectedSignal = correctedSignal(ind:end);

% Channel estimation
demodSignal = wlanNonHTDataRecover(offsetCorrectedSignal, cfg);
[eqSig, chEst] = wlanEqualize(offsetCorrectedSignal, demodSignal, cfg);

% Demodulate and decode the received signal
rxPSDU = wlanNonHTDataRecover(eqSig, cfg);

% Decode the received MAC frame
[rxPayload, rxFrameType] = wlanMACFrame(rxPSDU, cfg);

% Save the received signal from the simulation
save('simulated_received_signal.mat', 'filteredSignal');
