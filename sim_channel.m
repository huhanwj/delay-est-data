% Transmitter
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

[txPSDU,PSDUlen] = wlanMACFrame(payload, macConfig,"OutputFormat","bits"); % Generate the MAC frame
wifiSignal = wlanWaveformGenerator(txPSDU, cfg);% Generate the 802.11g waveform
% Raised Cosine Transmit Filter
rollOffFactor = 0.5;
txFilter = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', rollOffFactor, ...
    'FilterSpanInSymbols', 8, ...
    'OutputSamplesPerSymbol', 1);

% Apply pulse shaping
shapedSignal = txFilter(wifiSignal);

% Generate shaped signal for both antennas
shapedSignal1 = shapedSignal;
shapedSignal2 = shapedSignal;

% Simulate a more realistic channel
% Assume perfect channel knowledge at the receiver
KFactor = 4; % Rician K-factor
pathDelays = [0, 1.5e-8, 3e-8]; % Path delays
avgPathGains = [0, -3, -6]; % Average path gains in dB

channel1 = comm.RicianChannel(...
    'SampleRate', 20e6, ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', avgPathGains, ...
    'KFactor', KFactor, ...
    'MaximumDopplerShift', 0, ...
    'DirectPathInitialPhase', 0, ...
    'DirectPathDopplerShift', 0, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', seed, ...
    'PathGainsOutputPort', true);

channel2 = clone(channel1); % Create an identical channel for the second transmit antenna

[directPathSignal1, pathGains1] = channel1(shapedSignal1);
[directPathSignal2, pathGains2] = channel2(shapedSignal2);

% Combine the received signals from both transmit antennas
receivedSignal = directPathSignal1 + directPathSignal2;

% Receiver

% Raised Cosine Receive Filter
rxFilter = comm.RaisedCosineReceiveFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', rollOffFactor, ...
    'FilterSpanInSymbols', 8, ...
    'InputSamplesPerSymbol', 1, ...
    'DecimationFactor', 1);

% Apply matched filtering
filteredSignal = rxFilter(receivedSignal);
% Apply matched filtering just for saving
filteredSignal1 = rxFilter(directPathSignal1);
filteredSignal2 = rxFilter(directPathSignal2);

pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port'); % Create a phase frequency offset object
% Coarse frequency offset estimation and correction
coarseEst = wlanCoarseCFOEstimate(filteredSignal, 'CBW20');
correctedSignal = pfOffset(filteredSignal, -coarseEst);

% Fine frequency offset estimation and correction
fineEst = wlanFineCFOEstimate(correctedSignal, 'CBW20');
correctedSignal = pfOffset(correctedSignal, -fineEst);

% Symbol timing synchronization
[ind, symOffset] = wlanSymbolTimingEstimate(correctedSignal, 'CBW20');
offsetCorrectedSignal = correctedSignal(ind:end);

% Channel estimation
% Using perfect channel knowledge for simplicity
demodSignal = wlanNonHTDataRecover(offsetCorrectedSignal, cfg);
[eqSig, chEst] = wlanEqualize(offsetCorrectedSignal, demodSignal, cfg);

% Combine the channel estimates for both antennas
combinedChEst = pathGains1 + pathGains2;

% Apply the MIMO channel estimates to the equalized signal
eqSigWithMIMO = eqSig .* conj(combinedChEst);

% Demodulate and decode the received signal
rxPSDU = wlanNonHTDataRecover(eqSigWithMIMO, cfg);

% Decode the received MAC frame
[rxPayload, rxFrameType] = wlanMACFrame(rxPSDU, cfg);

% Save the received signals from both antennas
save('simulated_received_signal_antenna1.mat', 'filteredSignal1');
save('simulated_received_signal_antenna2.mat', 'filteredSignal2');
