% Receiver

% Parameters
centerFrequency = 2.2e9;
bandwidth = 20e6;
gain = 10; % Adjust according to your requirements
rollOffFactor = 0.5;
nSamples = 1024; % Adjust according to your WiFi signal parameters
captureTime = 10; % Time in seconds for capturing the signal

% Raised Cosine Receive Filter
rxFilter = comm.RaisedCosineReceiveFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', rollOffFactor, ...
    'FilterSpanInSymbols', 8, ... % Adjust according to your requirements
    'InputSamplesPerSymbol', 1, ... % Adjust according to your requirements
    'DecimationFactor', 1); % Adjust according to your requirements

% Create SDR receiver object
rx = sdrrx('ZC706 and FMCOMMS2/3/4', ...
           'CenterFrequency', centerFrequency, ...
           'GainSource', 'Manual', ...
           'Gain', gain, ...
           'SamplesPerFrame', nSamples, ...
           'OutputDataType', 'double');

% Create a non-HT (802.11g) format configuration object
cfg = wlanNonHTConfig('MCS', 0); % Default MCS (BPSK, rate 1/2)

% Receive the WiFi signal
endTime = datetime('now') + seconds(captureTime);
while datetime('now') < endTime
    rxSignal = rx();

    % Apply matched filtering
    filteredSignal = rxFilter(rxSignal);

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

    % Compare the received payload with the transmitted payload
    if isequal(rxPayload, payload)
        disp('Payload successfully received and decoded!');
    else
        disp('Payload not received correctly.');
    end
end
save('real_received_signal.mat', 'filteredSignal');

% Release the receiver
release(rx);
