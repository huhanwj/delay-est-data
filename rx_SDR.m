% Receiver

% Parameters
centerFrequency = 2.2e9;
bandwidth = 20e6;
gain = -10; % Adjust according to your requirements
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

%% Create SDR receiver object
rx = sdrrx('AD936x', ...
        'IPAddress', '192.168.3.2', ...
        'CenterFrequency', centerFrequency, ...
        'GainSource', 'Manual', ...
        'Gain', gain, ...
        'SamplesPerFrame', nSamples, ...
        'OutputDataType', 'double', ...
        'ChannelMapping', [1, 2]); % Enable both receive channels

% Receive the WiFi signal
endTime = datetime('now') + seconds(captureTime);
while datetime('now') < endTime
rxSignals = rx(); % rxSignals is a matrix where each column corresponds to an antenna
filteredSignal_save = cell(1,2);
% Process each received signal
for antennaIdx = 1:2
rxSignal = rxSignals(:, antennaIdx);

% Apply matched filtering
filteredSignal = rxFilter(rxSignal);
filteredSignal_save{antennaIdx} = filteredSignal;

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
% Save the filtered signals for both antennas
save('real_received_signal_antenna1.mat', 'filteredSignal_save{1}');
save('real_received_signal_antenna2.mat', 'filteredSignal_save{2}');

% Release the receiver
release(rx);