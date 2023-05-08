channel = "OverTheAir";
if strcmpi(channel, "OverTheAir")
    deviceName = "AD936x";
    channelNumber = 165;
    frequencyBand = 5;
    txGain = -10;
    rxGain = 10;
elseif strcmpi(channel, "GaussianNoise")
    % Specify the SNR of received signal for a simulated channel
    SNR = 10;
end

if strcmpi(channel,"OverTheAir")

    sdrReceiver = sdrrx(deviceName);
    sdrReceiver.BasebandSampleRate = sdrTransmitter.BasebandSampleRate;
    sdrReceiver.CenterFrequency = sdrTransmitter.CenterFrequency;
    sdrReceiver.OutputDataType = 'double';
    sdrReceiver.GainSource = 'Manual';
    sdrReceiver.Gain = rxGain;

    if ~strcmpi(deviceName,"Pluto")
        sdrReceiver.ShowAdvancedProperties = true;
        sdrReceiver.BypassUserLogic = true;
    end

    % Configure the capture length equivalent to twice the length of the
    % transmitted signal, this is to ensure that PSDUs are received in order.
    % On reception the duplicate MAC fragments are removed.
    sdrReceiver.SamplesPerFrame = 130000;
    fprintf('\nStarting a new RF capture.\n')
    % Set the capture duration and samples per frame
    captureDuration = 6; % in seconds
    samplesPerFrame = sdrReceiver.SamplesPerFrame;

    % Create a buffer to store the received waveform
    rxWaveformBuffer = [];

    % Capture and save the received waveform using a while loop
    startTime = tic;
    while toc(startTime) < captureDuration
        % Capture the received signal
        rxWaveform = capture(sdrReceiver, samplesPerFrame, 'Samples');

        % Concatenate the received waveform to the buffer
        rxWaveformBuffer = [rxWaveformBuffer; rxWaveform];
    end

    % Save the raw buffer waveform to a file
    save('rxWaveformBuffer.mat', 'rxWaveformBuffer');
% elseif strcmpi(channel,"GaussianNoise")
%     rxWaveform = awgn(txWaveform,SNR,'measured');
% else % No Impairments
%     rxWaveform = txWaveform;
end
%% if using SDR, release the objects
if strcmpi(channel,"OverTheAir")
    % release(sdrTransmitter);
    release(sdrReceiver);
end

% raised cosine filter parameters
% Design a custom raised cosine filter
rollOffFactor = 0.05;
span = 8; % Filter span in symbols
sps = 1; % Samples per symbol
numCoeffs = span*sps+1; % Number of filter coefficients
coeffs = rcosdesign(rollOffFactor, span, sps, 'sqrt'); % Filter coefficients

