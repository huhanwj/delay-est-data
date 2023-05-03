channel = "OverTheAir";
if strcmpi(channel, "OverTheAir")
    deviceName = "AD936x";
    channelNumber = 165;
    frequencyBand = 5;
    txGain = -10;
    RxGain = 10;
elseif strcmpi(channel, "GaussianNoise")
    % Specify the SNR of received signal for a simulated channel
    SNR = 10;
end
% Configure all the scopes and figures for the example
% Setup handle for image plot
if ~exist('imFig','var') || ~ishandle(imFig) %#ok<SUSENS> 
    imFig = figure;
    imFig.NumberTitle = 'off';
    imFig.Name = 'Image Plot';
    imFig.Visible = 'off';
else
    clf(imFig); % Clear figure
    imFig.Visible = 'off';
end

% Setup Spectrum viewer
spectrumScope = spectrumAnalyzer( ...
    SpectrumType='power-density', ...
    Title='Received Baseband WLAN Signal Spectrum', ...
    YLabel='Power spectral density', ...
    Position=[69 376 800 450]);

% Setup the constellation diagram viewer for equalized WLAN symbols
refQAM = wlanReferenceSymbols('64QAM');
constellation = comm.ConstellationDiagram(...
    Title='Equalized WLAN Symbols',...
    ShowReferenceConstellation=true,...
    ReferenceConstellation=refQAM,...
    Position=[878 376 460 460]);
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
    sdrReceiver.SamplesPerFrame = 30000;
    fprintf('\nStarting a new RF capture.\n')
    % Set the capture duration and samples per frame
    captureDuration = 10; % 10 seconds
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