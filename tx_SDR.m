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
%% Prepare data for transmission
fData = zeros(100,100,'uint8'); % should equal to a all black grayscale image
fData_t = fData(:); % Convert to vector
%% Fragment Transmit Data
msduLength = 2304; % Number of bytes in PSDU
numMSDUs = ceil(length(fData_t)/msduLength); % Number of required MPDUs
padZeros = msduLength - mod(length(fData_t),msduLength); % Number of zeros to pad
txData = [fData_t;zeros(padZeros,1)]; % Pad PSDU with zeros
txDataBits = double(int2bit(txData,8,false)); % Convert PSDU to bits
% Divide input data stream into fragments
bitsPerOctet = 8;
data = zeros(0,1);
for i=0:numMSDUs-1

    % Extract image data (in octets) for each MPDU
    frameBody = txData(i*msduLength+1:msduLength*(i+1),:);

    % Create MAC frame configuration object and configure sequence number
    cfgMAC = wlanMACFrameConfig(FrameType='Data',SequenceNumber=i);

    % Generate MPDU
    [psdu, lengthMPDU]= wlanMACFrame(frameBody,cfgMAC,OutputFormat='bits');

    % Concatenate PSDUs for waveform generation
    data = [data; psdu]; %#ok<AGROW>
end

%% Generate 802.11a Baseband WLAN Signal
nonHTcfg = wlanNonHTConfig;       % Create packet configuration
nonHTcfg.MCS = 6;                 % Modulation: 64QAM Rate: 2/3
nonHTcfg.NumTransmitAntennas = 1; % Number of transmit antenna
chanBW = nonHTcfg.ChannelBandwidth;
nonHTcfg.PSDULength = lengthMPDU; % Set the PSDU length
% Initialize the scrambler with a random integer for each packet
scramblerInitialization = randi([1 127],numMSDUs,1);

osf = 1.5;

sampleRate = wlanSampleRate(nonHTcfg); % Nominal sample rate in Hz
% Define idle time
idleTime = 20e-6; % 20 microseconds
% Calculate idle time samples based on the nominal sample rate
idleTimeSamples = round(idleTime * sampleRate * osf);

% % Generate baseband NonHT packets separated by idle time
% txWaveform = wlanWaveformGenerator(data,nonHTcfg, ...
%     NumPackets=numMSDUs,IdleTime=20e-6, ...
%     ScramblerInitialization=scramblerInitialization,...
%     OversamplingFactor=osf);

% Generate L-STF, L-LTF, and L-SIG fields
lstf = wlanLSTF(nonHTcfg, overSamplingFactor=osf);
lltf = wlanLLTF(nonHTcfg, overSamplingFactor=osf);
lsig = wlanLSIG(nonHTcfg, overSamplingFactor=osf);

% Generate the non-HT Data field for each packet and concatenate with idle time
waveform = [lstf; lltf; lsig];
for i = 1:numMSDUs
    psdu = data((i-1)*lengthMPDU*8+1:lengthMPDU*8*i);
    nonhtdata = wlanNonHTData(psdu, nonHTcfg, scramblerInitialization(i), OversamplingFactor = osf);
    
    % Concatenate the fields to create the complete waveform
    waveform = [waveform; nonhtdata];
    
    % Add idle time if not the last packet
    if i ~= numMSDUs
        waveform = [waveform; zeros(idleTimeSamples, 1)];
    end
end


% Design a custom raised cosine filter
rollOffFactor = 0.05;
span = 8; % Filter span in symbols
sps = 1; % Samples per symbol
numCoeffs = span*sps+1; % Number of filter coefficients
coeffs = rcosdesign(rollOffFactor, span, sps, 'sqrt'); % Filter coefficients

% Apply the filter to the waveform (pulse shaping)
shapedSignal = upfirdn(waveform, coeffs, sps);
if strcmpi(channel,"OverTheAir")

    % Transmitter properties
    sdrTransmitter = sdrtx(deviceName);
    sdrTransmitter.BasebandSampleRate = sampleRate*osf;
    sdrTransmitter.CenterFrequency = wlanChannelFrequency(channelNumber,frequencyBand);
    sdrTransmitter.Gain = txGain;
    
    % Pass the SDR I/O directly to host skipping FPGA on Zynq Radio or USRP
    % Embedded Series Radio
    if ~strcmpi(deviceName,"Pluto")
        sdrTransmitter.ShowAdvancedProperties = true;
        sdrTransmitter.BypassUserLogic = true;
    end

    fprintf('\nGenerating WLAN transmit waveform:\n')

    % Scale the normalized signal to avoid saturation of RF stages
    powerScaleFactor = 0.8;
    txWaveform = txWaveform.*(1/max(abs(txWaveform))*powerScaleFactor);

    % Transmit RF waveform
    transmitRepeat(sdrTransmitter,txWaveform);
end
%% if using SDR, release the objects
if strcmpi(channel,"OverTheAir")
    release(sdrTransmitter);
    % release(sdrReceiver);
end