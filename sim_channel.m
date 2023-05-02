% % Transmitter
% % Set the seed for the random number generator
% seed = 12345;
% rng(seed);
% % Generate a random payload
% payload = randi([0, 1], 1000, 1);

% % Create a non-HT (802.11g) format configuration object
% cfg = wlanNonHTConfig('MCS', 0); % Default MCS (BPSK, rate 1/2)

% % Generate the 802.11g frame
% macConfig = wlanMACFrameConfig(FrameType="Data");
% % macConfig.Address1 = [255,255,255,255,255,255]; % Broadcast
% Address2 = dec2hex(randi([0, 255],1,6))';
% Address3 = dec2hex(randi([0,255],1,6))';
% macConfig.Address2 = Address2(:)';
% macConfig.Address3 = Address3(:)';

% [txPSDU,PSDUlen] = wlanMACFrame(payload, macConfig,"OutputFormat","bits"); % Generate the MAC frame
% txLLTF = wlanLLTF(cfg); % Generate the non-HT long preamble
% wifiSignal = wlanWaveformGenerator(txPSDU, cfg);% Generate the 802.11g waveform
% % Raised Cosine Transmit Filter
% rollOffFactor = 0.5;
% txFilter = comm.RaisedCosineTransmitFilter(...
%     'Shape', 'Square root', ...
%     'RolloffFactor', rollOffFactor, ...
%     'FilterSpanInSymbols', 8, ...
%     'OutputSamplesPerSymbol', 1);

% % Apply pulse shaping
% shapedSignal = txFilter(wifiSignal);

% % Generate shaped signal for both antennas
% shapedSignal1 = shapedSignal;
% shapedSignal2 = shapedSignal;

% % Simulate a more realistic channel
% % Assume perfect channel knowledge at the receiver
% KFactor = 4; % Rician K-factor
% pathDelays = [0, 1.5e-8, 3e-8]; % Path delays
% avgPathGains = [0, -3, -6]; % Average path gains in dB

% channel1 = comm.RicianChannel(...
%     'SampleRate', 20e6, ...
%     'PathDelays', pathDelays, ...
%     'AveragePathGains', avgPathGains, ...
%     'KFactor', KFactor, ...
%     'MaximumDopplerShift', 0, ...
%     'DirectPathInitialPhase', 0, ...
%     'DirectPathDopplerShift', 0, ...
%     'RandomStream', 'mt19937ar with seed', ...
%     'Seed', seed, ...
%     'PathGainsOutputPort', true);

% channel2 = clone(channel1); % Create an identical channel for the second transmit antenna

% [directPathSignal1, pathGains1] = awgn(channel1(shapedSignal1), 10);
% [directPathSignal2, pathGains2] = awgn(channel2(shapedSignal2), 10);
% rxLLTF = awgn(channel1(txLLTF), 10); % pass the LLTF through the channel
% dLLTF = wlanLLTFDemodulate(rxLLTF,cfg);
% chEst = wlanLLTFChannelEstimate(dLLTF,cfg);
% % Combine the received signals from both transmit antennas
% receivedSignal = directPathSignal1 + directPathSignal2;

% % Receiver

% % Raised Cosine Receive Filter
% rxFilter = comm.RaisedCosineReceiveFilter(...
%     'Shape', 'Square root', ...
%     'RolloffFactor', rollOffFactor, ...
%     'FilterSpanInSymbols', 8, ...
%     'InputSamplesPerSymbol', 1, ...
%     'DecimationFactor', 1);

% % Apply matched filtering
% filteredSignal = rxFilter(receivedSignal);
% % Apply matched filtering just for saving
% filteredSignal1 = rxFilter(directPathSignal1);
% filteredSignal2 = rxFilter(directPathSignal2);

% pfOffset = comm.PhaseFrequencyOffset('SampleRate',20e6,'FrequencyOffsetSource','Input port'); % Create a phase frequency offset object
% % Coarse frequency offset estimation and correction
% coarseEst = wlanCoarseCFOEstimate(filteredSignal, 'CBW20');
% correctedSignal = pfOffset(filteredSignal, -coarseEst);

% % Fine frequency offset estimation and correction
% fineEst = wlanFineCFOEstimate(correctedSignal, 'CBW20');
% correctedSignal = pfOffset(correctedSignal, -fineEst);

% % Symbol timing synchronization
% [ind, symOffset] = wlanSymbolTimingEstimate(correctedSignal, 'CBW20');
% offsetCorrectedSignal = correctedSignal(ind:end);

% % Channel estimation
% % Using perfect channel knowledge for simplicity
% demodSignal = wlanNonHTDataRecover(offsetCorrectedSignal, chEst,0.1,cfg);
% % [eqSig, chEst] = wlanEqualize(offsetCorrectedSignal, demodSignal, cfg);

% % Combine the channel estimates for both antennas
% combinedChEst = pathGains1 + pathGains2;

% % Apply the MIMO channel estimates to the equalized signal
% eqSigWithMIMO = eqSig .* conj(combinedChEst);

% % Demodulate and decode the received signal
% rxPSDU = wlanNonHTDataRecover(eqSigWithMIMO, cfg);

% % Decode the received MAC frame
% [rxPayload, rxFrameType] = wlanMACFrame(rxPSDU, cfg);

% % Save the received signals from both antennas
% save('simulated_received_signal_antenna1.mat', 'filteredSignal1');
% save('simulated_received_signal_antenna2.mat', 'filteredSignal2');
channel = "GaussianNoise";
if strcmpi(channel, "OverTheAir")
    deviceName = "AD936x";
    channelNumber = 5;
    frequencyBand = 2.4;
    txGain = -10;
    RxGain = 10;
elseif strcmpi(channel, "GaussianNoise")
    % Specify the SNR of received signal for a simulated channel
    SNR = 20;
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
%% Fragment Transmit Data
msduLength = 2304; % Number of bytes in PSDU
numMSDUs = ceil(length(fData)/msduLength); % Number of required MPDUs
padZeros = msduLength - mod(length(fData),msduLength); % Number of zeros to pad
txData = [fData;zeros(padZeros,1)]; % Pad PSDU with zeros
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

% % Generate baseband NonHT packets separated by idle time
% txWaveform = wlanWaveformGenerator(data,nonHTcfg, ...
%     NumPackets=numMSDUs,IdleTime=20e-6, ...
%     ScramblerInitialization=scramblerInitialization,...
%     OversamplingFactor=osf);

% Generate L-STF, L-LTF, and L-SIG fields
lstf = wlanLSTF(nonHTcfg);
lltf = wlanLLTF(nonHTcfg);
lsig = wlanLSIG(nonHTcfg);

% Generate the non-HT Data field without guard interval
nonhtdata_nogi = wlanNonHTData(data, nonHTcfg, 'ScramblerInitialization', scramblerInitialization, 'GuardInterval', false);

% Add guard interval to each OFDM symbol
gi = wlanGIField(nonHTcfg);
numSymbols = size(nonhtdata_nogi, 1) / 80;  % 80 is the number of samples per symbol in 802.11a
nonhtdata = zeros(0, 1);
for k = 1:numSymbols
    symbol_start = (k-1) * 80 + 1;
    symbol_end = k * 80;
    symbol_with_gi = [gi; nonhtdata_nogi(symbol_start:symbol_end)];
    nonhtdata = [nonhtdata; symbol_with_gi];
end


% Concatenate the fields to create the complete waveform
waveform = [lstf; lltf; lsig; nonhtdata];

% Design a custom raised cosine filter
rollOffFactor = 0.5;
span = 8; % Filter span in symbols
sps = 1.5; % Samples per symbol
numCoeffs = span*sps+1; % Number of filter coefficients
coeffs = rcosdesign(rollOffFactor, span, sps, 'sqrt'); % Filter coefficients

% Apply the filter to the waveform (pulse shaping)
shapedSignal = upfirdn(waveform, coeffs, sps);



