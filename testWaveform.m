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
% Input an image file and convert to binary stream
fileTx = 'apply-hi.jpg';                   % Image file name
fData = imread(fileTx);                          % Read image data from file
scale = 0.2;                    % Image scaling factor
origSize = size(fData);                          % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:);               % Resize image
imsize = size(fData);                            % Store new image size
txImage = fData(:);

% Plot transmit image
% imFig.Visible = 'on';
% subplot(211);
% imshow(fData);
% title('Transmitted Image');
% subplot(212);
% title('Received image appears here...');
% set(gca,'Visible','off');
%% Fragment Transmit Data
msduLength = 2304; % MSDU length in bytes
numMSDUs = ceil(length(txImage)/msduLength);
padZeros = msduLength-mod(length(txImage),msduLength);
txData = [txImage;zeros(padZeros,1)];
txDataBits = double(int2bit(txData,8,false));
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

% Generate baseband NonHT packets separated by idle time
txWaveform = wlanWaveformGenerator(data,nonHTcfg, ...
    NumPackets=numMSDUs,IdleTime=20e-6, ...
    ScramblerInitialization=scramblerInitialization,...
    OversamplingFactor=osf);

waveform = [];
for i = 1:numMSDUs
    % Generate L-STF, L-LTF, and L-SIG fields
    lstf = wlanLSTF(nonHTcfg, overSamplingFactor=osf);
    lltf = wlanLLTF(nonHTcfg, overSamplingFactor=osf);
    lsig = wlanLSIG(nonHTcfg, overSamplingFactor=osf);
    
    % Generate the non-HT Data field for each packet and concatenate with idle time
    psdu = data((i-1)*lengthMPDU*8+1:lengthMPDU*8*i);
    nonhtdata = wlanNonHTData(psdu, nonHTcfg, scramblerInitialization(i), OversamplingFactor = osf);
    
    % Concatenate the fields to create the complete waveform
    waveform = [waveform; lstf; lltf; lsig; nonhtdata];
    
    % Add idle time
    waveform = [waveform; zeros(idleTimeSamples, 1)];
end

diff = abs(txWaveform - waveform);
% figure;
nonzero = find(diff);

figure;
plot(real(txWaveform));
hold on;
plot(real(waveform));
plot(nonzero, real(txWaveform(nonzero)), 'ro');
plot(nonzero, real(waveform(nonzero)), 'ko');
hold off;




