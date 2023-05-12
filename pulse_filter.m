%% Coeff for Tx Digital Filters before DAC
HB1 = [-53,0,313,0,-1155,0,4989,8192,4989,0,-1155,0,313,0,-53];
HB2 = [-9,0,73,128,73,0,-9];
HB3 = [1,2,1];
DEC3 = [36,-19,0,-156,-12,0,479,223,0,-1215,-993,0,3569,6277,8192,6277,3569,0,-993,-1215,0,223,479,0,-12,-156,0,-19,36];

% %% Frequency response of HB1
% figure;
% freqz(HB1,8192,'whole');
% 
% %% Frequency response of HB2
% figure;
% freqz(HB2,8192,'whole');
% 
% %% Frequency response of HB3/DEC3
% figure;
% freqz(HB3,8192,'whole');
% 
% figure;
% freqz(DEC3,8192,'whole');

%% Generate a random bit sequence
numBits = 1e5;
data = randi([0 3], numBits, 1);
%% BPSK modulation
modData = pskmod(data, 4, pi/4);
% Design a custom raised cosine filter
rollOffFactor = 0.25;
span = 16; % Filter span in symbols
sps = 4; % Samples per symbol
numCoeffs = span*sps+1; % Number of filter coefficients
coeffs = rcosdesign(rollOffFactor, span, sps, 'sqrt'); % Filter coefficients

modData = upfirdn(modData, coeffs, sps);

%% Apply filters
FilteredData = filter(HB1,1,modData);
FilteredData = filter(HB2,1,FilteredData);
FilteredData = filter(DEC3,1,FilteredData);
% % Time-domain signal
% scatterplot(modData);
% title('Original Signal');
% scatterplot(FilteredData);

% Spectrum analyzer
scope = spectrumAnalyzer;
scope(modData,FilteredData);
release(scope);
