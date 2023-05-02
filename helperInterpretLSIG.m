function [MCS,PSDULength,numRxSamples] = helperInterpretLSIG(recLSIGBits,sr)
% helperInterpretLSIG Interprets recovered L-SIG bits
%
%   [MCS,PSDULENGTH,NUMRXSAMPLES] = helperInterpretLSIG(RECLSIGBITS,SR)
%   returns the modulation and coding scheme, PSDU length and number of
%   samples within the packet given the recovered L-SIG bits and sample
%   rate in Hertz.   

% Copyright 2015 The MathWorks, Inc.

%#codegen

% Rate and length are determined from bits
rate = double(recLSIGBits(1:3));
length = double(recLSIGBits(5+(1:12)));

% MCS rate table for 802.11a
R = wlan.internal.nonHTRateSignalBits();
MCS = find(all(bsxfun(@eq,R(1:3,:),rate)))-1;
PSDULength = bi2de(length.');

% Get the indices for the non-HT data field (assuming 20MHz bandwidth)
cfgNonHT = wlanNonHTConfig('MCS',MCS,'PSDULength',PSDULength);    
nonHTDataInd = wlanFieldIndices(cfgNonHT,'NonHT-Data');

% Calculate the number of samples given the actual sampling rate
numRxSamples = double(nonHTDataInd(2))*sr/20e6;

end