%% Generate a random bit sequence
numBits = 1e4;
data = randi([0 1], numBits, 1);
%% BPSK modulation
modData = pskmod(data, 2);

%% Rayleigh and Rician channel
rayleighChannel = comm.RayleighChannel;
ricianChannel = comm.RicianChannel;
%% Pass the data
rayleighOutput = rayleighChannel(modData);
ricianOutput = ricianChannel(modData);

%% Comparison
[Pxx_input, F_input] = periodogram(modData);
[Pxx_Rayleigh, F_Rayleigh] = periodogram(rayleighOutput);
[Pxx_Rician, F_Rician] = periodogram(ricianOutput);

% Plot
% figure; 
% plot(F_input, 10*log10(Pxx_input), F_Rayleigh, 10*log10(Pxx_Rayleigh),F_Rician,10*log10(Pxx_Rician));
% xlabel("Normalized frequency");
% ylabel("PSD(dB)");
% legend("Input","Rayleigh","Rician");
t = (0:length(modData)-1)';

plot(t, rea(modData));
hold on;
plot(t, real(rayleighOutput));
plot(t, real(ricianOutput));
legend("Input","Rayleigh","Rician");

hold off;