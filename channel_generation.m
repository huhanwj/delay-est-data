clear;clc;



N_cp = 32;% length of CP



Half_Pulse_len = 8;% half length of pulse



Ts = 50;% sampling period

rolloff = 0.05;




N_sc = 64;% number of subcarriers

sc_index_used = [2:27, 39:64];% indices of used subcarriers


F = fft(eye(N_sc))/sqrt(N_sc);% normalized DFT matrix


[lltfLower,lltfUpper] = wlan.internal.lltfSequence();
LLTF = [zeros(6,1); lltfLower; 0; lltfUpper; zeros(5,1)];
sym_fre = LLTF([N_sc/2+1:N_sc,1:N_sc/2]);
sym_time = F'*sym_fre;
cp_sym_time = [sym_time(end-N_cp+1:end);sym_time];
cp_sym_time = [cp_sym_time;zeros(100,1)];% transmitted signal



L = 1;% number of paths
SNR = 20;



delay_ns = [24];% path delay

     
pathDelays = delay_ns*1e-9;
        

avgPathGains = SNR;


fs = 20000000;% sampling frequency
tgnChannel = comm.RayleighChannel('SampleRate',fs,'PathDelays',pathDelays, 'AveragePathGains',avgPathGains);
tgnChannel.PathGainsOutputPort = 1;
tgnChannel.NormalizePathGains = 0;

reset(tgnChannel); 
[~,PATHGAIN] = tgnChannel(cp_sym_time);


h_fit = zeros(N_sc,1);% channel impulse response

for ll = 1:L
    
    starting_time = floor(delay_ns(ll)/Ts)+1;

    delay_mod = mod(delay_ns(ll),Ts);

    delay_index = single(delay_ns(ll)/Ts+1);

    if delay_mod ~= 0
        indices = floor(delay_index)-Half_Pulse_len+1:ceil(delay_index)+Half_Pulse_len-1;
    else
        indices = floor(delay_index)-Half_Pulse_len+1:ceil(delay_index)+Half_Pulse_len;
    end

    h_fit(starting_time:starting_time+2*Half_Pulse_len-1) = h_fit(starting_time:starting_time+2*Half_Pulse_len-1) + (PATHGAIN(1,ll) * raisedcosine(indices-delay_index,rolloff)).';
    
end



H_noiseless = F(sc_index_used,:)*h_fit*sqrt(N_sc);% CSI

sigma2 = 1;
        
noise = normrnd(0,sqrt(sigma2/2),length(H_noiseless),1)+1j*normrnd(0,sqrt(sigma2/2),length(H_noiseless),1);
        
H = H_noiseless + noise;



