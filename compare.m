% Load the received signals
load('real_received_signal.mat');
load('simulated_received_signal.mat');

% Calculate the Mean Squared Error (MSE)
mse = immse(real_received_signal, simulated_received_signal);
disp(['MSE: ', num2str(mse)]);
