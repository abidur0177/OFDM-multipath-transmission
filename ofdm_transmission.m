% ------- First we will generate the QAM modulated symbols -------

% Modulation order
modOrder = 16;

% Bits per symbol calculation
bitsPerSymbol = log2(modOrder);

% Number of subcarriers
numCarr = 8192;

% Define the number of total bits to send
numBits = numCarr * bitsPerSymbol;

% Generate random bits with a length of the number of bits
if exist("numBits","var")
    % Generate random bits
    srcBits = randi([0,1], numBits, 1);
    % QAM modulation
    qamModOut = qammod(srcBits, modOrder, "InputType","bit", "UnitAveragePower",true);
    fprintf("Length of the QAM modulated signal: %d\n", length(qamModOut))
end
% Observe the constellation mapping
%scatterplot(qamModOut)

% ------- Our QAM symbols are ready -------





% ------- Generate the OFDM signal in time domain -------

% Length of the cyclic prefix
cpLen = 32;

% Create the OFDM data
ofdmModOut = ofdmmod(qamModOut, numCarr, cpLen);
fprintf("Length of the OFDM signal: %d\n", length(ofdmModOut))

% Observe the time domain OFDM signal with 
figure;
plot(real(ofdmModOut), 'b');
hold on;
plot(imag(ofdmModOut), 'r');
title("OFDM signal (Real and Imaginary)");
xlabel("Sample index");
ylabel("Amplitude");
legend("Real", "Imaginary");
grid on;

% -------- OFDM symbols are ready to transmit over the channel --------





% -------- Model the channel and impose on the transmitted signal -------

% Define the multipath channel impulse response
mpChan = [0.8; zeros(7,1); -0.5; zeros(7,1); -0.34];

% Define the SNR value
SNR = 15;

if exist("ofdmModOut", "var")
    % Filter the transmitted signal with the channel impulse response
    mpChanOut = filter(mpChan, 1, ofdmModOut);
    % AWGN noise is added on the filtered signal
    chanOut = awgn(mpChanOut, SNR, "measured");
    fprintf("Length of the received signal: %d\n", length(chanOut))
end

% ------- ChanOut is the received signal after noise addition -------





% ------- Demodulate the signal at the receiver side -------

% Discard CP and perform OFDM
ofdmDemodOut = ofdmdemod(chanOut, numCarr, cpLen, cpLen);
fprintf("Length of the demodulated signal: %d\n", length(ofdmDemodOut))
disp(size(ofdmDemodOut))

% Observe the demodulated signal on the constellation map
%scatterplot(ofdmDemodOut)

% ------- Demodulation is done -------





% ------- Now do the equalization -------

% Calculate the fft of the multipath channel
mpChanFreq = fftshift(fft(mpChan, numCarr));
disp(size(mpChanFreq))

% Equalize the signal in the frequency domain by doing dot division
eqOut = ofdmDemodOut ./ mpChanFreq;
fprintf("Length of the equalized signal: %d\n", length(eqOut))

% Observe the equalized signal in constellation map
scatterplot(eqOut)

if exist("eqOut", "var")
    % Recover the bits
    qamDemodOut = qamdemod(eqOut, modOrder, "OutputType","bit", "UnitAveragePower", true);
    disp(size(qamDemodOut))
    % Calculate the total error bits
    numBitErrors = nnz(srcBits ~= qamDemodOut);
    % Bit error rate (BER) calculation
    fprintf("BER of the system: %d\n", (numBitErrors / numBits)*100)
end

% ------- End of the OFDM transmission -------





% ------- Inspect the signal in terms of spectrum analysing -------

% Transmitted OFDM signal
tx_ofdm = ofdmModOut(cpLen+1:end);
% Calculate the fft of the transmitted OFDM signal
fft_tx_ofdm = fftshift(fft(tx_ofdm, numCarr));


% Define the frequency axis
freqAxis = linspace(-numCarr/2, numCarr/2, numCarr)/1e3;

% Compare the PSD of the received and eqalized signal respectively
figure;
plot(freqAxis, 20*log10(abs(fftshift(ofdmDemodOut))), 'r');
hold on;
plot(freqAxis, 20*log10(abs(fftshift(eqOut))), '--b');
xlabel("Frequency in KHz");
ylabel("Magnitude of demodulated and equalized signal");
legend("Demodulated","Equalized");
title("Comparison of Demodulated and Eqaulized OFDM Signals");
grid on;

% Compare the PSD of the transmitted signal and received signal
figure;
plot(freqAxis, 20*log10(abs(fft_tx_ofdm)), 'g');
hold on;
plot(freqAxis, 20*log10(abs(fftshift(ofdmDemodOut))), '--r');
xlabel("Frequency in KHz");
ylabel("Magnitude of transmitted and demodulated signal");
legend("Transmitted","Demodulated");
title("Comparison of Transmitted and Demodulated OFDM Signals");
grid on;

% Compare the PSD of the transmitted signal and equalized signal
figure;
plot(freqAxis, 20*log10(abs(fftshift(eqOut))), '--b');
hold on;
plot(freqAxis, 20*log10(abs(fft_tx_ofdm)), 'g');
xlabel("Frequency in KHz");
ylabel("Magnitude of equalized and transmitted signal");
legend("Equalized","Transmitted");
title("Comparison of Equalized and Transmitted OFDM Signals");
grid on;