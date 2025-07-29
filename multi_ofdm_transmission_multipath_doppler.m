clc; clear;

% Define static OFDM parameters
numScarr = 64; 
modOrder = 16; 
cpLen = 8;
bitsPerSym = log2(modOrder); 

% Pilot/data arrangement
numPilotSym = 7;          
numSym = 84;              % total OFDM symbols

% Pilot symbol indices evenly spaced across 84 symbols
pilotSymIdx = round(linspace(1, numSym, numPilotSym))
dataSymIdx = setdiff(1:numSym, pilotSymIdx)

% Number of bits
numDataCarr = numScarr; 
numDataBits = numDataCarr * length(dataSymIdx) * bitsPerSym; 

% Generate random data bits & map
dataBits = randi([0, 1], numDataBits, 1); 
dataSym = qammod(dataBits, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);

% Pilot symbols (fixed +1 for all subcarriers)
pilotSym = ones(numScarr, numPilotSym);

% Construct OFDM frame
ofdmFrame = zeros(numScarr, numSym); 
ofdmFrame(:, pilotSymIdx) = pilotSym;
ofdmFrame(:, dataSymIdx) = reshape(dataSym, numScarr, length(dataSymIdx));

% Channel model
sampleRate = 20e6;
maxDoppShift = 20;
pathDelay = [50, 150, 250, 350] * 1e-9;
pathGaindB = [0, -3, -6, -9];
mpChan = comm.RayleighChannel( ...
    "SampleRate", sampleRate, ...
    "AveragePathGains", pathGaindB, ...
    "PathDelays", pathDelay, ...
    "MaximumDopplerShift", maxDoppShift, ...
    "RandomStream", "mt19937ar with seed", ...
    "Seed", 10, ...
    "PathGainsOutputPort", true);

% SNR sweep
snrVec = 0:5:40; 
berVec = zeros(size(snrVec)); 

for snrIdx = 1:length(snrVec)
    snrdB = snrVec(snrIdx);
    
    % OFDM modulation with CP
    ofdmModOut = zeros(numScarr + cpLen, numSym); 
    for n = 1:numSym
        ifftOut = ifft(ofdmFrame(:, n), numScarr); 
        ofdmModOut(:, n) = [ifftOut(end - cpLen + 1:end); ifftOut]; 
    end
    
    % ---- Serialize for Rayleigh channel ----
    txSig = ofdmModOut(:);   % (72*84 x 1)

    % Pass through multipath Rayleigh channel
    [rxChanOut, ~] = mpChan(txSig);

    % Add AWGN
    sigPower = mean(abs(rxChanOut).^2); 
    snrLinear = 10^(snrdB / 10); 
    noisePower = sigPower / snrLinear; 
    noiseStd = sqrt(noisePower / 2); 
    awgnNoise = noiseStd * (randn(size(rxChanOut)) + 1i * randn(size(rxChanOut))); 
    rxSig = rxChanOut + awgnNoise;

    % ---- Reshape back to [72 x 84] ----
    chanOut = reshape(rxSig, numScarr + cpLen, numSym);

    % OFDM Demodulation (remove CP + FFT)
    ofdmDemodOut = zeros(numScarr, numSym); 
    for n = 1:numSym
        ofdmDemodOut(:, n) = fft(chanOut(cpLen + 1:end, n), numScarr); 
    end

    % Channel Estimation (using pilot symbols)
    pilotRx = ofdmDemodOut(:, pilotSymIdx); 
    pilotTx = ofdmFrame(:, pilotSymIdx); 
    pilotCSI = pilotRx ./ pilotTx;  % LS estimate per subcarrier
    
    % Interpolate CSI in time for data symbols
    csiEst = zeros(numScarr, numSym); 
    for k = 1:numScarr
        csiEst(k, :) = interp1(pilotSymIdx, pilotCSI(k, :), 1:numSym, 'linear', 'extrap');
    end

    % Equalization
    eqOut = ofdmDemodOut ./ csiEst; 
    dataSymRx = reshape(eqOut(:, dataSymIdx), [], 1); 
    qamDemodOut = qamdemod(dataSymRx, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);

    % BER calculation
    numBitErrors = nnz(dataBits ~= qamDemodOut); 
    ber = numBitErrors / numDataBits; 
    berVec(snrIdx) = ber;

    fprintf('SNR: %2d dB -> BER: %.5f\n', snrdB, ber);
end

% Plot BER vs SNR
figure;
semilogy(snrVec, berVec, '-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for OFDM (16-QAM, Block Pilots: 7 pilot symbols, total 84 symbols)');
