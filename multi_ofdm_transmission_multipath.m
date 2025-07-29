clc; clear;

% Define static OFDM parameters
numScarr = 72; 
modOrder = 16; 
cpLen = 24; 
bitsPerSym = log2(modOrder); 
numSym = 14; 
numPilotCarr = 8; 
numDataCarr = 64; 
numDataBits = numDataCarr*numSym*bitsPerSym; 
numPilotBits = numPilotCarr*numSym*1;

% Indexing
pilotIdx = 1:9:72; 
dataIdx = setdiff(1:numScarr, pilotIdx); 

% Generate random bits (same for all SNR values to ensure fair comparison)
dataBits = randi([0,1], numDataBits, 1); 
PilotBits = randi([0,1],numPilotBits, 1); 

% Generate symbols
dataSym = qammod(dataBits, modOrder, 'InputType','bit', 'UnitAveragePower', true); 
pilotSym = 2*PilotBits - 1; 

% Construct OFDM frame
ofdmFrame = zeros(numScarr, numSym); 
ofdmFrame(pilotIdx, :) = reshape(pilotSym, numPilotCarr, numSym); 
ofdmFrame(dataIdx, :) = reshape(dataSym, numDataCarr, numSym); 

% Channel model
mpChan = [0.8; zeros(2,1); -0.5; zeros(3,1); -0.34];
mpChan = mpChan / norm(mpChan);

% SNR sweep
snrVec = 0:5:100; 
berVec = zeros(size(snrVec)); 

for snrIdx = 1:length(snrVec)
    snrdB = snrVec(snrIdx);
    
    % OFDM modulation with CP
    ofdmModOut = zeros(numScarr+cpLen, numSym); 
    for n=1:numSym
        ifftOut = ifft(ofdmFrame(:, n), numScarr); 
        ofdmModOut(:, n) = [ifftOut(end-cpLen+1:end); ifftOut]; 
    end
    
    % Pass through multipath channel
    mpChanOut = zeros(size(ofdmModOut));
    for n=1:numSym
        mpChanOut(:, n) = filter(mpChan, 1, ofdmModOut(:, n)); 
    end

    % Add AWGN
    sigPower = mean(abs(mpChanOut).^2, 'all'); 
    snrLinear = 10^(snrdB/10); 
    noisePower = sigPower / snrLinear; 
    noiseStd = sqrt(noisePower/2); 
    awgn = noiseStd * (randn(size(mpChanOut)) + 1i * randn(size(mpChanOut))); 
    chanOut = mpChanOut + awgn;

    % OFDM Demodulation (remove CP + FFT)
    ofdmDemodOut = zeros(numScarr, numSym); 
    for n=1:numSym
        ofdmDemodOut(:, n) = fft(chanOut(cpLen+1:end, n), numScarr); 
    end

    % LS Channel Estimation using pilots
    pilotRx = ofdmDemodOut(pilotIdx, :); 
    pilotTx = ofdmFrame(pilotIdx, :); 
    pilotCSI = pilotRx ./ pilotTx; 
    csiEst = zeros(numScarr, numSym); 
    csiEst(pilotIdx,:) = pilotCSI;

    % Interpolate CSI at data subcarriers
    for n=1:numSym
        csiEst(:, n) = interp1(pilotIdx, csiEst(pilotIdx, n), (1:numScarr)', 'spline', 'extrap');
        csiEst(pilotIdx, n) = pilotCSI(:, n); % restore actual pilot CSI
    end

    % Equalization
    eqOut = ofdmDemodOut ./ csiEst; 
    dataSymRx = reshape(eqOut(dataIdx, :), [], 1); 
    qamDemodOut = qamdemod(dataSymRx, modOrder, 'OutputType','bit', 'UnitAveragePower', true);

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
title('BER vs SNR for OFDM (16-QAM, LS Estimation)');
