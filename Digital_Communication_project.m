%% Digital Communication Project

% Load the given Config
config = load("config.mat").config;
config.K = 10^4; % add number of bits (the dimention)
config.pulseshape = 0; % 0 means rect, 1 means sinc. will change in each question.
config.useSynch = 0; % To find the optimal idx with sync-signal("1") or not("0").
config.UpDownConvert = 0; % used in up-conversion in Q6. for simple case ("0").
config.useRxLPF = 0;
%% infobits

infobits = rand(1,config.K)>0.5;
%% Transmitor TX

function ChannelInVec = TX(config, infobits)
    infobits = infobits(:)'; % make sure the input vector is a row vector
    RPSKbits = 1 - 2 * infobits; % BPSK: 0 → +1, 1 → -1
    PulseTrain = upsample(RPSKbits, config.Ts); % Insert Ts-1 zeros
    
    % Choose pulse type
    if config.pulseshape == 0
        TXpulse = config.tpulserect;
    else
        TXpulse = config.tpulsesinc;
    end
    
    % Apply pulse shaping (convolution)
    TxfilterOutput = conv(TXpulse,PulseTrain);
    TxfilterOutput = TxfilterOutput(:)';

    % Silence (Zeros) addition
    paddingBefore = zeros(1, config.Fs);
    paddingAfter = zeros(1, 4 * config.Fs);

    % Add zero padding (1s silence at start, 4s at end)
    ChannelInVec = [paddingBefore, TxfilterOutput, paddingAfter];
    ChannelInVec = ChannelInVec(:)';
    if config.UpDownConvert == 1
        n = 0:length(ChannelInVec)-1; % Create index vector
        ChannelInVec = ChannelInVec .* cos(2 * pi * config.Fc / config.Fs * n);
    end
end
%% Question 2.a (for rect and sinc)

figure;
plot(config.tpulserect);
xlabel('Sample Index');
ylabel('Amplitude');
title('Q2.a: Rect signal');

figure;
plot(config.tpulsesinc);
xlabel('Sample Index');
ylabel('Amplitude');
title('Q2.c: Sinc signal');
grid on;
%% Q2.a,c part 2: Calculate the norm2 of the rect and sinc vector:

norm2_value_rect = sqrt(sum(config.tpulserect.^2));
disp("norm2 value of rect pulse is " + num2str(norm2_value_rect));

norm2_value_sinc = sqrt(sum(config.tpulserect.^2));
disp("norm2 value of sinc is pulse " + num2str(norm2_value_sinc));
%% Q2.b

clear All;

RPSKbits = 1-2*infobits; % BPSK mapping: 0 → +1, 1 → -1
PulseTrain = upsample(RPSKbits,config.Ts); % between each bit we add Ts-1 samples
samples = 32*config.Ts+1;  % number of samples need for 32 info bits (not just 32 samples)
for r_o_s = 0:1
    config.pulseshape = r_o_s;
    if config.pulseshape == 0  % prepare the pulsetype for the convolution:
        TXpulse=config.tpulserect;
    else
        TXpulse=config.tpulsesinc;
    end
    
    figure;
    ax1 = subplot(2,1,1); % 2 rows, 1 column, first subplot
    
    plot(1:samples, PulseTrain(1:samples));
    title('Q2.b Pulse Train');
    xlabel('Sample index');
    ylabel('Amplitude');
    
    TxfilterOutput = conv(TXpulse,PulseTrain);
    
    ax2 = subplot(2,1,2); % 2 rows, 1 column1, second subplot
    plot(1:samples, TxfilterOutput(1:samples));
    title('Q2.b Convolutioned Wave');
    xlabel('Sample index');
    ylabel('Amplitude');
    grid on;
end
%% Q2.3 The Channel TXRX

function ChannelOutVec = ChannelTXRX(config,ChannelInVec)
    % here we are going to first create the gausian noise,
    % add the input and SHALOM AL ISRAEL
    ChannelInVec = ChannelInVec(:)'; % make sure the input is in row vector shape:
    sigma = 10^(-config.snrdB/20);  % the variance - from snr=p/sigma^2
    noise_vec = randn(1,length(ChannelInVec))*sigma;
    ChannelOutVec = ChannelInVec + noise_vec;
    ChannelOutVec = ChannelOutVec(:)';
end
%% Q3 ploting the Channel effect's with different SNRdBs

SNRs_in_dB=[10,20,60,100];
config.pulsetype=1;  % makes sure that the pulse used is a sinc
ChannelInVec=TX(config,infobits);
samples = 4*config.Fs; % in order to see 4 seconds.
for i = 1:4
    figure;
    curr_SNR = SNRs_in_dB(i);
    ax1 = subplot(2,1,1); % 2 rows, 1 column, first subplot
    
    plot(1:samples, ChannelInVec(1:samples));
    title('Q3 - ChannelInVec');
    xlabel('Sample index');
    ylabel('Amplitude');
    
    config.snrdB=curr_SNR;
    ChannelOutVec = ChannelTXRX(config,ChannelInVec);
    
    ax2 = subplot(2,1,2); % 2 rows, 1 column1, second subplot
    plot(1:samples, ChannelOutVec(1:samples));
    title(['Q3 - ChannelOutVec (', int2str(curr_SNR), ')']);
    xlabel('Sample index');
    ylabel('Amplitude');
    linkaxes([ax1, ax2], 'x'); % Sync x-axes
end
%% 2.4 The reciever RX

% This reciever handling also sync-signal matched-filtering and convertion
% down + LPF
function rxbits = RX(config, ChannelOutVec)
    % make sure the input is in a row vector:
    ChannelOutVec = ChannelOutVec(:)';

    % Flip the transmission pulse for matched filtering
    if config.pulseshape == 0
        symbol_filter = flip(config.tpulserect);  % Rectangular pulse
    else
        symbol_filter = flip(config.tpulsesinc);  % Sinc pulse
    end

    % If we convert Down:
    if config.UpDownConvert == 1
        n = 0:length(ChannelOutVec)-1; % Create index vector
        ChannelOutVec = ChannelOutVec .* cos(2 * pi * config.Fc / config.Fs * n); %convert down
        ChannelOutVec = conv(ChannelOutVec,config.RxLPFpulse); % use LPF
        
        % === Channel Estimation for Beta (β) ===
        matched_output = conv(flip(config.synchsymbol),ChannelOutVec);
        Esynch = norm(config.synchsymbol)^2; % L2 norm (sqrt of sum of squares)
        % Find max absolute value of matched filter output
        [max_val,max_idx] = max(abs(matched_output));
        % Compute β/2 estimate
        beta_hat = max_val / Esynch;
        
        % Adjust sign of beta if original max_val was negative
        if matched_output(max_idx) < 0
            beta_hat = -beta_hat;
        end
        % disp("beta/2 is estimates as:"+num2str(beta_hat)); % used in Q6.b
        ChannelOutVec = ChannelOutVec ./ beta_hat; 
    end

    % Choose the optimal index - if we are using synch signal (question 5)
    % or not (previous questions)
    if config.useSynch==1
        % Use matched filter to find the start index:
        Synchfilter = conv(flip(config.synchsymbol),ChannelOutVec);  
        [~,idx_synch] = max(abs(Synchfilter(1:4*config.Fs)));
        opt_idx = idx_synch+1;  % we found that its the optimal index: sync+1
    else
        % if we know the length of silence before the infobits:
        opt_idx = config.Fs + length(symbol_filter);
    end

    indices = round(opt_idx + (0:config.K-1) * config.Ts);
    indices = indices(:);

    % Apply matched filtering (convolution)
    filt_res = conv(symbol_filter(:),ChannelOutVec(:));

    % Extract values at sampling points
    rxbits = filt_res(indices);

    % Decision rule: BPSK demodulation
    rxbits = (rxbits < 0);
    rxbits = rxbits(:)';
end
%% Question 4.a

matched_filter = flip(config.tpulserect);  % Rectangular pulse
start_index = config.Fs + length(matched_filter);

% Apply matched filtering (convolution)
filt_res = filter(matched_filter(:)',1,ChannelOutVec(:)');

% Plot a zoomed-in section around the start index
figure;
subplot(2,1,1);
plot(filt_res(-5000+start_index:start_index+5000));
title('Q4.a - Filtered signal');
xlabel('Bit samples');
ylabel('Values');
% Add a vertical line at 5000 with centered text
xline(5000, '-.r', "optimal point", 'LabelVerticalAlignment', 'middle');
grid on;

% Add a vertical line at 5000 with centered text
subplot(2,1,2);
hold on;
plot(filt_res)
title('Q4.a - Filtered signal (almost full size)');
xlabel('Bit samples');
ylabel('Values');
xline(start_index, '-.r');
hold off;
grid on;



%% Sanity check 4.b

% Use rect pulse
config.pulseshape = 0;

% Set SNR (100 dB)
config.snrdB = 100;

% Transmit signal

ChannelInVec = TX(config, infobits); % Using rect

% Pass through channel with noise
ChannelOutVec = ChannelTXRX(config, ChannelInVec);


% Receive signal
config.useSynch=0;
rxbits = RX(config, ChannelOutVec);

% Compute BER
BER = mean(infobits(:)'~=rxbits(:)');

disp("BER of 'Q4.b sanity check' is "+ num2str(BER));
%% Question 4.b

% load RefInputRect.mat
ref = load("RefInputRect.mat");
ChannelOutVec = ref.ChannelOutVec;
infobits = ref.infobits;
config.K=10^4;

config.useSynch=0; % use basic reciever
rxbits = RX(config, ChannelOutVec);
BER = mean(infobits(:)'~=rxbits(:)');
disp("BER of Q4.b is "+ num2str(BER));
%% Question 4.c

config = load("config.mat").config;
config.K = 10^4; % add number of bits (the dimension)
config.pulseshape = 1; % 1 means sinc (chage to 0 if you want it to be rect)
config.UpDownConvert = 0;
config.useRxLPF = 0;
config.useSynch = 0;

infobits = rand(1,config.K)>0.5;
SNRs_in_dB = -15:2:15;
BER = zeros(1,16);
for i = 1:16
    config.snrdB = SNRs_in_dB(i);
    ChannelInVec = TX(config,infobits);
    
    ChannelOutVec = ChannelTXRX(config,ChannelInVec);
    
    config.useSynch=0; % use basic reciever
    rxbits = RX(config,ChannelOutVec);
    curr_BER = mean(infobits(:)'~=rxbits(:)');
    BER(i)=curr_BER;
end
figure;
plot(SNRs_in_dB,BER,"o-");
title("Q4.c - BER vs SNRs (log scale) - sinc");
ylabel("Bit Error Rate");
xlabel("SNR (in dB)");
set(gca, 'YScale', 'log');  % Maybe you won't see the last BERs because its 0, which isn't defined in log-scale
grid on;

%% Question 5.c - new reciever

% i just updated the reciver to be able to use the sync signal
%% Q5.c check

config = load("config.mat").config;
config.K = 10^4; % add number of bits (the dimention)
config.pulseshape = 1; % 1 means sinc (chage to 0 if you want it to be rect)
config.UpDownConvert = 0;
config.useRxLPF = 0;
infobits = rand(1, config.K) > 0.5;

config.snrdB = 100;
ChannelInVec = TX(config,[config.synchbits(:)' ,infobits(:)']);
ChannelOutVec = ChannelTXRX(config,ChannelInVec);
config.useSynch=1; % use sync reciever
rxbits = RX(config,ChannelOutVec);
BER = mean(infobits(:)~=rxbits(:));
disp("BER of Q5.c is "+ num2str(BER));
%% Q5.d

% load RefInputRect.mat
ref = load("RefInputSynch.mat");
config.pulseshape = 1;
ChannelOutVec = ref.ChannelOutVec;
infobits = ref.infobits(:)';
config.K=10^4;
config.useSynch=1; % use sync reciever
rxbits = RX(config, ChannelOutVec);
BER = mean(infobits(:)~=rxbits(:));
disp("BER of Q5.d is "+ num2str(BER));
%% Q5.e

config = load("config.mat").config;
config.K = 10^4; % add number of bits (the dimention)
config.pulseshape = 1; % 1 means sinc (chage to 0 if you want it to be rect)
config.useSynch=1; % use sync reciever
config.useRxLPF=0;
config.UpDownConvert=0;
infobits = rand(1,config.K)>0.5;
SNRs_in_dB = -15:2:15;
BER = zeros(1,16);
for i = 1:16
    config.snrdB = SNRs_in_dB(i);
    ChannelInVec = TX(config,[config.synchbits(:)' ,infobits(:)']);
    ChannelOutVec = ChannelTXRX(config,ChannelInVec);
    rxbits = RX(config,ChannelOutVec);
    curr_BER = mean(infobits(:)'~=rxbits(:)');
    BER(i)=curr_BER;
end
figure;
plot(SNRs_in_dB,BER,"o-");
title("Q5.e - BER vs SNRs (log scale) - sinc");
ylabel("Bit Error Rate");
xlabel("SNR (in dB)");
set(gca, 'YScale', 'log');  % Maybe you won't see the last BERs because its 0, which isn't defined in log-scale
grid on;
%% Q6.a

config = load("config.mat").config;
config.K = 10^4; % add number of bits (the dimention)
config.pulseshape = 1; % 1 means sinc (chage to 0 if you want it to be rect)
config.UpDownConvert = 1;  % now we do want to use the convert
config.useRxLPF = 1;
infobits = rand(1, config.K) > 0.5;

config.snrdB = 100;
ChannelInVec = TX(config,[config.synchbits(:)' ,infobits]);
ChannelOutVec = ChannelTXRX(config,ChannelInVec);
config.useSynch=1; % use sync reciever
rxbits = RX(config,ChannelOutVec);
BER = mean(infobits(:)'~=rxbits(:)');
disp("BER of Q6.a is "+ num2str(BER));
%% Q6.b

% load RefInputRect.mat
ref = load("RefInputMod.mat");
config.pulseshape = 1;
config.useSynch = 1;
config.UpDownConvert = 1;  % now we do want to use the convert
config.useRxLPF = 1;

ChannelOutVec = ref.ChannelOutVec(:)';
infobits = ref.infobits(:)';
config.K=10^4;
config.snrdB = 100;
config.useSynch=1; % use sync reciever
rxbits = RX(config, ChannelOutVec);
BER = mean(infobits(:)'~=rxbits(:)');
disp("BER of Q6.b is "+ num2str(BER));
%% Q6.c

config.pulseshape = 1; %sinc
config.useSynch = 1;
config.UpDownConvert = 1;
config.useRxLPF = 1;
infobits = rand(1,config.K)>0.5;
config.K=10^4;

SNRs_in_dB = -15:2:15;
BER = zeros(1,16);
for i = 1:16
    config.snrdB = SNRs_in_dB(i);
    ChannelInVec = TX(config,[config.synchbits,infobits]);
    ChannelOutVec = ChannelTXRX(config,ChannelInVec);
    rxbits = RX(config,ChannelOutVec);
    curr_BER = mean(infobits(:)~=rxbits(:));
    BER(i)=curr_BER;
end
figure;
plot(SNRs_in_dB,BER,"o-");
title("Q6.c - BER vs SNRs (log scale) - sinc");
ylabel("Bit Error Rate");
xlabel("SNR (in dB)");
set(gca, 'YScale', 'log');  % Maybe you won't see the last BERs because its 0, which isn't defined in log-scale
grid on;
%% Q7.a
% Initialize flags & values
config.pulsetype = 1;
config.useSynch = 1;
config.UpDownConvert = 1;
config.useRxLPF = 1;
strvec = 'Communication is fun!!!  ';  % exactly 25 characters
infobits = (text2bitstream(strvec))';
config.K = length(infobits);  

infobits = (text2bitstream(strvec))';
ChannelInVec = TX(config,[config.synchbits infobits]);
% ChannelOutVec = ChannelTXRXAudio(config,ChannelInVec); - didn't work on
% my computer.. but the channel i implemented do work:
ChannelOutVec = ChannelTXRX(config,ChannelInVec);
rxbits = RX(config,ChannelOutVec);
strvec = bitstream2text(rxbits)


%% Q7.b

audio = load("RefInputAudio.mat");
config.K = 200;
ChannelOutVec = audio.ChannelOutVec';
rxbits = RX(config,ChannelOutVec);
messege = bitstream2text(rxbits);
disp(messege)
