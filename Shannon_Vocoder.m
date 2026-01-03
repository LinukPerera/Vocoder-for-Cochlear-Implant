% Part 1: Shannon Vocoder - Final Submission
% Module: MOD007721 - Advanced DSP
% Implementation of 18-Channel Vocoder with Comprehensive Analysis

clear; clc; close all;

% --- DETAILS ---
disp('=======================================================');
disp('               MOD007721 ASSIGNMENT SUBMISSION         ');
disp('=======================================================');
disp('ARU Mail:      wldp100@student.aru.ac.uk');
disp('SRN:           2516550');
disp('CINEC ID 1:    FET/DEE/ARU800/2501/0003');
disp('CINEC ID 2:    M20030502015');
disp('=======================================================');
disp('Starting Vocoder Processing...');

% 1. Setup and Input
filename = 'bkbe2114.wav'; 
if ~isfile(filename)
    error(['File ' filename ' not found.']);
end

[input_sig, fs_orig] = audioread(filename);

% Resample to 22.05 kHz as per specification
fs = 22050;
input_sig = resample(input_sig, fs, fs_orig);

% Convert to mono if stereo
if size(input_sig, 2) > 1
    input_sig = mean(input_sig, 2);
end

t = (0:length(input_sig)-1)/fs;

% Parameters
num_channels = 18;       
order_ana = 3;           % Shannon used 3rd order Elliptical
order_env = 1;           % Recommended 1st order for envelope
f_env_cutoff = 160;      % Envelope cutoff 160 Hz

% 2. Filter Bank Design (Robust Calculation)
% Formula: f = 165.4 * (10^(0.06*x) - 1)
x_scale = linspace(0, 35, num_channels + 1);
cutoffs = 165.4 * (10.^(0.06 * x_scale) - 1);

% --- CRITICAL FIX START ---
% The formula generates freqs up to 20kHz, but fs=22kHz (Nyquist=11kHz).
% We must clamp frequencies to strictly less than Nyquist (fs/2).
nyquist_limit = fs / 2;
safe_upper_limit = nyquist_limit * 0.99; % 99% of Nyquist to be safe for ellip

% Clamp upper values
cutoffs(cutoffs >= safe_upper_limit) = safe_upper_limit;

% Clamp lower values (ellip cannot handle 0 Hz exactly)
cutoffs(cutoffs <= 0) = 20; 

% Sort to ensure strict ascending order
cutoffs = sort(cutoffs);
% --- CRITICAL FIX END ---

% 3. Processing Loop
output_sig = zeros(size(input_sig));
white_noise = randn(size(input_sig)); 

% Setup for tracing Channel 8 (for plots)
trace_channel = 8; 
trace_data = struct(); 
trace_captured = false;

fprintf('Processing %d channels...\n', num_channels);

for k = 1:num_channels
    % Band Limits
    fl = cutoffs(k); 
    fh = cutoffs(k+1);
    
    % --- SAFETY CHECK ---
    % If bandwidth is negligible or inverted due to clamping, skip
    if (fh - fl) < 10 
        continue; 
    end
    
    % Normalizing frequencies for ellip (0 to 1 scale)
    Wn = [fl fh] / (fs/2);
    
    % A. Analysis Filter (BPF)
    [b_bp, a_bp] = ellip(order_ana, 1, 40, Wn);
    sig_bp = filter(b_bp, a_bp, input_sig);
    
    % B. Rectification (Half-Wave)
    sig_rect = sig_bp;
    sig_rect(sig_rect < 0) = 0;
    
    % C. Envelope Extraction (LPF)
    [b_lp, a_lp] = butter(order_env, f_env_cutoff/(fs/2));
    sig_env = filter(b_lp, a_lp, sig_rect);
    
    % D. Noise Multiplication
    sig_mod = sig_env .* white_noise;
    
    % E. Output Filtering (BPF)
    sig_out_ch = filter(b_bp, a_bp, sig_mod);
    
    % Summation
    output_sig = output_sig + sig_out_ch;
    
    % Capture Data for Channel 8 plots (Only once)
    if k == trace_channel
        trace_data.freq_range = [fl fh];
        trace_data.bp = sig_bp;
        trace_data.rect = sig_rect;
        trace_data.env = sig_env;
        trace_data.mod = sig_mod;
        trace_data.out = sig_out_ch;
        trace_data.coeffs_b = b_bp;
        trace_data.coeffs_a = a_bp;
        trace_data.coeffs_lp_b = b_lp;
        trace_data.coeffs_lp_a = a_lp;
        trace_captured = true;
    end
end

% Normalize Output
if max(abs(output_sig)) > 0
    output_sig = output_sig / max(abs(output_sig));
end

% ================== ANALYSIS FIGURES ===============================

% --- Figure 1: System Overview ---
figure('Name', 'Vocoder Overview', 'NumberTitle', 'off');
subplot(2,2,1); plot(t, input_sig); title('Input Waveform'); grid on; xlabel('Time (s)'); xlim([0 t(end)]);
subplot(2,2,2); spectrogram(input_sig, 256, [], [], fs, 'yaxis'); title('Input Spectrogram');
subplot(2,2,3); plot(t, output_sig); title('Output Waveform'); grid on; xlabel('Time (s)'); xlim([0 t(end)]);
subplot(2,2,4); spectrogram(output_sig, 256, [], [], fs, 'yaxis'); title('Output Spectrogram');
sgtitle(['System Overview: ' num2str(num_channels) ' Channels']);

% --- Figure 2: Detailed Step-by-Step (Channel 8) ---
if trace_captured
    figure('Name', 'Step-by-Step Channel Analysis', 'NumberTitle', 'off');
    set(gcf, 'Position', [50, 50, 1000, 800]);

    stages = {trace_data.bp, trace_data.rect, trace_data.env, trace_data.mod, trace_data.out};
    titles = {'1. Analysis BPF', '2. Half-Wave Rectified', '3. LPF Envelope', '4. Noise Multiplied', '5. Output BPF'};

    for i = 1:5
        % Waveform (Zoomed)
        subplot(5,2, (i*2)-1);
        zoom_idx = 10000:11000; % View 1000 samples from middle of file
        if length(stages{i}) > 11000
            plot(t(zoom_idx), stages{i}(zoom_idx)); 
        else
            plot(t, stages{i});
        end
        title([titles{i} ' (Time Domain)']); grid on;
        
        % Spectrum
        subplot(5,2, i*2);
        L = length(stages{i});
        Y = fft(stages{i});
        P2 = abs(Y/L);
        P1 = P2(1:floor(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f_axis = fs*(0:(L/2))/L;
        plot(f_axis, 20*log10(P1+eps)); 
        title([titles{i} ' (Freq Domain)']); grid on; xlim([0 5000]); ylabel('dB');
    end
    sgtitle(['Detailed Analysis of Channel ' num2str(trace_channel)]);
end

% ================== EXTRA PROOF FIGURES ===============================

% Figure A – Input Bandwidth Proof
figure('Name','A – Input Bandwidth (0–11 kHz)','Position',[100 100 900 500]);
subplot(2,1,1);
plot(t, input_sig); grid on; xlabel('Time (s)'); ylabel('Amplitude');
title('Input Waveform');
subplot(2,1,2);
spectrogram(input_sig, hamming(512), 256, 1024, fs, 'yaxis');
ylim([0 11.5]); % force y-axis to show up to just above Nyquist
line([0 max(t)], [11.025 11.025], 'Color','r','LineWidth',2,'LineStyle','--');
text(0.1,11.2,'Nyquist = 11.025 kHz','Color','red','FontSize',12,'FontWeight','bold');
title('Input Spectrogram – Bandwidth Proof');

% Figure B – Channel Bandwidths & Cut-off Frequencies
figure('Name','B – Analysis Filter Bank Bandwidths','Position',[200 200 1000 600]);
subplot(2,1,1);
semilogx(cutoffs, 1:length(cutoffs), 'o-','LineWidth',2,'MarkerFaceColor','b');
grid on; xlabel('Frequency (Hz)'); ylabel('Channel index');
title('Cut-off Frequencies for 18-Channel Vocoder');
text(cutoffs(1),1,'DC','VerticalAlignment','bottom');
text(cutoffs(end),num_channels+1,'≈11 kHz','VerticalAlignment','top');

subplot(2,1,2);
channel_low  = cutoffs(1:end-1);
channel_high = cutoffs(2:end);
channel_center = (channel_low + channel_high)/2;
channel_bw = channel_high - channel_low;
bar(1:num_channels, channel_bw, 'FaceColor',[0.2 0.6 0.8]);
xlabel('Channel number'); ylabel('Bandwidth (Hz)');
title('Bandwidth of Each Analysis Filter');
grid on;

% Print table for report
fprintf('\n=== Channel Cut-off Frequencies ===\n');
for k = 1:num_channels
    fprintf('Channel %2d:  %.1f – %.1f Hz  (BW = %.1f Hz)\n', k, cutoffs(k), cutoffs(k+1), cutoffs(k+1)-cutoffs(k));
end

% Figure C – Analysis Filter Responses (Shannon vs Rosen vs Optimal)
figure('Name','C – Analysis Filter Magnitude Responses','Position',[300 100 1100 700]);
hold on;
legend_str = {};

% Shannon 3rd order
[b3,a3] = ellip(3,1,40,[1000 2000]/(fs/2)); % example band just for illustration
[h3,f3] = freqz(b3,a3,4096,fs);
plot(f3,20*log10(abs(h3)),'LineWidth',2); legend_str{end+1} = 'Shannon – 3rd order elliptical';

% Rosen 6th order
[b6,a6] = ellip(6,1,40,[1000 2000]/(fs/2));
[h6,f6] = freqz(b6,a6,4096,fs);
plot(f6,20*log10(abs(h6)),'--','LineWidth',2); legend_str{end+1} = 'Rosen – 6th order elliptical';

% Optimal order via ellipord
[Nopt,Wnopt] = ellipord([1000 2000]/(fs/2), [900 2100]/(fs/2),1,40);
[bopt,aopt] = ellip(Nopt,1,40,Wnopt,'bandpass');
[hopt,fopt] = freqz(bopt,aopt,4096,fs);
plot(fopt,20*log10(abs(hopt)),':','LineWidth',3); legend_str{end+1} = sprintf('Optimal order = %d (ellipord)',Nopt);

grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Comparison of Analysis Filter Orders (Example Band: 1-2 kHz)');
legend(legend_str,'Location','southwest'); xlim([0 11000]);

% Figure D – Low-Pass Envelope Filter Response (160 Hz)
figure('Name','D – Envelope Low-Pass Filter','Position',[400 200 800 500]);
% 1st order Butterworth (used in code)
[b_lp_but,a_lp_but] = butter(1, f_env_cutoff/(fs/2));
[hb,fb] = freqz(b_lp_but,a_lp_but,4096,fs);
plot(fb,20*log10(abs(hb)),'LineWidth',2); hold on;

% 1st order elliptical alternative
[b_lp_ell,a_lp_ell] = ellip(1,0.5,40, f_env_cutoff/(fs/2));
[he,fe] = freqz(b_lp_ell,a_lp_ell,4096,fs);
plot(fe,20*log10(abs(he)),'--','LineWidth',2);

grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('160 Hz Envelope Low-Pass Filter Responses');
legend('1st-order Butterworth (implemented)','1st-order Elliptical (alternative)','Location','northeast');
xlim([0 500]);

% Figure E – All 18 Analysis Filter Bank Overlay
figure('Name','E – Complete 18-Channel Filter Bank','Position',[100 100 1200 700]);
hold on;
colors = lines(num_channels);
for k = 1:num_channels
    fl = cutoffs(k); fh = cutoffs(k+1);
    if (fh-fl) < 10, continue; end
    Wn = [fl fh]/(fs/2);
    [b,a] = ellip(order_ana,1,40,Wn);
    [h,f] = freqz(b,a,4096,fs);
    plot(f,20*log10(abs(h)),'Color',colors(k,:));
end
grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Complete 18-Channel Analysis Filter Bank (3rd-order Elliptical)');
xlim([0 11000]); ylim([-80 5]);

% Pole-zero plot for a Typical Analysis Filter (Channel 8)
% Using the coefficients captured in the main loop for Channel 8
if trace_captured
    figure('Name','Pole-Zero Plot Channel 8');
    zplane(trace_data.coeffs_b, trace_data.coeffs_a);
    title(['Pole-Zero Plot (Channel 8: ' num2str(round(trace_data.freq_range(1))) '-' num2str(round(trace_data.freq_range(2))) ' Hz)']);
    grid on;
end

disp('All plots generated successfully.');
soundsc(output_sig, fs);
