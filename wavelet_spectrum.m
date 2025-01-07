% Main script
% data_dir = '/home/msuclass9/study/Diane/SEDRP/data/PeriodMean/';  % Define the data directory
data_dir = '/Users/diane_wt/Downloads/2019_East_Tower-fft_20mAGL.xlsx';

% Call the function
[ww20, tt20, Sweep20, Ejection20, Outward20, Inward20, CumSec20] = get_data(data_dir);

% Function Definition (for example, if not already implemented)
function [ww, tt, Sweep, Ejection, Outward, Inward, CumSec] = get_data(dr)
     % Initialize output cell arrays
    ww = {};
    tt = {};
    Sweep = {};
    Ejection = {};
    Outward = {};
    Inward = {};
    CumSec = {};
    
    % Read the sheets
    df_pre = readtable(dr, 'Sheet', 'Pre-FFP', 'VariableNamingRule', 'preserve');
    df_dur = readtable(dr, 'Sheet', 'FFP', 'VariableNamingRule', 'preserve');
    df_post = readtable(dr, 'Sheet', 'Post-FFP', 'VariableNamingRule', 'preserve');
    
    % Combine the sheets into a cell array
    df = {df_pre, df_dur, df_post};
    
    % Loop through each phase (Pre, During, Post)
    for di = 1:3
        % Extract and store the relevant columns
        ww{end+1} = df{di}.("w' (m/s)");
        tt{end+1} = df{di}.("t' (C)");
        Sweep{end+1} = df{di}.("w't' (Sweep)");
        Ejection{end+1} = df{di}.("w't' (Ejection)");
        Outward{end+1} = df{di}.("w't' (Out Int)");
        Inward{end+1} = df{di}.("w't' (In Int)");
        CumSec{end+1} = df{di}.("Cum. Sec.");  % Assuming the first column contains this data
    end
end

%%
function [cospectrum, phase, freq] = cross_wavelet(a, b, scales)
    fc = centfrq('morl');
    cwtstruct_temp = cwtft(a, 'wavelet', 'morl', 'scales', scales);
    Wa = cwtstruct_temp.cfs;
    cwtstruct_temp = cwtft(b, 'wavelet', 'morl', 'scales', scales);
    Wb = cwtstruct_temp.cfs;
    freq = cwtstruct_temp.frequencies;

    cospectrum = squeeze(mean(real(Wa) .* real(Wb) + imag(Wa) .* imag(Wb), 2)) ./ (fc);
    quadrature = squeeze(imag(Wa) .* real(Wb) - real(Wa) .* imag(Wb)) ./ (fc);
    phase = atan2(quadrature, cospectrum);
end
%% calculate wavelet power spectrum (energy density)
Cg = pi;
fc = 0.251;
dt = 0.1;
npoints = 18000;
lenscales = floor(log2(18000));
scales = 2.^(-3:lenscales);

E_f_T = nan(length(scales), 3);
for di = 1:3
    df = tt20{di}(1:end-7);
    cwtstruct_temp = cwtft(df, 'wavelet', 'mexh', 'scales', scales);
    W_ab_T = cwtstruct_temp.cfs';
    frequencies = cwtstruct_temp.frequencies;
    E_f_T(:, di) = squeeze(mean(abs(W_ab_T).^2, 1)) ./ (Cg.*fc); % wavelet energy per unit time
end

%% 
figure('Position', [300 200 800 600])
colors = ['b', 'r', 'k'];
for di=1:3
    loglog(frequencies, frequencies.*E_f_T(:,di)', colors(di), 'LineWidth',2)
    hold on
end
loglog(frequencies, frequencies.^(-2/3) ./ 10^2, 'm-.', 'LineWidth',2)
grid on
ylim([10^(-5) 10])
xlabel('Frequency [Hz]')
ylabel('fE(f) [\circC^2]')
% ylabel('fE(f) [m^2 s^{-2}]')
title('Wavelet Energy Density @20m of East Tower')
legend('Pre-FFP', 'FFP', 'Post-FFP', '-2/3', 'box', 'off')

%% Cospectrum
cospectrum = cell(3, 1);
quadrature = cell(3, 1);
cospectrum_freq = cell(3, 1);
phase = cell(3, 1);
for di=1:3
    [cospectrum{di}, phase{di}, cospectrum_freq{di}] = cross_wavelet(tt20{di}(1:end-7),ww20{di}(1:end-7), scales);
    % [wcoh,wcs,f] = wcoherence(tt20{di}(1:end-7),ww20{di}(1:end-7),10);
    % cospectrum{di} = squeeze(mean(real(wcs), 2)) ./ (fc.*Cg);
    % quadrature{di} = squeeze(mean(imag(wcs), 2)) ./ (fc.*Cg);
    % cospectrum_freq{di} = f;
    % phase{di} = atan2(quadrature{di}, cospectrum{di});
end

%%
figure('Position', [300 200 800 600])
colors = ['b', 'r', 'k'];
for di=1:3
    loglog(cospectrum_freq{di}, cospectrum_freq{di}.*cospectrum{di}', colors(di), 'LineWidth',2)
    hold on
end
loglog(cospectrum_freq{1}, cospectrum_freq{1}.^(-4/3) ./ 10^3, 'm-.', 'LineWidth',2)
grid on
ylim([10^(-5) 10])
xlabel('Frequency [Hz]')
ylabel('fE(f) [m s^{-1} \circC]')
title('wt Wavelet Co-secptrum @20m of East Tower')
legend('Pre-FFP', 'FFP', 'Post-FFP', '-4/3', 'box', 'off')