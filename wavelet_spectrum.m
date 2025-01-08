clc
clear
%%
% Main script
data_dir = '/Users/diane_wt/Downloads/work/PreMean/'; % 2019_East_Tower-fft_20mAGL
heights = ['3m', '10m', '20m'];
towers = ['East', 'West', 'Flux', 'South_Mobile', 'North'];
% Call the function
[ww20, tt20, CumSec20] = get_data(data_dir);

% Function Definition (for example, if not already implemented)
function [ww, tt, CumSec] = get_data(dr)
     % Initialize output cell arrays
    ww = cell{3, 4};
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
        CumSec{end+1} = df{di}.("Cum. Sec.");  % Assuming the first column contains this data
    end
end

%%
function [cospectrum, phase, freq] = cross_wavelet(a, b, scales)
    fc = centfrq('morl');
    Cg = 0.776; % Torrence et al., 1998
    cwtstruct_temp = cwtft(a, 'wavelet', 'morl', 'scales', scales);
    Wa = cwtstruct_temp.cfs;
    cwtstruct_temp = cwtft(b, 'wavelet', 'morl', 'scales', scales);
    Wb = cwtstruct_temp.cfs;
    freq = cwtstruct_temp.frequencies;

    cospectrum = squeeze(mean(real(Wa) .* real(Wb) + imag(Wa) .* imag(Wb), 2)) ./ (fc.*Cg);
    quadrature = squeeze(mean(imag(Wa) .* real(Wb) - real(Wa) .* imag(Wb), 2)) ./ (fc.*Cg);
    phase = atan2(quadrature, cospectrum).*(180/pi());
end

%%
function [flag1, flag2, flag3, flag4] = get_quadrant(phi)
    % Logical flags for each quadrant
    flag1 = (phi >= -45) & (phi < 45);
    flag2 = (phi >= 45) & (phi < 135);
    flag3 = (phi >= 135) | (phi < -135);
    flag4 = (phi >= -135) & (phi < -45);
end

%% calculate wavelet power spectrum (energy density)
Cg = pi;
fc = 0.251;
dt = 0.1;
npoints = 18000;
lenscales = floor(log2(18000));
s0 = 2*dt; % following Torrence et al. 1998
scales = s0 .* 2.^(-3:lenscales);

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
colors = {'b', 'r', 'k'};
for di=1:3
    loglog(cospectrum_freq{di}, cospectrum_freq{di}.*cospectrum{di}', colors{di}, 'LineWidth',2)
    hold on
end
loglog(cospectrum_freq{1}, cospectrum_freq{1}.^(-4/3) ./ 10^3, 'm-.', 'LineWidth',2)
grid on
ylim([10^(-5) 10])
xlabel('Frequency [Hz]')
ylabel('fE(f) [m s^{-1} \circC]')
title('wt Wavelet Co-secptrum @20m of East Tower')
legend('Pre-FFP', 'FFP', 'Post-FFP', '-4/3', 'box', 'off')

%%
flags = cell(3, 4);
for di=1:3
    [flags{di, 1}, flags{di, 2}, flags{di, 3}, flags{di, 4}] = get_quadrant(phase{di});
end

figure('Position', [300 200 800 600])
linestyles = {'-o', ':*', '-.d', '--^'};
di=1;
for fi=1:4
    x = cospectrum_freq{di};
    y = cospectrum_freq{di}.*cospectrum{di}';
    semilogx(x(flags{di, fi}), y(flags{di, fi}),linestyles{fi}, 'Color',colors{di},'LineWidth',2)
    hold on
end
yline(0, 'Color','k', 'LineStyle','--', 'HandleVisibility','off')
xlabel('Frequency [Hz]')
ylabel('fE(f) [m s^{-1} \circC]')
title('wt Wavelet Co-secptrum @20m of East Tower')
legend('0', '90', '180', '270', 'box', 'off')
