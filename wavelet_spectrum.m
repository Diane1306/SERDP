clc
clear
%%
% Main script
data_dir = '/Users/diane_wt/Downloads/work/PreMean/'; % 2019_East_Tower-fft_20mAGL
towers = {'East', 'West', 'Flux', 'North', 'South_Mobile'};
% Call the function
% [ww20, tt20, CumSec20] = get_data(data_dir, towers, '20m');
[ww10, tt10, CumSec10] = get_data(data_dir, {towers{1:4}}, '10m');
% [ww3, tt3, CumSec3] = get_data(data_dir, towers, '3m');

% Function Definition (for example, if not already implemented)
function [ww, tt, CumSec] = get_data(dr, towers, height)
     % Initialize output cell arrays
    ww = cell(3, 5);
    tt = cell(3, 5);
    CumSec = cell(3, 5);
    
    for ti=1:length(towers)
        % Read the sheets
        df = cell(3,1);
        df{1} = readtable(strcat(dr,height,'/2019_',towers{ti},'_Tower-fft_',height,'AGL.xlsx'), 'Sheet', 'Pre-FFP', 'VariableNamingRule', 'preserve');
        df{2} = readtable(strcat(dr,height,'/2019_',towers{ti},'_Tower-fft_',height,'AGL.xlsx'), 'Sheet', 'FFP', 'VariableNamingRule', 'preserve');
        df{3} = readtable(strcat(dr,height,'/2019_',towers{ti},'_Tower-fft_',height,'AGL.xlsx'), 'Sheet', 'Post-FFP', 'VariableNamingRule', 'preserve');     
        % Loop through each phase (Pre, During, Post)
        for di = 1:3
            % Extract and store the relevant columns
            ww{di, ti} = df{di}.("w' (m/s)");
            tt{di, ti} = df{di}.("t' (C)");
            CumSec{di, ti} = df{di}.("Cum. Sec.");  % Assuming the first column contains this data
        end
        clear df
    end
end

%%
function [cospectrum, freq] = cross_wavelet(a, b, scales)
    fc = centfrq('morl');
    Cg = 0.776; % Torrence et al., 1998
    cwtstruct_temp = cwtft(a, 'wavelet', 'morl', 'scales', scales);
    Wa = cwtstruct_temp.cfs;
    cwtstruct_temp = cwtft(b, 'wavelet', 'morl', 'scales', scales);
    Wb = cwtstruct_temp.cfs;
    freq = cwtstruct_temp.frequencies';
    
    covar = mean((a-mean(a)).*(b-mean(b)));
    cospectrum = squeeze(mean(real(Wa) .* real(Wb) + imag(Wa) .* imag(Wb), 2)) ./ (fc.*Cg.*covar);
end


%% calculate wavelet power spectrum (energy density)
E_f_T_m = nan(length(scales), 3, 3);
E_f_w_m = nan(length(scales), 3, 3);
cospectrum_m = nan(length(scales), 3, 3);
[E_f_T_m(:, :, 3), E_f_w_m(:, :, 3), frequencies, cospectrum_m(:, :, 3), cospectrum_freq] = get_power(ww3, tt3, 5);
[E_f_T_m(:, :, 2), E_f_w_m(:, :, 2), frequencies, cospectrum_m(:, :, 2), cospectrum_freq] = get_power(ww10, tt10, 4);
[E_f_T_m(:, :, 1), E_f_w_m(:, :, 1), frequencies, cospectrum_m(:, :, 1), cospectrum_freq] = get_power(ww20, tt20, 5);

function [E_f_T_m, E_f_w_m, frequencies, cospectrum_m, cospectrum_freq] = get_power(ww, tt, lentowers)
    Cg = pi;
    fc = 0.251;
    dt = 0.1;
    npoints = 18000;
    lenscales = floor(log2(npoints));
    s0 = 2*dt; % following Torrence et al. 1998
    scales = s0 .* 2.^(-3:lenscales);
    
    E_f_T = nan(length(scales), 3, lentowers);
    E_f_w = nan(length(scales), 3, lentowers);
    for ti=1:lentowers
        for di = 1:3
            df = rmmissing(tt{di, ti});
            cwtstruct_temp = cwtft(df, 'wavelet', 'mexh', 'scales', scales);
            W_ab_T = cwtstruct_temp.cfs';
            frequencies = cwtstruct_temp.frequencies';
            E_f_T(:, di, ti) = squeeze(mean(abs(W_ab_T).^2, 1)) ./ (Cg.*fc.*var(df)); % wavelet energy per unit time normalized by variance
            
            df = rmmissing(ww{di, ti});
            cwtstruct_temp = cwtft(df, 'wavelet', 'mexh', 'scales', scales);
            W_ab_w = cwtstruct_temp.cfs';
            E_f_w(:, di, ti) = squeeze(mean(abs(W_ab_w).^2, 1)) ./ (Cg.*fc.*var(df)); % wavelet energy per unit time normalized by variance 
        end
    end
    E_f_T_m = mean(E_f_T, 3);
    E_f_w_m = mean(E_f_w, 3);
    
    cospectrum = nan(length(scales), 3, lentowers);
    for ti=1:lentowers
        for di=1:3
            [cospectrum(:, di, ti), cospectrum_freq] = cross_wavelet(rmmissing(tt{di, ti}),rmmissing(ww{di, ti}), scales);
            
        end
    end
    cospectrum_m = mean(cospectrum, 3);
end


%%
vars = {E_f_T_m, E_f_w_m, cospectrum_m};
ylabels = {"f\timesE_{T'}(f) / T'^2", "f\timesE_{w'}(f) / w'^2", "f\timesCO_{w'T'}(f) / w'T'"};
heights = {'19 m', '10 m', '3 m'};
colors = ['b', 'r', 'k'];
figure('Position', [200 100 1200 800])
for hi=1:3
    for vi=1:3
        subaxis(3, 3, (hi-1)*3+vi, 'sh', 0.032, 'sv', 0.016, 'padding', 0, 'ML', 0.08, 'MB', 0.08, 'MR', 0.07, 'MT', 0.05);
        for di=1:3
            loglog(frequencies, frequencies.*vars{vi}(:, di, hi), 'marker', 'o', 'Color', colors(di), 'LineWidth',2)
            hold on
        end
        if vi<3
            loglog(frequencies, frequencies.^(-2/3) ./ 10^2.5, 'm-.', 'LineWidth',2)
        else
            loglog(frequencies, nan*frequencies.^(-2/3) ./ 10^2.5, 'm-.', 'LineWidth',2)
        end
        hold on
        if vi>2
            loglog(cospectrum_freq, cospectrum_freq.^(-4/3) ./ 10^4, 'c-.', 'LineWidth',2)
        else
            loglog(cospectrum_freq, nan*cospectrum_freq.^(-4/3) ./ 10^4, 'c-.', 'LineWidth',2)
        end
        grid on
        ylim([10^(-5) 10])
        if hi==1 && vi==2
            legend('Pre-FFP', 'FFP', 'Post-FFP', '-2/3', '-4/3', 'box', 'off', 'Orientation','horizontal')
        end
        if vi<3
            xlim([0.0001 2])
        else
            xlim([0.0001 0.4])
        end
        ylim([10^(-4) 1])

        if vi>1
            set(gca,'YTickLabel',[]);
        end

        if hi<3
            set(gca,'XTickLabel',[]);
        else
            xlabel('Frequency [Hz]')
        end

        if hi==2
            ylabel(ylabels{vi})
        end

        if vi==3
            xL=xlim;
            yL=ylim;
            text(xL(2)*5, yL(2)/100, heights{hi},'HorizontalAlignment','right','VerticalAlignment','middle', "FontSize",24,"FontWeight","bold")
        end

    end
end

%%
% flags = cell(3, 4);
% for di=1:3
%     [flags{di, 1}, flags{di, 2}, flags{di, 3}, flags{di, 4}] = get_quadrant(phase{di});
% end
% 
% figure('Position', [300 200 800 600])
% linestyles = {'-o', ':*', '-.d', '--^'};
% di=1;
% for fi=1:4
%     x = cospectrum_freq{di};
%     y = cospectrum_freq{di}.*cospectrum{di}';
%     semilogx(x(flags{di, fi}), y(flags{di, fi}),linestyles{fi}, 'Color',colors{di},'LineWidth',2)
%     hold on
% end
% yline(0, 'Color','k', 'LineStyle','--', 'HandleVisibility','off')
% xlabel('Frequency [Hz]')
% ylabel('fE(f) [m s^{-1} \circC]')
% title('wt Wavelet Co-secptrum @20m of East Tower')
% legend('0', '90', '180', '270', 'box', 'off')

% %%
% function [flag1, flag2, flag3, flag4] = get_quadrant(phi)
%     % Logical flags for each quadrant
%     flag1 = (phi >= -45) & (phi < 45);
%     flag2 = (phi >= 45) & (phi < 135);
%     flag3 = (phi >= 135) | (phi < -135);
%     flag4 = (phi >= -135) & (phi < -45);
% end
