clc
clear
%%
% Main script
data_dir = '/Users/diane_wt/Downloads/work/PreMean/'; % 2019_East_Tower-fft_20mAGL
towers = {'Flux', 'North', 'West', 'East', 'South_Mobile'};
% Call the function
[ww20, tt20, CumSec20] = get_data(data_dir, towers, '20m');
[ww10, tt10, CumSec10] = get_data(data_dir, towers(1:4), '10m');
[ww3, tt3, CumSec3] = get_data(data_dir, towers, '3m');

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

%% create standard CumSec
CumSec_n = CumSec20;
CumSec_n{3, 5} = CumSec3{3, 5};
time_temp = CumSec20{2, 5}(1):0.1:CumSec20{2, 5}(1)+(12000-1)*0.1;
CumSec_n{2, 5} = time_temp';

%%
function [cospectrum, freq] = cross_wavelet(a, b, scales)
    Cg = 0.776; % Torrence et al., 1998
    cwtstruct_temp = cwtft(a, 'wavelet', 'morl', 'scales', scales);
    Wa = cwtstruct_temp.cfs;
    cwtstruct_temp = cwtft(b, 'wavelet', 'morl', 'scales', scales);
    Wb = cwtstruct_temp.cfs;
    freq = cwtstruct_temp.frequencies';
    fc = freq(1).*scales(1);

    covar = mean((a-mean(a)).*(b-mean(b)));
    cospectrum = squeeze(mean(real(Wa) .* real(Wb) + imag(Wa) .* imag(Wb), 2)) ./ (fc.*Cg.*covar);
end

%%
% Helper function to process a single array
function [filteredData, originalIndices] = processArray(data)
    % Identify missing data (assume missing data is NaN)
    missingIdx = isnan(data);

    % If missing data is at the beginning or end, remove those elements
    startIdx = find(~missingIdx, 1, 'first');
    endIdx = find(~missingIdx, 1, 'last');

    if isempty(startIdx) || isempty(endIdx)
        % All data is missing, return empty results
        filteredData = [];
        originalIndices = [];
        return;
    end

    % Remove data before startIdx and after endIdx
    data = data(startIdx:endIdx);
    originalIndices = startIdx:endIdx;

    % Update missingIdx for the trimmed data
    missingIdx = isnan(data);

    % Perform linear interpolation for missing values in the middle
    x = 1:length(data);
    data(missingIdx) = interp1(x(~missingIdx), data(~missingIdx), x(missingIdx));

    % Return the filtered data and the final indices
    filteredData = data;
end

function [filteredData1, filteredData2, commonIndices] = fillAndFilterDataTwoArrays(data1, data2)
    % Process both arrays
    [filteredData1, originalIndices1] = processArray(data1);
    [filteredData2, originalIndices2] = processArray(data2);

    % Find common indices
    commonIndices = intersect(originalIndices1, originalIndices2);

    % Map filtered data to common indices
    filteredData1 = filteredData1(ismember(originalIndices1, commonIndices));
    filteredData2 = filteredData2(ismember(originalIndices2, commonIndices));
end


%% calculate wavelet power spectrum (energy density)
dt = 0.1;
npoints = 18000;
lenscales = floor(log2(npoints));
s0 = 2*dt; % following Torrence et al. 1998
scales = s0 .* 2.^(-3:lenscales);
E_f_T_m = nan(length(scales), 3, 3);
E_f_w_m = nan(length(scales), 3, 3);
cospectrum_m = nan(length(scales), 3, 3);
[E_f_T_m(:, :, 3), E_f_w_m(:, :, 3), ~, cospectrum_m(:, :, 3), ~] = get_power(ww3, tt3, 5, scales, CumSec3, CumSec_n);
[E_f_T_m(:, :, 2), E_f_w_m(:, :, 2), ~, cospectrum_m(:, :, 2), ~] = get_power(ww10, tt10, 4, scales, CumSec10, CumSec_n);
[E_f_T_m(:, :, 1), E_f_w_m(:, :, 1), frequencies, cospectrum_m(:, :, 1), cospectrum_freq] = get_power(ww20, tt20, 5, scales, CumSec20, CumSec_n);

function data_interpolated = interpolate(data, timei, timen)
    if length(timei)~=length(timen) 
        if length(timei)>length(timen)
            data_interpolated = data(1:length(timen));
        else
            data_interpolated = interp1(timei, data, timen, 'linear', 'extrap');
        end
    else
        data_interpolated = data;
    end
end

function [E_f_T_m, E_f_w_m, frequencies, cospectrum_m, cospectrum_freq] = get_power(ww, tt, lentowers, scales, Cumseci, Cumsecn)
    Cg = pi;
    fc = 0.251;
    E_f_T = nan(length(scales), 3, lentowers);
    E_f_w = nan(length(scales), 3, lentowers);
    for ti=1:lentowers
        for di = 1:3
            df = interpolate(tt{di, ti}, Cumseci{di, ti}, Cumsecn{di, ti});
            [df, ~] = processArray(df);
            cwtstruct_temp = cwtft(df, 'wavelet', 'mexh', 'scales', scales);
            W_ab_T = cwtstruct_temp.cfs';
            frequencies = cwtstruct_temp.frequencies';
            var_df = mean((df-mean(df)).^2);
            E_f_T(:, di, ti) = squeeze(mean(abs(W_ab_T).^2, 1)) ./ (Cg.*fc.*var_df); % frequency-dependent wavelet energy 
            
            df = interpolate(ww{di, ti}, Cumseci{di, ti}, Cumsecn{di, ti});
            [df, ~] = processArray(df);
            cwtstruct_temp = cwtft(df, 'wavelet', 'mexh', 'scales', scales);
            W_ab_w = cwtstruct_temp.cfs';
            var_df = mean((df-mean(df)).^2);
            E_f_w(:, di, ti) = squeeze(mean(abs(W_ab_w).^2, 1)) ./ (Cg.*fc.*var_df); % frequency-dependent wavelet energy
        end
    end
    E_f_T_m = mean(E_f_T, 3);
    E_f_w_m = mean(E_f_w, 3);
    
    cospectrum = nan(length(scales), 3, lentowers);
    for ti=1:lentowers
        for di=1:3
            df_tt = interpolate(tt{di, ti}, Cumseci{di, ti}, Cumsecn{di, ti});
            [df_tt, ~] = processArray(df_tt);
            df_ww = interpolate(ww{di, ti}, Cumseci{di, ti}, Cumsecn{di, ti});
            [df_ww, ~] = processArray(df_ww);
            [cospectrum(:, di, ti), cospectrum_freq] = cross_wavelet(df_tt, df_ww, scales);
            
        end
    end
    cospectrum_m = mean(cospectrum, 3);
end


%%
% Set font properties
set(0, 'DefaultAxesFontWeight', 'bold');
set(0, 'DefaultAxesFontSize', 20);

vars = {E_f_T_m, E_f_w_m, cospectrum_m};
ylabels = {"f\timesE_{T'}(f) / T'^2 ", "f\timesE_{w'}(f) / w'^2 ", "f\timesCO_{w'T'}(f) / w'T' "};
heights = {'19 m', '10 m', '3 m'};
colors = ['b', 'r', 'k'];
figure('Position', [200 100 1200 800])
for hi=1:3
    for vi=1:3
        subaxis(3, 3, (hi-1)*3+vi, 'sh', 0.032, 'sv', 0.02, 'padding', 0, 'ML', 0.08, 'MB', 0.08, 'MR', 0.07, 'MT', 0.05);
        for di=1:3
            if vi<3
                loglog(frequencies, frequencies.*vars{vi}(:, di, hi), 'marker', 'o', 'Color', colors(di), 'LineWidth',2)
            else
                loglog(cospectrum_freq, cospectrum_freq.*vars{vi}(:, di, hi), 'marker', 'o', 'Color', colors(di), 'LineWidth',2)
            end
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
        if hi==1 && vi==2
            legend('Pre-FFP', 'FFP', 'Post-FFP', '-2/3', '-4/3', 'box', 'off', 'Orientation','horizontal')
        end
        xlim([0.0001 2])
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
            text(xL(2)*6, yL(2)/100, heights{hi},'HorizontalAlignment','right','VerticalAlignment','middle', "FontSize",24,"FontWeight","bold")
        end

    end
end

%%
tt = {tt3, tt10, tt20};
ww = {ww3, ww10, ww20};
Cumsec = {CumSec3, CumSec10, CumSec20};
towerstitle = {'Flux', 'North', 'West', 'East', 'South'};
figure('Position', [200 100 1400 800])

for ti=1:5
    for hi=1:3
        if ~(hi==2 && ti==5)
            subaxis(5, 3, (ti-1)*3+hi, 'sh', 0.003, 'sv', 0.028, 'padding', 0, 'ML', 0.05, 'MB', 0.04, 'MR', 0.05, 'MT', 0.04);
            
            data1 = [processArray(interpolate(ww{hi}{1,ti}, Cumsec{hi}{1, ti}, CumSec_n{1, ti}));...
                processArray(interpolate(ww{hi}{2,ti}, Cumsec{hi}{2, ti}, CumSec_n{2, ti}));...
                processArray(interpolate(ww{hi}{3,ti}, Cumsec{hi}{3, ti}, CumSec_n{3, ti}))];
            data2 = [processArray(interpolate(tt{hi}{1,ti}, Cumsec{hi}{1, ti}, CumSec_n{1, ti}));...
                processArray(interpolate(tt{hi}{2,ti}, Cumsec{hi}{2, ti}, CumSec_n{2, ti}));...
                processArray(interpolate(tt{hi}{3,ti}, Cumsec{hi}{3, ti}, CumSec_n{3, ti}))];
            [filteredData1, filteredData2, commonIndices] = fillAndFilterDataTwoArrays(data1, data2);
            wcoherence(filteredData1, filteredData2,seconds(0.1),PhaseDisplayThreshold=0.6);
            hold on
            xline(30*60-1, 'Color','r','LineWidth',2,'LineStyle','--')
            hold on
            if ti==2 || ti==3 || ti==5
                xline((30+20)*60-1, 'Color','r','LineWidth',2,'LineStyle','--')
                hold on
            elseif ti==1 || ti==4
                xline((30+15)*60-1, 'Color','r','LineWidth',2,'LineStyle','--')
                hold on
            end
            if ti==4 && hi==2
                cb = colorbar;
                cb.Location = "southoutside";
                cb.Position = [0.36 0.18 0.28 0.01]; % [left, bottom, width, height]
            else
                colorbar('off'); % Turn off colorbar visibility for this subplot
            end
            xlabel('')
            xL=xlim;
            yL=ylim;
            if ti==1
                if hi<3
                    xticks([0, 30*60-1, (30+15)*60-1]);
                    xticklabels({'14:55', '15:25', '15:40'});
                else
                    xticks([0, 30*60-1, (30+15)*60-1, xL(2)]);
                    xticklabels({'14:55', '15:25', '15:40', '16:10'});
                end
            elseif ti==2
                if hi<3
                    xticks([0, 30*60-1, (30+20)*60-1]);
                    xticklabels({'15:43', '16:13', '16:33'});
                else
                    xticks([0, 30*60-1, (30+20)*60-1, (30+20+30)*60-1]);
                    xticklabels({'15:43', '16:13', '16:33', '17:03'});
                end
            elseif ti==3
                if hi<3
                    xticks([0, 30*60-1, (30+20)*60-1]);
                    xticklabels({'14:55', '15:25', '15:45'});
                else
                    xticks([0, 30*60-1, (30+20)*60-1, (30+20+30)*60-1]);
                    xticklabels({'14:55', '15:25', '15:45', '16:15'});
                end
            elseif ti==4
                if hi<3
                    xticks([0, 30*60-1, (30+15)*60-1]);
                    xticklabels({'15:08', '15:38', '15:53'});
                else
                    xticks([0, 30*60-1, (30+15)*60-1, (30+15+30)*60-1]);
                    xticklabels({'15:08', '15:38', '15:53', '16:23'});
                end
            elseif ti==5
                xticks([0, 30*60-1, (30+20)*60-1, (30+20+30)*60-1]);
                xticklabels({'14:25', '14:55', '15:15', '15:45'});
            end
            
            if hi>1
                ylabel('')
                yticklabels('')
            else
                yticks(0:2:11);
                yticklabels(2.^(0:2:11));
                ylabel('Period [s]')
            end
            
            if ti>1
                title('')
            else
                title(heights{3-hi+1}, "FontSize",24)
            end
            
            if hi==3
                text(xL(2)*1.008, yL(2)/2, towerstitle{ti},'HorizontalAlignment','left','VerticalAlignment','middle', "FontSize",22,"FontWeight","bold")
            end
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
