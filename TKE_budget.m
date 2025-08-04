clc
clear
%%
% Main script
data_dir = '/Users/diane_wt/Downloads/work/PreMean/'; % 2019_East_Tower-fft_20mAGL
towers = {'Flux', 'North', 'West', 'East', 'South_Mobile'};

% Call the function
[T20, U20, V20, W20, wt20, uw20, vw20, ww20, w20, tke20, usw20, Usw20, CumSec20] = get_data(data_dir, towers, '20m');
[T10, U10, V10, W10, wt10, uw10, vw10, ww10, w10, tke10, usw10, Usw10, CumSec10] = get_data(data_dir, towers(1:4), '10m');
[T3, U3, V3, W3, wt3, uw3, vw3, ww3, w3, tke3, usw3, Usw3, CumSec3] = get_data(data_dir, towers, '3m');

% Function Definition (for example, if not already implemented)
function [T, U, V, W, wt, uw, vw, ww, w, tke, usw, Usw, CumSec] = get_data(dr, towers, height)
    lentowers = length(towers);
    lenperiods = 3;
    % Initialize output cell arrays
    T = cell(lenperiods, lentowers);
    U = cell(lenperiods, lentowers);
    V = cell(lenperiods, lentowers);
    W = cell(lenperiods, lentowers);
    wt = cell(lenperiods, lentowers);
    uw = cell(lenperiods, lentowers);
    vw = cell(lenperiods, lentowers);
    ww = cell(lenperiods, lentowers);
    w = cell(lenperiods, lentowers);
    tke = cell(lenperiods, lentowers);
    usw = cell(lenperiods, lentowers);
    Usw = cell(lenperiods, lentowers);
    CumSec = cell(lenperiods, lentowers);
    
    for ti=1:length(towers)
        % Read the sheets
        df = cell(3,1);
        df{1} = readtable(strcat(dr,height,'/2019_',towers{ti},'_Tower-fft_',height,'AGL.xlsx'), 'Sheet', 'Pre-FFP', 'VariableNamingRule', 'preserve');
        df{2} = readtable(strcat(dr,height,'/2019_',towers{ti},'_Tower-fft_',height,'AGL.xlsx'), 'Sheet', 'FFP', 'VariableNamingRule', 'preserve');
        df{3} = readtable(strcat(dr,height,'/2019_',towers{ti},'_Tower-fft_',height,'AGL.xlsx'), 'Sheet', 'Post-FFP', 'VariableNamingRule', 'preserve');     
        % Loop through each phase (Pre, During, Post)
        for di = 1:3
            % Extract and store the relevant columns
            T{di, ti} = mean(df{di}.("T (C)"), 'omitmissing');
            U{di, ti} = mean(df{di}.("U (m/s)"), 'omitmissing');
            V{di, ti} = mean(df{di}.("V (m/s)"), 'omitmissing');
            W{di, ti} = mean(df{di}.("W (m/s)"), 'omitmissing');
            wt{di, ti} = df{di}.("w't' (mC/s)");
            uw{di, ti} = df{di}.("u'w' (m2/s2)");
            vw{di, ti} = df{di}.("v'w' (m2/s2)");
            ww{di, ti} = df{di}.("w'2 (m2/s2)");
            w{di, ti} = df{di}.("w' (m/s)");
            tke{di, ti} = df{di}.("TKE (m2/s2)");
            usw{di, ti} = df{di}.("S' (m/s)");
            Usw{di, ti} = mean(df{di}.("Horizontal (Streamwise) Velocity (m/s)"), 'omitmissing');
            CumSec{di, ti} = df{di}.("Cum. Sec.");  
        end
        clear df
    end
end

%% create standard CumSec
CumSec_n = CumSec20;
CumSec_n{3, 5} = CumSec3{3, 5};
time_temp = CumSec20{2, 5}(1):0.1:CumSec20{2, 5}(1)+(12000-1)*0.1;
CumSec_n{2, 5} = time_temp';

%% dissipation rate
height = [3, 10, 19];
samplingfreq = 10;
buoy = cell([size(CumSec10),3]);
shear = cell([size(CumSec10),3]);
transport = cell([size(CumSec10),3]);
dissipation_rate = cell([size(CumSec10),3]);
time_ave = cell(size(CumSec10));

for di=1:3
    for ti=1:4
        buoy{di, ti, 1} = calc_1min_buoy(T3{di,ti}, wt3{di,ti}, samplingfreq, CumSec3{di, ti}, CumSec_n{di, ti});
        shear{di, ti, 1} = calc_1min_shear("forward", uw3{di,ti}, vw3{di,ti}, ww3{di,ti}, U3{di,ti}, U10{di,ti}, U20{di,ti}, V3{di,ti}, V10{di,ti}, V20{di,ti}, W3{di,ti}, W10{di,ti}, W20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec_n{di, ti});
        transport{di, ti, 1} = calc_1min_transport("forward", w3{di,ti}, w10{di,ti}, w20{di,ti}, tke3{di,ti}, tke10{di,ti}, tke20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec10{di, ti}, CumSec20{di, ti}, CumSec_n{di, ti});
        dissipation_rate{di, ti, 1} = calc_1min_disp(usw3{di,ti}, Usw3{di,ti}, samplingfreq, CumSec3{di, ti}, CumSec_n{di, ti});
        
        buoy{di, ti, 2} = calc_1min_buoy(T10{di,ti}, wt10{di,ti}, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
        shear{di, ti, 2} = calc_1min_shear("center", uw10{di,ti}, vw10{di,ti}, ww10{di,ti}, U3{di,ti}, U10{di,ti}, U20{di,ti}, V3{di,ti}, V10{di,ti}, V20{di,ti}, W3{di,ti}, W10{di,ti}, W20{di,ti}, height, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
        transport{di, ti, 2} = calc_1min_transport("center", w3{di,ti}, w10{di,ti}, w20{di,ti}, tke3{di,ti}, tke10{di,ti}, tke20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec10{di, ti}, CumSec20{di, ti}, CumSec_n{di, ti});
        dissipation_rate{di, ti, 2} = calc_1min_disp(usw10{di,ti}, Usw10{di,ti}, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
        
        buoy{di, ti, 3} = calc_1min_buoy(T20{di,ti}, wt20{di,ti}, samplingfreq, CumSec20{di, ti}, CumSec_n{di, ti});
        shear{di, ti, 3} = calc_1min_shear("backward", uw20{di,ti}, vw20{di,ti}, ww20{di,ti}, U3{di,ti}, U10{di,ti}, U20{di,ti}, V3{di,ti}, V10{di,ti}, V20{di,ti}, W3{di,ti}, W10{di,ti}, W20{di,ti}, height, samplingfreq, CumSec20{di, ti}, CumSec_n{di, ti});
        transport{di, ti, 3} = calc_1min_transport("backward", w3{di,ti}, w10{di,ti}, w20{di,ti}, tke3{di,ti}, tke10{di,ti}, tke20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec10{di, ti}, CumSec20{di, ti}, CumSec_n{di, ti});
        dissipation_rate{di, ti, 3} = calc_1min_disp(usw20{di,ti}, Usw20{di,ti}, samplingfreq, CumSec20{di, ti}, CumSec_n{di, ti});
        
        time_ave{di, ti} = calc_1min_ave(CumSec_n{di, ti}, samplingfreq);
    end
end

function term = calc_1min_ave(term_temp, samplingfreq)
    len_term = ceil(length(term_temp)/(samplingfreq*60));
    term = zeros(len_term, 1);
    for i=1:len_term-1
        term(i) = mean(term_temp((i-1)*samplingfreq*60+1:i*samplingfreq*60), "omitmissing");
    end
    term(len_term) = mean(term_temp((len_term-1)*samplingfreq*60+1:end), "omitmissing");
end

function buoy = calc_1min_buoy(Tmi, wti, samplingfreq, timei, timen)
    if length(timei)~=length(timen)
        if length(timei)>length(timen)
            wti = wti(1:length(timen));
        else
            wti = interp1(timei, wti, timen, 'linear', 'extrap');
        end
    end
    g = 9.81;  %m/s^2 gravity
    buoy_temp = g ./ (Tmi + 273.15) .* wti;

    buoy = calc_1min_ave(buoy_temp, samplingfreq);
end

function shear = calc_1min_shear(grad_method, uwi, vwi, wwi, u1, u2, u3, v1, v2, v3, w1, w2, w3, height, samplingfreq, timei, timen)
    if length(timei)~=length(timen) 
        if length(timei)>length(timen)
            uwi = uwi(1:length(timen));
            vwi = vwi(1:length(timen));
            wwi = wwi(1:length(timen));
        else
            uwi = interp1(timei, uwi, timen, 'linear', 'extrap');
            vwi = interp1(timei, vwi, timen, 'linear', 'extrap');
            wwi = interp1(timei, wwi, timen, 'linear', 'extrap');
        end
    end

    if strcmp(grad_method, "center") % 10m
        % Weighted Finite Difference Method
        weight1 = (height(3)-height(2)) / (height(3)-height(1));
        weight3 = (height(2)-height(1)) / (height(3)-height(1));
        grad_U = weight1 .* (u2-u1)./(height(2)-height(1)) + weight3 .* (u3-u2)./(height(3)-height(2));
        grad_V = weight1 .* (v2-v1)./(height(2)-height(1)) + weight3 .* (v3-v2)./(height(3)-height(2));
        grad_W = weight1 .* (w2-w1)./(height(2)-height(1)) + weight3 .* (w3-w2)./(height(3)-height(2));
    elseif strcmp(grad_method, "forward") % 3m
        grad_U = (u2-u1)./(height(2)-height(1));
        grad_V = (v2-v1)./(height(2)-height(1));
        grad_W = (w2-w1)./(height(2)-height(1));
    elseif strcmp(grad_method, "backward") % 19m
        grad_U = (u3-u2)./(height(3)-height(2));
        grad_V = (v3-v2)./(height(3)-height(2));
        grad_W = (w3-w2)./(height(3)-height(2));
    end

    shear_temp = -uwi .* grad_U -vwi .* grad_V -wwi .* grad_W;
    shear = calc_1min_ave(shear_temp, samplingfreq);
end

function transport = calc_1min_transport(grad_method, w1, w2, w3, tke1, tke2, tke3, height, samplingfreq, time1, time2, time3, timen)
    if length(time1)~=length(timen)
        if length(time1)>length(timen)
            w1 = w1(1:length(timen));
            tke1 = tke1(1:length(timen));
        else
            w1 = interp1(time1, w1, timen, 'linear', 'extrap');
            tke1 = interp1(time1, tke1, timen, 'linear', 'extrap');
        end
    end
    if length(time2)~=length(timen)
        if length(time2)>length(timen)
            w2 = w2(1:length(timen));
            tke2 = tke2(1:length(timen));
        else
            w2 = interp1(time2, w2, timen, 'linear', 'extrap');
            tke2 = interp1(time2, tke2, timen, 'linear', 'extrap');
        end
    end
    if length(time3)~=length(timen)
        if length(time3)>length(timen)
            w3 = w3(1:length(timen));
            tke3 = tke3(1:length(timen));
        else
            w3 = interp1(time3, w3, timen, 'linear', 'extrap');
            tke3 = interp1(time3, tke3, timen, 'linear', 'extrap');
        end
    end

    if strcmp(grad_method, "center") % 10m
        % Weighted Finite Difference Method
        weight1 = (height(3)-height(2)) / (height(3)-height(1));
        weight3 = (height(2)-height(1)) / (height(3)-height(1));
        transport_temp = weight1 .* (w2.*tke2 - w1.*tke1)./(height(2)-height(1)) + weight3 .* (w3.*tke3 - w2.*tke2)./(height(3)-height(2));
    elseif strcmp(grad_method, "forward")
        transport_temp = (w2.*tke2 - w1.*tke1)./(height(2)-height(1));
    elseif strcmp(grad_method, "backward")
        transport_temp = (w3.*tke3 - w2.*tke2)./(height(3)-height(2));
    end
    transport = calc_1min_ave(transport_temp, samplingfreq);
end

function dissipation_rate = calc_1min_disp(ui, Umi, samplingfreq, timei, timen)
    if length(timei)~=length(timen)
        if length(timei)>length(timen)
            ui = ui(1:length(timen));
        else
            ui = interp1(timei, ui, timen, 'linear', 'extrap');
        end
    end
    structure_function = calc_structure_function(ui);
    epsilon = calc_dissipation_rate(structure_function, Umi, length(ui), 1/samplingfreq, 'structure');
    % [~, spectra_function, ~] = calc_spectrum(ui, samplingfreq, 2);
    % epsilon = calc_dissipation_rate(spectra_function, Umi, length(ui), 1/samplingfreq, 'spectra');
    dissipation_rate = -epsilon;
end

function epsilon = calc_dissipation_rate(structure_function, U, nPoints, dt, method)
    % U is the mean velocity to convert time to distance
    if strcmp(method, 'structure')
        y1 = structure_function;
        timemax_insec = nPoints*dt/2;
        r_value = linspace(dt.*U, timemax_insec.*U, floor(nPoints/2));
        y2 = r_value.^(2/3);
        [~, I] = min(abs(log(y2(1:20))-log(y1(1:20))));
        c_A2 = (y2 ./ (y2(I) ./ (y1(I)))) ./ y2;
        epsilon_array = (c_A2 .^ (3/2)) .* 0.35;
        epsilon = epsilon_array(1);
    elseif strcmp(method, 'spectra')
        y1 = structure_function';
        n = 2 .* length(structure_function);  % length of dataset
        Nf = floor(n / 2);  % Nyquist frequency
        freq = fftshift(fftfreq(n, 1 / 20));  % frequency corresponding to sampling rate of fs, centered at 0
        freq = freq(Nf+1:end);  % only positive frequencies

        [x1, y1] = smooth2E_AddHighFreqPoints(freq, y1, 2048);
        y2 = x1 .^ (-5/3);
        istart=7;
        [~, I] = min(abs(log(y2(istart:end-3))-log(y1(istart:end-3))));
        epsilon = 0.49 .* ((y2(istart)./ (y2(istart-1+I) ./ y1(istart-1+I))).^(3/2)) .* (x1(istart).^(5/2));
    end
end

function structure_function = calc_structure_function(u)
    % u is the velocity time series
    n = length(u);
    r = 1:floor(n/2);
    structure_function = zeros(size(r));
    for i = 1:length(r)
        ri = r(i);
        vel_diff = (u(1+ri:end) - u(1:end-ri)).^2;
        structure_function(i) = mean(vel_diff(:), 'omitmissing'); % average over the grid
    end
end

function [freq, spectrum, sum_spectra] = calc_spectrum(data, fs, return_spectra)
    % Calculate various spectra from input data
    % Input:
    %   data: input data vector
    %   fs: sampling frequency
    %   return_spectra: 1 for energy spectrum, 2 for spectral density, 3 for frequency weighted spectral density
    % Output:
    %   freq: frequency vector
    %   spectrum: spectrum based on return_spectra choice
    
    n = length(data);  % length of dataset
    Nf = floor(n / 2);  % Nyquist frequency
    freq = fftshift(fftfreq(n, 1 / fs));  % frequency corresponding to sampling rate of fs, centered at 0
    freq = freq(Nf+1:end);  % only positive frequencies
    FFT = fft(data) / n;  % complex Fourier transform
    E_spectra = zeros(1, Nf);  % energy spectrum
    
    if mod(n, 2)
        E_spectra = 2 * (abs(FFT(2:Nf+1)) .^ 2);  % for odd n
    else
        E_spectra(1:Nf-1) = 2 * (abs(FFT(2:Nf)) .^ 2);  % for even n
        E_spectra(Nf) = abs(FFT(Nf+1)) .^ 2;
    end
    
    sum_spectra = sum(E_spectra);

    delta_f = 1 / n;
    S_spectra = E_spectra / delta_f;  % convert energy spectrum to spectral density
    fw_S_spectra = freq .* S_spectra;  % frequency weighted spectral density

    if return_spectra == 1
        spectrum = E_spectra;
    elseif return_spectra == 2
        spectrum = S_spectra;
    elseif return_spectra == 3
        spectrum = fw_S_spectra;
    end
end

function [fs, ps, binindices] = smooth2E_AddHighFreqPoints(f, p, limit)
    if nargin < 3
        limit = 2048;
    end
    
    l = length(f);
    binindices = [1];
    
    if l < limit
        lr = l;
    else
        lr = l - limit;
    end
    i = 1;
    while binindices(i) < lr
        if 2^i <= limit
            binindices = [binindices, 2^i];
        else
            binindices = [binindices, binindices(i) + limit];
        end
        i = i + 1;
    end
    
    if binindices(i) < l
        binindices(i) = l;
    end

    fs = [f(1), f(2), arrayfun(@(i) mean(f(binindices(i)+1:binindices(i+1))), 2:length(binindices)-1)];
    ps = [p(1), p(2), arrayfun(@(i) mean(p(binindices(i)+1:binindices(i+1))), 2:length(binindices)-1)];
end
%%
xvalues = 1:80;
timeticks = {{' ', '15:05', ' ', '15:25', ' ', '15:45', ' ', '16:05', ' '},...
    {' ', '15:53', ' ', '16:13', ' ', '16:33', ' ', '16:53', ' '},...
    {' ', '15:05', ' ', '15:25', ' ', '15:45', ' ', '16:05', ' '},...
    {' ', '15:18', ' ', '15:38', ' ', '15:58', ' ', '16:18', ' '}};
% Set font properties
set(0, 'DefaultAxesFontWeight', 'bold');
set(0, 'DefaultAxesFontSize', 20);
figure('Position', [200 100 1200 800])
for ti=1:4
    for hi=1:3
        if ti<3
            pi = ((ti-1)*3+hi)*2-1;
        else
            pi = ((ti-1)*3+hi)*2-12;
        end
        if pi<7
            mt_value = 0.05;
        else
            mt_value = 0.08;
        end
        subaxis(6, 2, pi, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'ML', 0.06, 'MB', 0.08, 'MR', 0.01, 'MT', mt_value);
        plot([buoy{1, ti, 3-hi+1}; buoy{2,ti, 3-hi+1}; buoy{3,ti, 3-hi+1}], 'r', 'LineWidth',2)
        hold on
        plot([shear{1,ti, 3-hi+1}; shear{2,ti, 3-hi+1}; shear{3,ti, 3-hi+1}], 'b', 'LineWidth',2)
        hold on
        plot([transport{1,ti, 3-hi+1}; transport{2,ti, 3-hi+1}; transport{3,ti, 3-hi+1}], 'c', 'LineWidth',2)
        hold on
        plot([dissipation_rate{1,ti, 3-hi+1}.*ones(size(buoy{1,ti, 3-hi+1})); dissipation_rate{2,ti, 3-hi+1}.*ones(size(buoy{2,ti, 3-hi+1})); dissipation_rate{3,ti, 3-hi+1}.*ones(size(buoy{3,ti, 3-hi+1}))], 'k--', 'LineWidth',2)
        hold on
        ylim([-0.4 0.75])
        xlim([0 80])
        xticks(0:10:80)
        if ti>2
            yticklabels('')
        end
        if ~ismember(pi, [5, 6, 11, 12])
            xticklabels('')
        else
            xticklabels(timeticks{ti})
            xlabel('EDT')
        end
        grid on
        xL=xlim;
        yL=ylim;
        if hi==1
            text(xL(2), yL(2), towers{ti},'HorizontalAlignment','right','VerticalAlignment','top', "FontSize",20,"FontWeight","bold")
        end
        text(xL(1), yL(2), [num2str(height(3-hi+1)),' m'],'HorizontalAlignment','left','VerticalAlignment','top', "FontSize",20,"FontWeight","bold")
        
        if ti==1 || ti==4
            xline(30, 'm--', 'LineWidth',2)
            hold on
            xline(45, 'm--', 'LineWidth',2)
        else
            xline(30, 'm--', 'LineWidth',2)
            hold on
            xline(50, 'm--', 'LineWidth',2)
        end
        if pi==1
            legend('buoyancy', 'shear', 'transport', 'dissipation', 'box', 'off', 'Orientation','horizontal', "FontSize",20)
        end
        if pi==5
            ylabel('TKE Budget Components [m^2 s^{-3}]', "FontSize",20)
        end
    end
end