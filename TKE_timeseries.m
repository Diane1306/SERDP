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
buoy = cell(size(CumSec10));
shear = cell(size(CumSec10));
transport = cell(size(CumSec10));
dissipation_rate = cell(size(CumSec10));
time_ave = cell(size(CumSec10));

for di=1:3
    for ti=1:4
        buoy{di, ti} = calc_1min_buoy(T10{di,ti}, wt10{di,ti}, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
        shear{di, ti} = calc_1min_shear(uw10{di,ti}, vw10{di,ti}, ww10{di,ti}, U3{di,ti}, U10{di,ti}, U20{di,ti}, V3{di,ti}, V10{di,ti}, V20{di,ti}, W3{di,ti}, W10{di,ti}, W20{di,ti}, height, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
        transport{di, ti} = calc_1min_transport(w3{di,ti}, w10{di,ti}, w20{di,ti}, tke3{di,ti}, tke10{di,ti}, tke20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec10{di, ti}, CumSec20{di, ti}, CumSec_n{di, ti});
        dissipation_rate{di, ti} = calc_1min_disp(usw10{di,ti}, Usw10{di,ti}, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
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

function buoy = calc_1min_buoy(Tm2, wt2, samplingfreq, time2, timen)
    if length(time2)~=length(timen)
        if length(time2)>length(timen)
            wt2 = wt2(1:length(timen));
        else
            wt2 = interp1(time2, wt2, timen, 'linear', 'extrap');
        end
    end
    g = 9.81;  %m/s^2 gravity
    buoy_temp = g ./ (Tm2 + 273.15) .* wt2;

    buoy = calc_1min_ave(buoy_temp, samplingfreq);
end

function shear = calc_1min_shear(uw2, vw2, ww2, u1, u2, u3, v1, v2, v3, w1, w2, w3, height, samplingfreq, time2, timen)
    if length(time2)~=length(timen) 
        if length(time2)>length(timen)
            uw2 = uw2(1:length(timen));
            vw2 = vw2(1:length(timen));
            ww2 = ww2(1:length(timen));
        else
            uw2 = interp1(time2, uw2, timen, 'linear', 'extrap');
            vw2 = interp1(time2, vw2, timen, 'linear', 'extrap');
            ww2 = interp1(time2, ww2, timen, 'linear', 'extrap');
        end
    end
    % Weighted Finite Difference Method
    weight1 = (height(3)-height(2)) / (height(3)-height(1));
    weight3 = (height(2)-height(1)) / (height(3)-height(1));
    grad_U = weight1 .* (u2-u1)./(height(2)-height(1)) + weight3 .* (u3-u2)./(height(3)-height(2));
    grad_V = weight1 .* (v2-v1)./(height(2)-height(1)) + weight3 .* (v3-v2)./(height(3)-height(2));
    grad_W = weight1 .* (w2-w1)./(height(2)-height(1)) + weight3 .* (w3-w2)./(height(3)-height(2));

    shear_temp = -uw2 .* grad_U -vw2 .* grad_V -ww2 .* grad_W;
    shear = calc_1min_ave(shear_temp, samplingfreq);
end

function transport = calc_1min_transport(w1, w2, w3, tke1, tke2, tke3, height, samplingfreq, time1, time2, time3, timen)
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
    % Weighted Finite Difference Method
    weight1 = (height(3)-height(2)) / (height(3)-height(1));
    weight3 = (height(2)-height(1)) / (height(3)-height(1));
    transport_temp = weight1 .* (w2.*tke2 - w1.*tke1)./(height(2)-height(1)) + weight3 .* (w3.*tke3 - w2.*tke2)./(height(3)-height(2));
    transport = calc_1min_ave(transport_temp, samplingfreq);
end

function dissipation_rate = calc_1min_disp(u2, Um2, samplingfreq, time2, timen)
    if length(time2)~=length(timen)
        if length(time2)>length(timen)
            u2 = u2(1:length(timen));
        else
            u2 = interp1(time2, u2, timen, 'linear', 'extrap');
        end
    end
    structure_function = calc_structure_function(u2);
    epsilon = calc_dissipation_rate(structure_function, Um2, length(u2), 1/samplingfreq);
    dissipation_rate = -epsilon;
end

function epsilon = calc_dissipation_rate(structure_function, U, nPoints, dt)
    % U is the mean velocity to convert time to distance
    y1 = structure_function;
    timemax_insec = nPoints*dt/2;
    r_value = linspace(dt.*U, timemax_insec.*U, floor(nPoints/2));
    y2 = r_value.^(2/3);
    [~, I] = min(abs(log(y2(1:20))-log(y1(1:20))));
    c_A2 = (y2 ./ (y2(I) ./ (y1(I)))) ./ y2;
    epsilon_array = (c_A2 .^ (3/2)) .* 0.35;
    epsilon = epsilon_array(1);
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

%%
xvalues = 1:80;
% Set font properties
set(0, 'DefaultAxesFontWeight', 'bold');
set(0, 'DefaultAxesFontSize', 20);
figure('Position', [200 100 1200 800])
for ti=1:4
    subaxis(4, 1, ti, 'sh', 0, 'sv', 0.01, 'padding', 0, 'ML', 0.08, 'MB', 0.08, 'MR', 0.07, 'MT', 0.05);
    plot([buoy{1,ti}; buoy{2,ti}; buoy{3,ti}], 'r', 'LineWidth',2)
    hold on
    plot([shear{1,ti}; shear{2,ti}; shear{3,ti}], 'b', 'LineWidth',2)
    hold on
    plot([transport{1,ti}; transport{2,ti}; transport{3,ti}], 'c', 'LineWidth',2)
    hold on
    plot([dissipation_rate{1,ti}.*ones(size(buoy{1,ti})); dissipation_rate{2,ti}.*ones(size(buoy{2,ti})); dissipation_rate{3,ti}.*ones(size(buoy{3,ti}))], 'k--', 'LineWidth',2)
    hold on
    ylim([-0.3 0.6])
    xlim([0 80])
    if ti~=4
        xticklabels('')
    end
    grid on
    xL=xlim;
    yL=ylim;
    text(xL(2), yL(2), towers{ti},'HorizontalAlignment','right','VerticalAlignment','top', "FontSize",20,"FontWeight","bold")
    if ti==1 || ti==4
        xline(31, 'm--', 'LineWidth',2)
        hold on
        xline(45, 'm--', 'LineWidth',2)
    else
        xline(31, 'm--', 'LineWidth',2)
        hold on
        xline(50, 'm--', 'LineWidth',2)
    end
    if ti==1
        legend('buoyancy', 'shear', 'transport', 'dissipation', 'box', 'off', 'Orientation','horizontal')
    end
    if ti==2
        ylabel('TKE Budget Components [m^2 s^{-3}]')
    end
    if ti==4
        xlabel('Minutes')
    end
end