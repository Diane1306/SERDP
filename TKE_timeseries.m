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
%% dissipation rate
height = [3, 10, 19];
di=1;
ti=1;
samplingfreq = 10;
buoy = calc_1min_buoy(T10{di,ti}, wt10{di,ti}, samplingfreq, ~isnan(CumSec10{di,ti}));
shear = calc_1min_shear(uw10{di,ti}, vw10{di,ti}, ww10{di,ti}, U3{di,ti}, U10{di,ti}, U20{di,ti}, V3{di,ti}, V10{di,ti}, V20{di,ti}, W3{di,ti}, W10{di,ti}, W20{di,ti}, height, samplingfreq, ~isnan(CumSec10{di,ti}));
transport = calc_1min_transport(W3{di,ti}, W10{di,ti}, W20{di,ti}, tke3{di,ti}, tke10{di,ti}, tke20{di,ti}, height, samplingfreq, ~isnan(CumSec10{di,ti}));
dissipation_rate = calc_1min_disp(usw10{di,ti}, Usw10{di,ti}, samplingfreq, ~isnan(CumSec10{di,ti}));

function term = calc_1min_ave(term_temp, samplingfreq)
    len_term = ceil(length(term_temp)/(samplingfreq*60));
    term = zeros(len_term, 1);
    for i=1:len_term-1
        term(i) = mean(term_temp((i-1)*samplingfreq*60+1:i*samplingfreq*60), "omitmissing");
    end
    term(len_term) = mean(term_temp((len_term-1)*samplingfreq*60+1:end), "omitmissing");
end

function buoy = calc_1min_buoy(Tm, wt, samplingfreq, flag)
    % flag is indicated by ~isnan(CumSec)
    wt = wt(flag);
    g = 9.81;  %m/s^2 gravity
    buoy_temp = g ./ Tm .* wt;

    buoy = calc_1min_ave(buoy_temp, samplingfreq);
end

function shear = calc_1min_shear(uw, vw, ww, u1, u2, u3, v1, v2, v3, w1, w2, w3, height, samplingfreq, flag)
    uw = uw(flag);
    vw = vw(flag);
    ww = ww(flag);
    
    % Weighted Finite Difference Method
    weight1 = (height(3)-height(2)) / (height(3)-height(1));
    weight3 = (height(2)-height(1)) / (height(3)-height(1));
    grad_U = weight1 .* (u2-u1)./(height(2)-height(1)) + weight3 .* (u3-u2)./(height(3)-height(2));
    grad_V = weight1 .* (v2-v1)./(height(2)-height(1)) + weight3 .* (v3-v2)./(height(3)-height(2));
    grad_W = weight1 .* (w2-w1)./(height(2)-height(1)) + weight3 .* (w3-w2)./(height(3)-height(2));

    shear_temp = -uw .* grad_U -vw .* grad_V -ww .* grad_W;
    shear = calc_1min_ave(shear_temp, samplingfreq);
end

function transport = calc_1min_transport(w1, w2, w3, tke1, tke2, tke3, height, samplingfreq, flag)
    w1 = w1(flag);
    w2 = w2(flag);
    w3 = w3(flag);
    tke1 = tke1(flag);
    tke2 = tke2(flag);
    tke3 = tke3(flag);

    % Weighted Finite Difference Method
    weight1 = (height(3)-height(2)) / (height(3)-height(1));
    weight3 = (height(2)-height(1)) / (height(3)-height(1));
    transport_temp = weight1 .* (w2.*tke2 - w1.*tke1)./(height(2)-height(1)) + weight3 .* (w3.*tke3 - w2.*tke2)./(height(3)-height(2));
    transport = calc_1min_ave(transport_temp, samplingfreq);
end

function dissipation_rate = calc_1min_disp(u, Um, samplingfreq, flag)
    u = u(flag);
    structure_function = calc_structure_function(u);
    epsilon = calc_dissipation_rate(structure_function, Um, length(u), 1/samplingfreq);

    dissipation_rate = calc_1min_ave(epsilon, samplingfreq);
end

function epsilon = calc_dissipation_rate(structure_function, U, nPoints, dt)
    % U is the mean velocity to convert time to distance
    y1 = structure_function';
    timemax_insec = nPoints*dt/2;
    r_value = linspace(dt.*U, timemax_insec.*U, floor(nPoints/2));
    y2 = r_value.^(2/3);
    [~, I] = min(abs(log(y2(1:20))-log(y1(1:20))));
    c_A2 = (y2 ./ (y2(I) / (y1(I)))) ./ y2;
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
