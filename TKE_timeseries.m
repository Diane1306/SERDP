clc
clear
%%
% Main script
data_dir = '/Users/diane_wt/Downloads/work/PreMean/'; % 2019_East_Tower-fft_20mAGL
towers = {'Flux', 'North', 'West', 'East', 'South_Mobile'};

% Call the function
[T20, U20, V20, W20, wt20, uw20, vw20, ww20, w20, tke20, usw20, Usw20, CumSec20] = get_data(data_dir, towers(1), '20m');
[T10, U10, V10, W10, wt10, uw10, vw10, ww10, w10, tke10, usw10, Usw10, CumSec10] = get_data(data_dir, towers(1), '10m');
[T3, U3, V3, W3, wt3, uw3, vw3, ww3, w3, tke3, usw3, Usw3, CumSec3] = get_data(data_dir, towers(1), '3m');

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
            T{di, ti} = df{di}.("Average T")(1);
            U{di, ti} = df{di}.("Average U")(1);
            V{di, ti} = df{di}.("Average V")(1);
            W{di, ti} = df{di}.("Average W")(1);
            wt{di, ti} = df{di}.("w't' (mC/s)");
            uw{di, ti} = df{di}.("u'w' (m2/s2)");
            vw{di, ti} = df{di}.("v'w' (m2/s2)");
            ww{di, ti} = df{di}.("w'2 (m2/s2)");
            w{di, ti} = df{di}.("w' (m/s)");
            tke{di, ti} = df{di}.("TKE (m2/s2)");
            usw{di, ti} = df{di}.("S' (m/s)");
            Usw{di, ti} = df{di}.("Mean Horizontal (Streamwise) Velocity (m/s)")(1);
            CumSec{di, ti} = df{di}.("Cum. Sec.");  
        end
        clear df
    end
end
%% dissipation rate

function epsilon = calc_dissipation_rate(structure_function, U, nPoints, dt)
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
    n = length(u);
    r = 1:floor(n/2);
    structure_function = zeros(size(r));
    for i = 1:length(r)
        ri = r(i);
        vel_diff = (u(1+ri:end) - u(1:end-ri)).^2;
        structure_function(i) = nanmean(vel_diff(:)); % average over the grid
    end
end
