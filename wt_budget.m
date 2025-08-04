clc
clear
%%
% Main script
data_dir = '/Users/diane_wt/Downloads/work/PreMean/'; % 2019_East_Tower-fft_20mAGL
towers = {'Flux', 'North', 'West', 'East', 'South_Mobile'};

% Call the function
[T20, ww20, t20, w20, CumSec20] = get_data(data_dir, towers, '20m');
[T10, ww10, t10, w10, CumSec10] = get_data(data_dir, towers(1:4), '10m');
[T3, ww3, t3, w3, CumSec3] = get_data(data_dir, towers, '3m');

% Function Definition (for example, if not already implemented)
function [T, ww, t, w, CumSec] = get_data(dr, towers, height)
    lentowers = length(towers);
    lenperiods = 3;
    % Initialize output cell arrays
    T = cell(lenperiods, lentowers);
    ww = cell(lenperiods, lentowers);
    t = cell(lenperiods, lentowers);
    w = cell(lenperiods, lentowers);
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
            ww{di, ti} = df{di}.("w'2 (m2/s2)");
            t{di, ti} = df{di}.("t' (C)");
            w{di, ti} = df{di}.("w' (m/s)");
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

%% budget
height = [3, 10, 19];
samplingfreq = 10;
buoy = cell([size(CumSec10),3]);
shear = cell([size(CumSec10),3]);
transport = cell([size(CumSec10),3]);
time_ave = cell(size(CumSec10));

for di=1:3
    for ti=1:4
        buoy{di, ti, 1} = calc_1min_buoy(T3{di,ti}, t3{di,ti}, samplingfreq, CumSec3{di, ti}, CumSec_n{di, ti});
        shear{di, ti, 1} = calc_1min_shear("forward", ww3{di,ti}, T3{di,ti}, T10{di,ti}, T20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec_n{di, ti});
        transport{di, ti, 1} = calc_1min_transport("forward", w3{di,ti}, w10{di,ti}, w20{di,ti}, t3{di,ti}, t10{di,ti}, t20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec10{di, ti}, CumSec20{di, ti}, CumSec_n{di, ti});
        
        buoy{di, ti, 2} = calc_1min_buoy(T10{di,ti}, t10{di,ti}, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
        shear{di, ti, 2} = calc_1min_shear("center", ww10{di,ti}, T3{di,ti}, T10{di,ti}, T20{di,ti}, height, samplingfreq, CumSec10{di, ti}, CumSec_n{di, ti});
        transport{di, ti, 2} = calc_1min_transport("center", w3{di,ti}, w10{di,ti}, w20{di,ti}, t3{di,ti}, t10{di,ti}, t20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec10{di, ti}, CumSec20{di, ti}, CumSec_n{di, ti});
        
        buoy{di, ti, 3} = calc_1min_buoy(T20{di,ti}, t20{di,ti}, samplingfreq, CumSec20{di, ti}, CumSec_n{di, ti});
        shear{di, ti, 3} = calc_1min_shear("backward", ww20{di,ti}, T3{di,ti}, T10{di,ti}, T20{di,ti}, height, samplingfreq, CumSec20{di, ti}, CumSec_n{di, ti});
        transport{di, ti, 3} = calc_1min_transport("backward", w3{di,ti}, w10{di,ti}, w20{di,ti}, t3{di,ti}, t10{di,ti}, t20{di,ti}, height, samplingfreq, CumSec3{di, ti}, CumSec10{di, ti}, CumSec20{di, ti}, CumSec_n{di, ti});
        
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

function shear = calc_1min_shear(grad_method, wwi, T1, T2, T3, height, samplingfreq, timei, timen)
    if length(timei)~=length(timen) 
        if length(timei)>length(timen)
            wwi = wwi(1:length(timen));
        else
            wwi = interp1(timei, wwi, timen, 'linear', 'extrap');
        end
    end

    if strcmp(grad_method, "center") % 10m
        % Weighted Finite Difference Method
        weight1 = (height(3)-height(2)) / (height(3)-height(1));
        weight3 = (height(2)-height(1)) / (height(3)-height(1));
        grad_T = weight1 .* (T2-T1)./(height(2)-height(1)) + weight3 .* (T3-T2)./(height(3)-height(2));
    elseif strcmp(grad_method, "forward") % 3m
        grad_T = (T2-T1)./(height(2)-height(1));
    elseif strcmp(grad_method, "backward") % 19m
        grad_T = (T3-T2)./(height(3)-height(2));
    end

    shear_temp = -wwi .* grad_T;
    shear = calc_1min_ave(shear_temp, samplingfreq);
end


function transport = calc_1min_transport(grad_method, w1, w2, w3, t1, t2, t3, height, samplingfreq, time1, time2, time3, timen)
    if length(time1)~=length(timen)
        if length(time1)>length(timen)
            w1 = w1(1:length(timen));
            t1 = t1(1:length(timen));
        else
            w1 = interp1(time1, w1, timen, 'linear', 'extrap');
            t1 = interp1(time1, t1, timen, 'linear', 'extrap');
        end
    end
    if length(time2)~=length(timen)
        if length(time2)>length(timen)
            w2 = w2(1:length(timen));
            t2 = t2(1:length(timen));
        else
            w2 = interp1(time2, w2, timen, 'linear', 'extrap');
            t2 = interp1(time2, t2, timen, 'linear', 'extrap');
        end
    end
    if length(time3)~=length(timen)
        if length(time3)>length(timen)
            w3 = w3(1:length(timen));
            t3 = t3(1:length(timen));
        else
            w3 = interp1(time3, w3, timen, 'linear', 'extrap');
            t3 = interp1(time3, t3, timen, 'linear', 'extrap');
        end
    end

    if strcmp(grad_method, "center") % 10m
        % Weighted Finite Difference Method
        weight1 = (height(3)-height(2)) / (height(3)-height(1));
        weight3 = (height(2)-height(1)) / (height(3)-height(1));
        transport_temp = weight1 .* (w2.*w2.*t2 - w1.*w1.*t1)./(height(2)-height(1)) + weight3 .* (w3.*w3.*t3 - w2.*w2.*t2)./(height(3)-height(2));
    elseif strcmp(grad_method, "forward")
        transport_temp = (w2.*w2.*t2 - w1.*w1.*t1)./(height(2)-height(1));
    elseif strcmp(grad_method, "backward")
        transport_temp = (w3.*w3.*t3 - w2.*w2.*t2)./(height(3)-height(2));
    end
    transport = calc_1min_ave(transport_temp, samplingfreq);
end


function buoy = calc_1min_buoy(Tmi, ti, samplingfreq, timei, timen)
    if length(timei)~=length(timen)
        if length(timei)>length(timen)
            ti = ti(1:length(timen));
        else
            ti = interp1(timei, ti, timen, 'linear', 'extrap');
        end
    end
    
    g = 9.81;  %m/s^2 gravity
    buoy_temp = g ./ (Tmi + 273.15) .* ti .* ti;

    buoy = calc_1min_ave(buoy_temp, samplingfreq);
end

%% plot wt budget
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
            mt_value = 0.1;
        end
        subaxis(6, 2, pi, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'ML', 0.06, 'MB', 0.08, 'MR', 0.01, 'MT', mt_value);
        plot([buoy{1, ti, 3-hi+1}; buoy{2,ti, 3-hi+1}; buoy{3,ti, 3-hi+1}], 'r', 'LineWidth',2)
        hold on
        plot([shear{1,ti, 3-hi+1}; shear{2,ti, 3-hi+1}; shear{3,ti, 3-hi+1}], 'b', 'LineWidth',2)
        hold on
        plot([transport{1,ti, 3-hi+1}; transport{2,ti, 3-hi+1}; transport{3,ti, 3-hi+1}], 'c', 'LineWidth',2)
        hold on
        ylim([-2 10])
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
            legend('buoyancy', 'stratification', 'transport', 'box', 'off', 'Orientation','horizontal', "FontSize",20)
        end
        if pi==5
            ylabel("{\boldmath$\mathrm{\overline{w'T'}\hspace{1.7mm}Budget\hspace{1.7mm}Components\hspace{1.7mm}[K\hspace{1.7mm}m\hspace{1.7mm}s^{-2}]}$}", "FontSize",20, 'Interpreter','latex')
        end
    end
end