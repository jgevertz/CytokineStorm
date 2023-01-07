%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% Implementation of Cytokine Storm Model Local Sensitivity Analysis     %
% Written by Irina Kareva and Jana Gevertz. Last Update: 1/6/2023       %
% - Local sensivity analysis: vary all parameters up/down by vary*100%, %
%   where vary is a set of values defined by user on line -             %
% - Default initial conditions and parameter values are set from lines  %
%   23 to 50. Any can be changed by user.                               %
% - Dose/schedule are set at lines 52-56. Any can be changed by user,   %
%   though data in paper corresonds to:                                 %
%           p.doseD1 = 10; p.dosenumD1 = 3; p.intvlD1 = 24;             %
% - Outputs the relative impact of changing each parameter up and down  %
%   by vary*100% in a plot. Impact is studied on both tumor size at     %
%   t = 400, and max cytokine value over a 400-time-unit window         %
% - Then a final plot is made that shows the ranking (from 1, most      %
%   sensitive to 13, least sensitive) of each parameter, for each value %
%   of vary.                                                            %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc; 

%% Initial conditions
p.X0=0; 
p.Y0=0; 
p.Cin0=0; 
p.T0=0.1; 
in_cond=[p.X0 p.Y0 p.Cin0 p.T0]; 

%% Fixed params
% Immune cells: x equation 
p.y1 = 1; % y_L
p.y2 = 6; % y_U
p.b1 = 0.01; %b_yx
p.k1 = 2; % k_x, growth stimulated by killing tumor
p.kim = 0.693/7; % immune decay
p.k01 = 0.1; %drug inflow

% Cytokines: y equation 
p.Yin=1e-4; 
p.b2 = 0.8; %b_xy
p.c2 = 1;
p.k2 = .42;

% Tumor: T equation 
p.r = 0.06;
p.K = 1; 
p.c1 = 1.0;
p.k3 = 0.05; % k_T
p.allee = 0.05; % m

%% Drug dosing/schedule settings
p.doseD1 = 20; % Dose of drug
p.dosenumD1 = 1; % Number of doses to give
p.intvlD1 = 24; % Spacing between doses (assume time in hours) 
tStart1=200; % When to start treatment

%% Setting up ode solver
options=odeset('RelTol',1.0e-5);
tmax=400; 
last_dose = tStart1+p.dosenumD1*p.intvlD1;
if(last_dose>tmax)
    fprintf('Code will not work since trying to give drug until time %d but tmax = %d\n',...
        last_dose,tmax);
    stop
end

%% Baseline case
[time, Y]=cytokinestorm_v2(p,in_cond,tmax,tStart1);
T_baseline = Y(end,4); %tumor cells
Ymax_baseline = max(Y(:,2));
fprintf('Baseline scenario has T(400) = %f and max(Y) = %f\n',T_baseline,...
    Ymax_baseline); 

%% Local sensitivity analysis about baseline parameter values
vary = [0.01 0.5]; 
numParamsSens = 13; 
perc_change_T = zeros(length(vary),2*numParamsSens);
perc_change_Ymax = zeros(length(vary),2*numParamsSens);
Tidx = zeros(length(vary),2*numParamsSens);
Tsrt = zeros(length(vary),2*numParamsSens);
Ymaxidx = zeros(length(vary),2*numParamsSens);
Ymaxsrt = zeros(length(vary),2*numParamsSens);
ParamSens_T = zeros(length(vary),numParamsSens);
ParamSens_Ymax = zeros(length(vary),numParamsSens);
margins = [0.1,0.05]; % for subplot
for i = 1:length(vary)
    [perc_change_T(i,:), perc_change_Ymax(i,:)] = sensitivity(p,in_cond,tmax,...
        tStart1,vary(i),T_baseline,Ymax_baseline);
    %% Plot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.02, 0.7, 0.97]);
    subplot_tight(2,1,1,margins);
    bar(perc_change_T(i,:))
    xticks(1:1:length(perc_change_T(i,:)))
    xticklabels({'b_x_y^-','b_x_y^+','r^-','r^+','k_y^-','k_y^+','k_T^-',...
        'k_T^+','k_i_m^-','k_i_m^+','c_1^-','c_1^+','k_0_1^-','k_0_1^+',...
        'c_2^-','c_2^+','y_U^-','y_U^+','y_L^-','y_L^+','b_y_x^-','b_y_x^+',...
        'k_x^-','k_x^+','Y_i_n^-','Y_i_n^+'});
    ylabel('Relative change in T(400)')
    ylim([-inf,8])
    set(gca,'FontSize',16)

    subplot_tight(2,1,2,margins);
    bar(perc_change_Ymax(i,:))
    xticks(1:1:length(perc_change_Ymax(i,:)))
    xticklabels({'b_x_y^-','b_x_y^+','r^-','r^+','k_y^-','k_y^+','k_T^-',...
        'k_T^+','k_i_m^-','k_i_m^+','c_1^-','c_1^+','k_0_1^-','k_0_1^+',...
        'c_2^-','c_2^+','y_U^-','y_U^+','y_L^-','y_L^+','b_y_x^-','b_y_x^+',...
        'k_x^-','k_x^+','Y_i_n^-','Y_i_n^+'});
    ylabel('Relative change in max(Y)')
    set(gca,'FontSize',16)
    sgtitle(['Vary parameter by \pm' num2str(vary(i))],'FontSize',18,'FontWeight','bold');
    
    %% Rank parameters by most important
    [Tsrt(i,:), Tidx(i,:)] = sort(abs(perc_change_T(i,:)),'descend');
    [Ymaxsrt(i,:), Ymaxidx(i,:)] = sort(abs(perc_change_Ymax(i,:)),'descend');
    
    %% Rank by considering first appearance of parameter on ordered list
    ParamSens_T(i,:) = rank_parameters(Tidx(i,:));
    ParamSens_Ymax(i,:) = rank_parameters(Ymaxidx(i,:));
end

%% Now give a each parameter a label from 1-13
ParamSens_T = ceil(ParamSens_T/2);
ParamSens_Ymax = ceil(ParamSens_Ymax/2);
ParamRanking_T = zeros(size(ParamSens_T));
for i = 1:size(ParamSens_T,2)
    for j = 1:size(ParamSens_T,1)
        [val,index] = find(ParamSens_T(j,:)==i,1);
        ParamRanking_T(j,i) = index;
    end
end

ParamRanking_Ymax = zeros(size(ParamSens_Ymax));
for i = 1:size(ParamSens_Ymax,2)
    for j = 1:size(ParamSens_Ymax,1)
        [val,index] = find(ParamSens_Ymax(j,:)==i,1); 
        ParamRanking_Ymax(j,i) = index;
    end
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.75, 0.6]);
subplot(1,2,1)
b = bar(ParamRanking_T');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','Color','b')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','Color','r')
xticks(1:1:size(ParamRanking_T,2))
xticklabels({'b_x_y','r','k_y','k_T','k_i_m','c_1','k_0_1',...
        'c_2','y_U','y_L','b_y_x','k_x','Y_i_n'});
xlabel('Parameter','FontSize',16)
ylabel('Rank: Tumor Sensitivity','FontSize',16)
a = get(gca,'YTickLabel'); % or use XTickLabel, does same thing 
set(gca,'YTickLabel',a,'FontSize',14)
legend('1% variation','50% variation','Location','NorthWest','FontSize',16);

subplot(1,2,2)
b = bar(ParamRanking_Ymax');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','Color','b')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','Color','r')
xticks(1:1:size(ParamRanking_Ymax,2))
xticklabels({'b_x_y','r','k_y','k_T','k_i_m','c_1','k_0_1',...
        'c_2','y_U','y_L','b_y_x','k_x','Y_i_n'});
xlabel('Parameter','FontSize',16)
ylabel('Rank: Maximum Cytokine Sensitivity','FontSize',16)
a = get(gca,'YTickLabel'); % or use XTickLabel, does same thing 
set(gca,'YTickLabel',a,'FontSize',14)
sgtitle('Parameters Ranked from Most (Rank 1) to Least (Rank 13) Sensitive',...
    'FontSize',16,'FontWeight','bold')
fname_fig = 'parameter_sensitivity';
saveas(gcf,[fname_fig,'.fig']);
saveas(gcf,[fname_fig,'.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ParamSens = rank_parameters(Idx)
    ParamSens = []; remove = []; 
    count = 1;
    for j = 1:length(Idx)
        curr = Idx(j); 
        A = find(remove == curr, 1);
        if(isempty(A)==1) % not in remove list
            ParamSens(count) = curr;
            if mod(curr,2) == 1 % odd
                remove(count) = curr+1;
                count = count+1;
            else % even
                remove(count) = curr-1; 
                count = count+1;
            end
        end
    end
end

function [perc_change_T,perc_change_Ymax] = sensitivity(p,in_cond,tmax,...
    tStart1,vary,T_baseline,Ymax_baseline)
    Tfinal = [];
    Ymax_final = [];
    count = 1;

    %% Parameters not studied: m, K
    b2_baseline = p.b2;  %b_ic = b_xy
    p.b2 = b2_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With b2 = %f, Tfinal = %f, Ymax_final = %f\n',p.b2,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.b2 = b2_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With b2 = %f, Tfinal = %f, Ymax_final = %f\n',p.b2,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.b2 = 0.8; % reset

    r_baseline = p.r;
    p.r = r_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With r = %f, Tfinal = %f, Ymax_final = %f\n',p.r,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.r = r_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With r = %f, Tfinal = %f, Ymax_final = %f\n',p.r,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.r = 0.06; % reset

    k2_baseline = p.k2; %k_y
    p.k2 = k2_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k2 = %f, Tfinal = %f, Ymax_final = %f\n',p.k2,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.k2 = k2_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k2 = %f, Tfinal = %f, Ymax_final = %f\n',p.k2,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.k2 = .42; % reset

    k3_baseline = p.k3; %k_T
    p.k3 = k3_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k3 = %f, Tfinal = %f, Ymax_final = %f\n',p.k3,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.k3 = k3_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k3 = %f, Tfinal = %f, Ymax_final = %f\n',p.k3,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.k3 = 0.05; % reset

    kim_baseline = p.kim; 
    p.kim = kim_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With kim = %f, Tfinal = %f, Ymax_final = %f\n',p.kim,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.kim = kim_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With kim = %f, Tfinal = %f, Ymax_final = %f\n',p.kim,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.kim = 0.693/7;% reset

    c1_baseline = p.c1; 
    p.c1 = c1_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With c1 = %f, Tfinal = %f, Ymax_final = %f\n',p.c1,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.c1 = c1_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With c1 = %f, Tfinal = %f, Ymax_final = %f\n',p.c1,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.c1 = 1.0; % reset

    k01_baseline = p.k01; 
    p.k01 = k01_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k01 = %f, Tfinal = %f, Ymax_final = %f\n',p.k01,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.k01 = k01_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k01 = %f, Tfinal = %f, Ymax_final = %f\n',p.k01,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.k01 = 0.1; % reset

    c2_baseline = p.c2; 
    p.c2 = c2_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With c2 = %f, Tfinal = %f, Ymax_final = %f\n',p.c2,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.c2 = c2_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With c2 = %f, Tfinal = %f, Ymax_final = %f\n',p.c2,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.c2 = 1; % reset

    y2_baseline = p.y2; %y_U
    p.y2 = y2_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With y2 = %f, Tfinal = %f, Ymax_final = %f\n',p.y2,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.y2 = y2_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With y2 = %f, Tfinal = %f, Ymax_final = %f\n',p.y2,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.y2 = 6; % reset

    y1_baseline = p.y1; % yL
    p.y1 = y1_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With y1 = %f, Tfinal = %f, Ymax_final = %f\n',p.y1,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.y1 = y1_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With y1 = %f, Tfinal = %f, Ymax_final = %f\n',p.y1,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.y1 = 1;  % reset

    b1_baseline = p.b1;   %b_ci = b_yx
    p.b1 = b1_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With b1 = %f, Tfinal = %f, Ymax_final = %f\n',p.b1,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.b1 = b1_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With b1 = %f, Tfinal = %f, Ymax_final = %f\n',p.b1,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.b1 = 0.01; % reset

    k1_baseline = p.k1; %k_x
    p.k1 = k1_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k1 = %f, Tfinal = %f, Ymax_final = %f\n',p.k1,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.k1 = k1_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With k1 = %f, Tfinal = %f, Ymax_final = %f\n',p.k1,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.k1 = 2; % reset

    Yin_baseline = p.Yin; 
    p.Yin = Yin_baseline*(1-vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With Yin = %f, Tfinal = %f, Ymax_final = %f\n',p.Yin,Tfinal(count),...
        Ymax_final(count));
    count = count+1; 
    p.Yin = Yin_baseline*(1+vary); 
    [time,Y] = cytokinestorm_v2(p,in_cond,tmax,tStart1);
    Tfinal(count) = Y(end,4); %tumor cells
    Ymax_final(count) = max(Y(:,2));
    fprintf('With Yin = %f, Tfinal = %f, Ymax_final = %f\n',p.Yin,Tfinal(count),...
        Ymax_final(count));
    count = count+1;
    p.Yin=1e-4; % reset
    
    perc_change_T = (Tfinal-T_baseline)/T_baseline;
    perc_change_Ymax = (Ymax_final-Ymax_baseline)/Ymax_baseline;
end 

function [Time, X] = cytokinestorm_v2(p,y0,tspan,tStart1)
    %% Model
    function dydt = ode(t,y0)
        x   = y0(1); 
        y   = y0(2); 
        Cin = y0(3);
        T   = y0(4); 
      
        % immune cells
        dx = p.k01*Cin...
            +p.k1*p.k3*x*T/(p.c1+x)...
            +p.b1*(x/(p.c2+x))*(y-p.y1)*(p.y2-y)... %or x^2(y-1)(6-y)
            -p.kim*x;
      
        dy   = p.Yin-p.k2*y+p.b2*y*(x/(p.c2+x)); % cytokines
        dCin = -p.k01*Cin; % drug
        dT = p.r*T*(T-p.allee)*(1-T/p.K)- p.k3*x*T/(p.c1+x); % tumor
        
        dydt = [dx;dy;dCin;dT];
        
    end
    %% Dosing
    ivdoseD1 = p.doseD1; % Drug 1
    % drug 2 - not used here
    tStart2 = 0.001;
    p.doseD2 = 0;% mpk
    p.dosenumD2 = 1;        % dosenumD2 here, was 5
    p.intvlD2 = 168;%hours   % intvlD2 here, was 24
    ivdoseD2 = p.doseD2;

    %% Auxiliary "keystone" variable from the solution
    xx0  = y0(1); 
    yy0  = y0(2);
    Cin0 = y0(3);
    T0 = y0(4);
    for i = 1:p.dosenumD1
        tDose1(i) = p.intvlD1*(i-1) + tStart1;
        cDose1(i) = ivdoseD1;
    end
        
    % If you change p.dosenumD2 to something other than
    % p.dosenumD1, then the lengths of tDose2 and cDose2
    % are different from tDose1 and cDose1
    for i = 1:p.dosenumD2
        tDose2(i) = p.intvlD2*(i-1) + tStart2;
        cDose2(i) = ivdoseD2;
    end

    % Your time interval is the UNION of both tDose1 and tDose2
    % Put the 2 sets together, in order, no duplicates
    tUnion = union( tDose1, tDose2, 'sorted');
    tUnion = [0 tUnion];
    tLength = length(tUnion);
    ti = [tUnion max(tspan)];
    
    % Run the for the union of time; step through each interval
    tPoints = []; cPoints = [];
    for i = 1:tLength 
        % Given a TIME, does that TIME exist in tDose1 array?
        [lia1,locb1] = ismember( tUnion(i), tDose1 );
        if lia1 == 1
            Cin0 = Cin0 + cDose1(locb1);
        end
        
   
        % Build the tSpan and run the Diff Eq
        tSpan = [ti(i) ti(i+1)];
        [tOut,cOut] = ode23s(@ode,tSpan,[xx0 yy0 Cin0 T0]);
        
        % Concatenate the results of the ode45 work
        tPoints = [tPoints; tOut];
        cPoints = [cPoints; cOut];
        
        % Set new init conditions for ode45 and repeat the loop
        xx0 = cOut(end,1); yy0 = cOut(end,2); ...
        Cin0 = cOut(end,3); T0 = cOut(end,4); 
    end

    Time = tPoints;
    X = cPoints;
end