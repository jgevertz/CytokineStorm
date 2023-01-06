%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% Implementation of Cytokine Storm Model and Analysis                   %
% Written by Irina Kareva and Jana Gevertz. Last Update: 1/6/2023       %
% - Default initial conditions and parameter values are set from lines  %
%   33 to 60. Any can be changed by user.                               %
% - Dose/schedule are set at lines 62-66. Any can be changed by user.   %
% - Enter 0 at prompt to simply solve the model at the entered          %
%   parameter values, ICs, and dosing schedule. Plots of the dynamics   %
%   of each variable will be created.                                   %
% - Enter 1 at prompt to perform a  sweep over rate of tumor killing    %
%   (k_T) and rate of cytokine stimulation (b_xy) parameter space. For  %
%   each parameterization tested:                                       %
%   1) Assess if tumor was effectively treated. We say treatment is     %
%      effective if tumor size is <=0.1*carrying capacity = 0.1.        %
%   2) Assess if a cytokine storm occurred using one of two             %
%   definitions. Definition D1: a storm occurs if the absolute value of %
%   cytokines levels crosses the threshold value of y_U. Definition D2: %
%   a storm occurs if the maximal rate of change of cytokine level      %
%   crosses the threshold value of storm_cutoff_rate = 1.4.             %
% - If 1 is entered at prompt, the code will also analyze a single      %
%   paremeterization that falls in the "effective treatment with        %
%   cytokine storm" regime and ask: can we find a different             %
%   administration schedule that preserves efficacy while not causing a %
%   storm?                                                              %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc; tic;
prompt = "Do you want to solve for just one set of parameters (enter 0), or do a parameter sweep (enter 1)? ";
decision = input(prompt);

%% Initial conditions
p.X0=0; 
p.Y0=0; 
p.Cin0=0; 
p.T0=0.1; 
in_cond=[p.X0 p.Y0 p.Cin0 p.T0]; %initial conditions 

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
p.doseD1 = 10; % Dose of drug
p.dosenumD1 = 1; % Number of doses to give
p.intvlD1 = 24; % Spacing between doses (assume time in hours) 
tStart1=200; % When to start treatment

%% Setting up ode solver
options=odeset('RelTol',1.0e-5);
tmax=3600; 
last_dose = tStart1+p.dosenumD1*p.intvlD1;
if(last_dose>tmax)
    fprintf('Code will not work since trying to give drug until time %d but tmax = %d\n',...
        last_dose,tmax);
    stop
end

if decision == 0 % Just plots solution curves
    fprintf('Running %d doseages of dose = %f, with spacing between doses of %d\n',...
       p.dosenumD1, p.doseD1, p.intvlD1); 
    [time, Y]=cytokinestorm_v2(p,in_cond,tmax,tStart1);
    immune   = Y(:,1);  
    cytokine = Y(:,2);  
    Cin      = Y(:,3);
    T        = Y(:,4); %tumor cells
    
    %% Plot solution
    figure; 
    subplot(3,2,1) % immune 
    plot(time,immune,'linewidth',2); 
    xlabel('time')
    ylabel('immune')
    title('immune')

    subplot(3,2,2) % cytokine
    yline(p.y1,'r:','linewidth',2); hold on;
    yline(p.y2,'k:','linewidth',2); 
    plot(time,cytokine,'linewidth',2); hold off;
    xlabel('time')
    ylabel('cytokine')
    title('cytokines')
    legend('y1','y2','cytokine')

    subplot(3,2,3) % tumor
    plot(time,T,'linewidth',2); 
    xlabel('time')
    ylabel('tumor')

    subplot(3,2,4) % phase parameter plot
    yline(p.y1,'r:','linewidth',2); hold on;
    yline(p.y2,'k:','linewidth',2);
    plot(immune,cytokine,'linewidth',2); hold off;
    xlabel('immune')
    ylabel('cytokine')
    title ('phase parameter portrait')
    legend('y1','y2','ppp')

    subplot(3,2,5) % drug
    plot(time,Cin,'linewidth',2); 
    xlabel('time')
    ylabel('C_i_n')

elseif decision == 1
    %% Sets parameter range to search, and discretization
    num_pts = 20; % can sample more or less points in range by changing this
	effective_cutoff = 0.1*p.K; % cutoff set to 1/10 of carrying capacity
	storm_cutoff = max(p.y1,p.y2); % above larger threshold is "storm"
    storm_cutoff_rate = 1.4;  % determined looking at min value of rate on cytokine heatmap
                              % for the "surface" that formed in space      
    
    % Rate of tumor killing, k3 = k_T, varies from 0.004 to 0.164                         
    k3_min = 0.004; k3_max = 0.164; 
    k3_step = (k3_max-k3_min)/num_pts;
    k3 = linspace(k3_min,k3_max,num_pts+1);   
    
    % Rate of cytokine stimulation, b2 = b_xy, varies from 0.4 to 2
    b2_min = 0.4; b2_max = 2; 
	b2_step = (b2_max-b2_min)/num_pts;
	b2 = linspace(b2_min,b2_max,num_pts+1); 
        
    Tfinal2 = zeros(num_pts,num_pts);
    max_cytokine2 = zeros(num_pts,num_pts); % For storm D1
    param_regime2 = zeros(num_pts,num_pts); % For storm D1
    
    max_rate_cyotkine2 = zeros(num_pts,num_pts); % For storm D3
    param_regime_rate2 = zeros(num_pts,num_pts); % For storm D3
    
    %% Loop through all parameter combinations of b2 and k3
    for i=1:length(b2)
        for j = 1:length(k3)
            p.b2 = b2(i);
            p.k3 = k3(j);
           
            [time, Y]=cytokinestorm_v2(p,in_cond,tmax,tStart1);
            immune   = Y(:,1);  
            cytokine = Y(:,2);  
            Cin      = Y(:,3);
            T        = Y(:,4); %tumor cells
            Tfinal2(i,j) = T(end); % final tumor size
            
            %% Definition D1 of Cytokine Storm: Max absolute value of cytokines
            max_cytokine2(i,j) = max(cytokine); 
            fprintf('Cytokine stimulation rate of = %f, rate of tumor killing = %f has:\n',...
                p.b2,p.k3);
            fprintf('\tTfinal = %f\n\tmax(cytokine) = %f\n',...
                Tfinal2(i,j),max_cytokine2(i,j)); 
            
            % Determine which regime we are in for D1
            if( (Tfinal2(i,j)<=effective_cutoff) && (max_cytokine2(i,j)<=storm_cutoff) )
                % P1) Effective treatment with no cytokine storm (best case)
                param_regime2(i,j) = 1; 
            elseif( (Tfinal2(i,j)<=effective_cutoff) && (max_cytokine2(i,j)>storm_cutoff) )
                % P2) Effective treatment but cytokine storm occurs 
                param_regime2(i,j) = 2;
            elseif( (Tfinal2(i,j)>effective_cutoff) && (max_cytokine2(i,j)<=storm_cutoff) )
                % P3) Ineffective treatment, no cytokine storm   
                param_regime2(i,j) = 3;
            else
                % P4) Ineffective treatment, cyokine storm occurs (worst case)
                param_regime2(i,j) = 4;
            end
            fprintf('\tD1) Parameter Regime: %d\n',param_regime2(i,j));
            
            % Definition D3 of Cytokine Storm: Max rate of change of cytokines
            cytokine_rates = diff(cytokine)./diff(time);
            max_rate_cytokine2(i,j) = max(cytokine_rates); 
            fprintf('\tTfinal2 = %f\n\tmax(cytokine_rate) = %f\n',...
                Tfinal2(i,j),max_rate_cytokine2(i,j)); 
            
            % Determine which regime we are in for D3                                % space
            if( (Tfinal2(i,j)<=effective_cutoff) && (max_rate_cytokine2(i,j)<=storm_cutoff_rate) )
                % P1) Effective treatment with no cytokine storm (best case)
                param_regime_rate2(i,j) = 1; 
            elseif( (Tfinal2(i,j)<=effective_cutoff) && (max_rate_cytokine2(i,j)>storm_cutoff_rate) )
                % P2) Effective treatment but cytokine storm occurs 
                param_regime_rate2(i,j) = 2;
            elseif( (Tfinal2(i,j)>effective_cutoff) && (max_rate_cytokine2(i,j)<=storm_cutoff_rate) )
                % P3) Ineffective treatment, no cytokine storm   
                param_regime_rate2(i,j) = 3;
            else
                % P4) Ineffective treatment, cyokine storm occurs (worst case)
                param_regime_rate2(i,j) = 4;
            end
            fprintf('\tD3) Parameter Regime: %d\n',param_regime_rate2(i,j));
                 
        end
    end
    
    %% Make heatmap of final tumor size 
    figure;
    h1 = heatmap(b2,k3,Tfinal2');
    xlabel('Rate of cytokine stimulation'); % b2
    ylabel('Rate of tumor killing');        % k3
    grid off;
    h1.Title= 'Final tumor size';
    
    %% Make heatmap of max cytokine value
    figure;
    h1 = heatmap(b2,k3,max_cytokine2');
    xlabel('Rate of cytokine stimulation'); % b2
    ylabel('Rate of tumor killing');        % k3
    grid off;
    h1.Title= 'Maximum cytokine value';
    
    %% Make heatmap of parameter regime according to max cytokine value
    map = [1 1 1 % white
    1 1 0  % yellow
    1 0 0 % red
    0 0 0]; % black    
    figure;
    imagesc(b2,k3,param_regime2')
    colormap(map)
    cbh = colorbar(); 
    set(cbh, 'YTick', [1, 2, 3, 4], ...
    'YTickLabel', {'P1', 'P2', 'P3', 'P4'})
    xlabel('Rate of cytokine stimulation'); % b2
    ylabel('Rate of tumor killing');        % k3
    title('D1 Parameter Regime');

    %% Make heatmap of max cytokine rate
    figure;
    h1 = heatmap(b2,k3,max_rate_cytokine2');
    xlabel('Rate of cytokine stimulation'); % b2 = b_xy
    ylabel('Rate of tumor killing');        % k3 = k_T
    grid off;
    h1.Title= 'Maximum rate of cytokine change';

    %% Make heatmap of parameter regime according to max cytokine rate
    figure;
    imagesc(b2,k3,param_regime_rate2')
    colormap(map)
    cbh = colorbar(); 
    set(cbh, 'YTick', [1, 2, 3, 4], ...
    'YTickLabel', {'P1', 'P2', 'P3', 'P4'})
    xlabel('Rate of cytokine stimulation'); % b2
    ylabel('Rate of tumor killing');        % k3
    title('D3 Parameter Regime');    
    
    %% Now pick a point in P2 (effective treatment, storm): done manually
    p.b2 = 2; 
    p.k3 = 0.15;
    
    %% Preserve total adminstered dose at 10, and vary number and spacing of doses
    total_drug = p.doseD1;
    num_pts2 = 7; 
    doseNum_min = 2;
    doseNum_step = doseNum_min;
    doseNum_max = num_pts2*doseNum_step+doseNum_min;
    doseNum = linspace(doseNum_min,doseNum_max,num_pts2+1); 
    
    interval_min = 24;
    interval_step =interval_min;
    interval_max = num_pts2*interval_step+interval_min;
    interval = linspace(interval_min,interval_max,num_pts2+1); 
    Tfinal_control2 = zeros(num_pts2,num_pts2);
    max_cytokine_control2 = zeros(num_pts2,num_pts2); % For storm D1
    param_regime_control2 = zeros(num_pts2,num_pts2); % For storm D1
    
    %% Loop through set of dose numbers and spacings
    for i=1:length(doseNum)
        for j = 1:length(interval)
            p.intvlD1 = interval(j);%hours   % intvlD1 here
            p.dosenumD1 = doseNum(i);        % dosenumD1 here
            p.doseD1 = total_drug/p.dosenumD1;
            fprintf('i = %d, j = %d: doseNum = %d, interval = %d, dose = %f has:\n',...
                i,j,doseNum(i),interval(j),p.doseD1); 
            [time, Y]=cytokinestorm_v2(p,in_cond,tmax,tStart1);
            immune   = Y(:,1);  
            cytokine = Y(:,2);  
            Cin      = Y(:,3);
            T        = Y(:,4); %tumor cells
            Tfinal_control2(i,j) = T(end); 
            
            %% Definition D1 of Cytokine Storm: Max absolute value of cytokines
            max_cytokine_control2(i,j) = max(cytokine); 
            fprintf('\tTfinal = %f\n\tmax(cytokine) = %f\n',...
                Tfinal_control2(i,j),max_cytokine_control2(i,j)); 
            
            if( (Tfinal_control2(i,j)<=effective_cutoff) && (max_cytokine_control2(i,j)<=storm_cutoff) )
                % P1) Effective treatment with no cytokine storm (best case)
                param_regime_control2(i,j) = 1; 
            elseif( (Tfinal_control2(i,j)<=effective_cutoff) && (max_cytokine_control2(i,j)>storm_cutoff) )
                % P2) Effective treatment but cytokine storm occurs 
                param_regime_control2(i,j) = 2;
            elseif( (Tfinal_control2(i,j)>effective_cutoff) && (max_cytokine_control2(i,j)<=storm_cutoff) )
                % P3) Ineffective treatment, no cytokine storm   
                param_regime_control2(i,j) = 3;
            else
                % P4) Ineffective treatment, cyokine storm occurs (worst case)
                param_regime_control2(i,j) = 4;
            end
            fprintf('\tD1) Parameter Regime: %d\n',param_regime_control2(i,j));
        end
    end
    
    %% Make heatmap of final tumor size in dose number-spacing space 
    figure;
    h1 = heatmap(doseNum,interval,Tfinal_control2');
    xlabel('Number of Doses'); 
    ylabel('Spacing between doses (hours)');                      
    grid off;
    h1.Title= 'Final tumor size';
    
    %% Make heatmap of max cytokine value in dose number-spacing space 
    figure;
    h1 = heatmap(doseNum,interval,max_cytokine_control2');
    xlabel('Number of Doses'); 
    ylabel('Spacing between doses (hours)');          
    grid off;
    h1.Title= 'Maximum cytokine value';
      
    %% Make heatmap of parameter regime according to max cytokine value in dose number-spacing space 
    figure; 
    imagesc(doseNum,interval,param_regime_control2')
    colormap(map)
    cbh = colorbar(); 
    set(cbh, 'YTick', [1, 2, 3, 4], ...
    'YTickLabel', {'P1', 'P2', 'P3', 'P4'})    
    xlabel('Number of Doses'); 
    ylabel('Spacing between doses (hours)');
    title('At fixed cumulative dose, can we control storm but preserve efficacy?');
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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