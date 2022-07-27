%% mt_Econ501B_HW1.m
% is a MATLAB code for Question 4 of the Econ 501B's assignment 1.
% 
% ------------------------------------------------------------------------
% Copyright (c) 2021 Masaki Tanaka (Washington University in St Louis)
% All rights researved.
%
% This source code or any portion thereof must not be reproduced or used
% in any manner whatsoever.
%
% Create: Nov 234, 2021
% Revised: Nov 26, 2021
%  -- Debug
% Revised: Nov 27, 2021
%  -- Modified some typos
%


%% 0. Preamble ***********************************************************
close all;
clear;
clc;

global l_vec r w beta  a_vec V0_mat l_j Pi_j a_k 

%% 1. Set Parameters *****************************************************

% set model parameters
l_vec         = [0.1, 1.0, 1.25]';
r_original    = 0.0505;
w             = 0.0560;
beta          = 0.95;
a_lb_original = 0.10; 
               % max borrowing
               % caution: declare as a positive value. The min of a will be
               % set to - a_lb.
a_limit = 5.00; % computational upper limit of asset holding

Pi = [0.4, 0.5, 0.1;
      0.2, 0.7, 0.1;
      0.2, 0.1, 0.7]; % transision matrix

% simulation setting
r_vec    = [0.0305, 0.0705];
a_lb_vec = [0.2 , 0.0];

% set computation parameters
nGrid   = 500     ; % # of grids for a (current asset holdings)
tol     = 10^(-5) ; % tolerance level
maxiter = 500     ; % maximum number of iterations
options = optimset('Display','off'); % Options for fminbnd

% set graphics parameters
ax_min  = -0.3;
ax_max  =  0.5;

%% 2. Value-function iterations (VFIs): baseline case  *******************
r    = r_original;    % set r to the baseline value
a_lb = a_lb_original; % set a_lb to the baseline value

a_vec  = linspace(-a_lb, a_limit, nGrid); % Samples of a
V0_mat = zeros(length(l_vec), nGrid); % the initial values for VFIs

V1_mat    = zeros(length(l_vec), nGrid);
a_prm_mat = zeros(length(l_vec), nGrid);

dif = 1; % difference between the updated V and the old V.
i   = 1; % iteration counter
while (dif > tol) & (i < maxiter)
    fprintf('\n')
    fprintf('%s:: Iteration %d \n', mfilename, i)

    for j = 1 : length(l_vec) % for-loop wrt current l
        l_j = l_vec(j);
        Pi_j = Pi(j ,:);

        for k = 1 : length(a_vec) % for-loop wrt current a
            % fminbnd is a MATLAB built-in optimization routine
            % which mimimizes a single-value function on a fixed interval.
            % So, we cannot do optimization on all the grid points at
            % the same time, and have to use for-loop. 
            a_k = a_vec(k);
            a_ub = (1 + r) * a_k + w * l_j;
            if a_ub > a_limit
                a_ub = a_limit;
            end
            [a_prm_jk, minusV_jk] = fminbnd(@ValFunc_Aiyagari, -a_lb, ...
                                    a_ub, options);
            a_prm_mat(j, k)       = a_prm_jk;
            V1_mat(j, k)          = - minusV_jk;
        end
    end
    dif   = max(max(abs(V1_mat - V0_mat))); % Difference betweeen the updated 
                                            % V and the old V.
                                            % Differences are calculated on
                                            % each (a, l). And we take the max. 
    V0_mat = V1_mat; % Overwrite the old V with the updated V
    i      = i + 1;  % Increase the iteration counter.
end

if i < maxiter
    fprintf('\n')
    disp([mfilename, ':: The value function iterations stopped because ',...
        'the value function converged.'])
    fprintf('\n')
else
    fprintf('\n')
    disp([mfilename, ':: The value function iterations stopped because ',...
        'the number of trials reached the max of iterations.'])
    fprintf('\n')
end

V0_mat(a_prm_mat > 4.99)    = NaN      ; % boundary solution -> NaN
a_prm_mat(a_prm_mat > 4.99) = NaN      ; % boundary solution -> NaN

V_baseline     = V0_mat   ; % store the result
a_prm_baseline = a_prm_mat; % store the result
a_baseline     = a_vec    ; % store the result

%% 3. Value-function iterations (VFIs): with different values of r ********
a_lb        = a_lb_original; % set a_lb to the baseline value
V_r_chg     = zeros(length(l_vec), nGrid, length(r_vec));
a_prm_r_chg = zeros(length(l_vec), nGrid, length(r_vec));
a_r_chg     = zeros(nGrid, length(r_vec));

for m = 1 : length(r_vec)
    r = r_vec(m); % set r to an alternative value

    a_vec     = linspace(-a_lb, a_limit, nGrid); % Samples of a
    V0_mat    = zeros(length(l_vec), nGrid); % the initial values for VFIs

    V1_mat    = zeros(length(l_vec), nGrid);
    a_prm_mat = zeros(length(l_vec), nGrid);
    dif       = 1; % difference between the updated V and the old V.
    i         = 1; % iteration counter

    while (dif > tol) & (i < maxiter)
        fprintf('\n')
        fprintf('%s:: Iteration %d \n', mfilename, i)

        for j = 1 : length(l_vec) % for-loop wrt current l
            l_j = l_vec(j);
            Pi_j = Pi(j ,:);

            for k = 1 : length(a_vec) % for-loop wrt current a
                a_k = a_vec(k);

                a_ub = (1 + r) * a_k + w * l_j;
                if a_ub > a_limit
                    a_ub = a_limit;
                end

                if a_ub > -a_lb
                    [a_prm_jk, minusV_jk] = fminbnd(@ValFunc_Aiyagari,...
                                        -a_lb, a_ub, options);
                    a_prm_mat(j, k) = a_prm_jk;
                    V1_mat(j, k)    = - minusV_jk;
                else
                    a_prm_mat(j, k) = -Inf;
                    V1_mat(j, k)    = -Inf;
                end
             end
        end

        dif   = max(max(abs(V1_mat - V0_mat)));
        V0_mat = V1_mat; % Overwrite the old V with the updated V
        i      = i + 1;  % Increase the iteration counter.
    end

    if i < maxiter
        fprintf('\n')
        disp([mfilename, ':: The value function iterations stopped because ',...
            'the value function converged.'])
        fprintf('\n')
    else
        fprintf('\n')
        disp([mfilename, ':: The value function iterations stopped because ',...
            'the number of trials reached the max of iterations.'])
        fprintf('\n')
    end

    a_prm_mat(V0_mat == -Inf) = NaN      ; % penalty value -> NaN
    V0_mat(V0_mat == -Inf)    = NaN      ; % penalty value -> NaN
    V0_mat(a_prm_mat > 4.99)    = NaN      ; % boundary solution -> NaN
    a_prm_mat(a_prm_mat > 4.99) = NaN      ; % boundary solution -> NaN

    V_r_chg(:, :, m)     = V0_mat   ; % store the result
    a_prm_r_chg(:, :, m) = a_prm_mat; % store the result
    a_r_chg(:, m)        = a_vec    ; % store the result
end

%% 4. Value-function iterations (VFIs): with different values of a_lb *****
r              = r_original; % set r to the baseline value
V_a_lb_chg     = zeros(length(l_vec), nGrid, length(r_vec));
a_prm_a_lb_chg = zeros(length(l_vec), nGrid, length(r_vec));
a_a_lb_chg     = zeros(nGrid, length(r_vec));

for m = 1 : length(a_lb_vec)
    a_lb = a_lb_vec(m); % set a_lb to an alternative value

    a_vec  = linspace(-a_lb, a_limit, nGrid); % Samples of a
    V0_mat = zeros(length(l_vec), nGrid); % the initial values for VFIs

    V1_mat    = zeros(length(l_vec), nGrid);
    a_prm_mat = zeros(length(l_vec), nGrid);
    dif = 1; % difference between the updated V and the old V.
    i   = 1; % iteration counter

    while (dif > tol) & (i < maxiter)
        fprintf('\n')
        fprintf('%s:: Iteration %d \n', mfilename, i)

        for j = 1 : length(l_vec) % for-loop wrt current l
            l_j = l_vec(j);
            Pi_j = Pi(j ,:);

            for k = 1 : length(a_vec) % for-loop wrt current a
                a_k = a_vec(k);
                a_ub = (1 + r) * a_k + w * l_j;
                if a_ub > a_limit
                    a_ub = a_limit;
                end

                if a_ub > -a_lb
                    [a_prm_jk, minusV_jk] = fminbnd(@ValFunc_Aiyagari,...
                                        -a_lb, a_ub, options);
                    a_prm_mat(j, k) = a_prm_jk;
                    V1_mat(j, k)    = - minusV_jk;
                else
                    a_prm_mat(j, k) = -Inf;
                    V1_mat(j, k)    = -Inf;
                end
            end
        end

        dif   = max(max(abs(V1_mat - V0_mat)));
        V0_mat = V1_mat; % Overwrite the old V with the updated V
        i      = i + 1;  % Increase the iteration counter.
    end

    if i < maxiter
        fprintf('\n')
        disp([mfilename, ':: The value function iterations stopped because ',...
            'the value function converged.'])
        fprintf('\n')
    else
        fprintf('\n')
        disp([mfilename, ':: The value function iterations stopped because ',...
            'the number of trials reached the max of iterations.'])
        fprintf('\n')
    end

    a_prm_mat(V0_mat == -Inf) = NaN      ; % penalty value -> NaN
    V0_mat(V0_mat == -Inf)    = NaN      ; % penalty value -> NaN
    V0_mat(a_prm_mat > 4.99)    = NaN      ; % boundary solution -> NaN
    a_prm_mat(a_prm_mat > 4.99) = NaN      ; % boundary solution -> NaN


    V_a_lb_chg(:, :, m)      = V0_mat   ; % store the result
    a_prm_a_lb_chg(:, :, m)  = a_prm_mat; % store the result
    a_a_lb_chg(:, m)         = a_vec    ; % store the result
end


%% 5. Graphics ***********************************************************

% 5-1) baseline case:
gcf = figure;

subplot(2, 1, 1)
hold on
plot(a_baseline, V_baseline(1, :), 'red')
plot(a_baseline, V_baseline(2, :), 'blue')
plot(a_baseline, V_baseline(3, :), 'green')
xlabel('a')
ylabel('V')
xlim([ax_min, ax_max]);
legend({'l = 0.1','l = 1.0', 'l = 1.25'}, 'Location','northwest')
title('Value function')
hold off

subplot(2, 1, 2)
hold on
plot(a_baseline, a_prm_baseline(1, :), 'red')
plot(a_baseline, a_prm_baseline(2, :), 'blue')
plot(a_baseline, a_prm_baseline(3, :), 'green')
xlabel('a')
ylabel('a prime')
xlim([ax_min, ax_max]);
ylim([ax_min, ax_max]);
legend({'l = 0.1','l = 1.0', 'l = 1.25'},'Location','northwest')
title('Policy function')
hold off

print(gcf, 'Result_baseline','-dpng') %'-dpng':PNG, '-dpdf':PDF
close(gcf)    


% 5-2) interest rate:
gcf = figure;
subplot(3, 1, 1)
hold on
plot(a_baseline, a_prm_baseline(1, :), 'red')
plot(a_r_chg(:,1), a_prm_r_chg(1,:,1), 'blue')
plot(a_r_chg(:,2), a_prm_r_chg(1,:,2), 'green')
xlabel('a')
ylabel('a prime')
xlim([ax_min, ax_max]);
ylim([ax_min, ax_max]);
legend({sprintf('r = %.4f', r_original),...
        sprintf('r = %.4f', r_vec(1)),...
        sprintf('r = %.4f', r_vec(2))}, 'Location','northwest')
title(sprintf('Policy function: l = %.2f', l_vec(1)))
hold off

subplot(3, 1, 2)
hold on
plot(a_baseline, a_prm_baseline(2, :), 'red')
plot(a_r_chg(:,1), a_prm_r_chg(2,:,1), 'blue')
plot(a_r_chg(:,2), a_prm_r_chg(2,:,2), 'green')
xlabel('a')
ylabel('a prime')
xlim([ax_min, ax_max]);
ylim([ax_min, ax_max]);
legend({sprintf('r = %.4f', r_original),...
        sprintf('r = %.4f', r_vec(1)),...
        sprintf('r = %.4f', r_vec(2))}, 'Location','northwest')
title(sprintf('Policy function: l = %.2f', l_vec(2)))
hold off

subplot(3, 1, 3);
hold on
plot(a_baseline, a_prm_baseline(3,:), 'red')
plot(a_r_chg(:,1), a_prm_r_chg(3,:,1),'blue')
plot(a_r_chg(:,2), a_prm_r_chg(3,:,2), 'green')
xlabel('a')
ylabel('a prime')
xlim([ax_min, ax_max]);
ylim([ax_min, ax_max]);
legend({sprintf('r = %.4f', r_original),...
        sprintf('r = %.4f', r_vec(1)),...
        sprintf('r = %.4f', r_vec(2))}, 'Location','northwest')
title(sprintf('Policy function: l = %.2f', l_vec(3)))
hold off

print(gcf, 'Result_r_changes','-dpng') %'-dpng':PNG, '-dpdf':PDF
close(gcf)  


% 5-3) interest rate:
gcf = figure;
subplot(3, 1, 1)
hold on
plot(a_baseline, V_baseline(1, :), 'red')
plot(a_r_chg(:,1), V_r_chg(1,:,1), 'blue')
plot(a_r_chg(:,2), V_r_chg(1,:,2), 'green')
xlabel('a')
ylabel('V')
xlim([ax_min, ax_max]);
legend({sprintf('r = %.4f', r_original),...
        sprintf('r = %.4f', r_vec(1)),...
        sprintf('r = %.4f', r_vec(2))}, 'Location','northwest')
title(sprintf('Value function: l = %.2f', l_vec(1)))
hold off

subplot(3, 1, 2)
hold on
plot(a_baseline, V_baseline(2, :), 'red')
plot(a_r_chg(:,1), V_r_chg(2,:,1), 'blue')
plot(a_r_chg(:,2), V_r_chg(2,:,2), 'green')
xlabel('a')
ylabel('V')
xlim([ax_min, ax_max]);
legend({sprintf('r = %.4f', r_original),...
        sprintf('r = %.4f', r_vec(1)),...
        sprintf('r = %.4f', r_vec(2))}, 'Location','northwest')
title(sprintf('Value function: l = %.2f', l_vec(2)))
hold off

subplot(3, 1, 3);
hold on
plot(a_baseline, V_baseline(3,:), 'red')
plot(a_r_chg(:,1), V_r_chg(3,:,1),'blue')
plot(a_r_chg(:,2), V_r_chg(3,:,2), 'green')
xlabel('a')
ylabel('V')
xlim([ax_min, ax_max]);
legend({sprintf('r = %.4f', r_original),...
        sprintf('r = %.4f', r_vec(1)),...
        sprintf('r = %.4f', r_vec(2))}, 'Location','northwest')
title(sprintf('Value function: l = %.2f', l_vec(3)))
hold off

print(gcf, 'Result_r_changes_V','-dpng') %'-dpng':PNG, '-dpdf':PDF
close(gcf)  


% 5-4) Borrowing contraint:
gcf = figure;
subplot(3, 1, 1)
hold on
plot(a_a_lb_chg(:,1), a_prm_a_lb_chg(1,:,1), 'blue')
plot(a_baseline, a_prm_baseline(1, :), 'red')
plot(a_a_lb_chg(:,2), a_prm_a_lb_chg(1,:,2), 'green')
xlabel('a')
ylabel('a prime')
xlim([ax_min, ax_max]);
ylim([ax_min, ax_max]);
legend({sprintf('a lb = %.1f', a_lb_vec(1)),...
        sprintf('a lb = %.1f', r_original),...
        sprintf('a lb = %.1f', a_lb_vec(2))}, 'Location','northwest')
title(sprintf('Policy function: l = %.2f', l_vec(1)))
hold off

subplot(3, 1, 2)
hold on
plot(a_a_lb_chg(:,1), a_prm_a_lb_chg(2,:,1), 'blue')
plot(a_baseline, a_prm_baseline(2, :),'red')
plot(a_a_lb_chg(:,2), a_prm_a_lb_chg(2,:,2),'green')
xlabel('a')
ylabel('a prime')
xlim([ax_min, ax_max]);
ylim([ax_min, ax_max]);
legend({sprintf('a lb = %.1f', a_lb_vec(1)),...
        sprintf('a lb = %.1f', r_original),...
        sprintf('a lb = %.1f', a_lb_vec(2))}, 'Location','northwest')
title(sprintf('Policy function: l = %.2f', l_vec(2)))
hold off

subplot(3, 1, 3)
hold on
plot(a_a_lb_chg(:,1), a_prm_a_lb_chg(3,:,1), '-blue')
plot(a_baseline, a_prm_baseline(3, :), 'red')
plot(a_a_lb_chg(:,2), a_prm_a_lb_chg(3,:,2), '-green')
xlabel('a')
ylabel('a prime')
xlim([ax_min, ax_max]);
ylim([ax_min, ax_max]);
legend({sprintf('a lb = %.1f', a_lb_vec(1)),...
        sprintf('a lb = %.1f', r_original),...
        sprintf('a lb = %.1f', a_lb_vec(2))}, 'Location','northwest')
title(sprintf('Policy function: l = %.2f', l_vec(3)))
hold off

print(gcf, 'Result_a_lb_changes','-dpng') %'-dpng':PNG, '-dpdf':PDF
close(gcf)  


% 5-5) Borrowing contraint:
gcf = figure;
subplot(3, 1, 1)
hold on
plot(a_a_lb_chg(:,1), V_a_lb_chg(1,:,1), 'blue')
plot(a_baseline, V_baseline(1, :), 'red')
plot(a_a_lb_chg(:,2), V_a_lb_chg(1,:,2), 'green')
xlabel('a')
ylabel('V')
xlim([ax_min, ax_max]);
legend({sprintf('a lb = %.1f', a_lb_vec(1)),...
        sprintf('a lb = %.1f', r_original),...
        sprintf('a lb = %.1f', a_lb_vec(2))}, 'Location','northwest')
title(sprintf('Value function: l = %.2f', l_vec(1)))
hold off

subplot(3, 1, 2)
hold on
plot(a_a_lb_chg(:,1), V_a_lb_chg(2,:,1), 'blue')
plot(a_baseline, V_baseline(2, :),'red')
plot(a_a_lb_chg(:,2), V_a_lb_chg(2,:,2),'green')
xlabel('a')
ylabel('V')
xlim([ax_min, ax_max]);
legend({sprintf('a lb = %.1f', a_lb_vec(1)),...
        sprintf('a lb = %.1f', r_original),...
        sprintf('a lb = %.1f', a_lb_vec(2))}, 'Location','northwest')
title(sprintf('Value function: l = %.2f', l_vec(2)))
hold off

subplot(3, 1, 3)
hold on
plot(a_a_lb_chg(:,1), V_a_lb_chg(3,:,1), '-blue')
plot(a_baseline, V_baseline(3, :), 'red')
plot(a_a_lb_chg(:,2), V_a_lb_chg(3,:,2), '-green')
xlabel('a')
ylabel('V')
xlim([ax_min, ax_max]);
legend({sprintf('a lb = %.1f', a_lb_vec(1)),...
        sprintf('a lb = %.1f', r_original),...
        sprintf('a lb = %.1f', a_lb_vec(2))}, 'Location','northwest')
title(sprintf('Value function: l = %.2f', l_vec(3)))
hold off

print(gcf, 'Result_a_lb_changes_V','-dpng') %'-dpng':PNG, '-dpdf':PDF
close(gcf)  

save consumer_problem_result;

%% ================== Optimized function ==================================
%%  
function fval = ValFunc_Aiyagari(a_prm)

global  l_j r w beta Pi_j a_vec a_k V0_mat

% (1) Calculate continuous V(a') -----------------------------------------
% We prepare values of value function only on the grid points wrt a.
% However, MATLAB optimization routine can assign values not in our
% grid points to a'. So, we have to calcluate continuous V(a') from
% our discrete V(a)(=V(a')). Here, I do so by linear interporation. 

a_below_id = sum(a_prm > a_vec);

V_a_prm = V0_mat(:, a_below_id) ...
         + (V0_mat(:, a_below_id + 1 ) - V0_mat(:, a_below_id))...
         ./(a_vec(a_below_id + 1) - a_vec(a_below_id)) ...
         * (a_prm - a_vec(a_below_id));

V_a_prm(isnan(V_a_prm)) = -Inf;

% (2) Assign the function value -------------------------------------------
fval = log((1 + r) * a_k + w * l_j - a_prm) + beta * Pi_j * V_a_prm;

fval = -fval; % make it negative since MATLAB minimizes this function,
              % not maximize.
end