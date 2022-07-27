%% mt_Econ501B_Final_Q6.m
%  is the MATLAB program for Question 6 in the final exam of Economics
%  501B.
%
% ------------------------------------------------------------------------
% Copyright (c) 2021 Masaki Tanaka (Washington University in St Louis)
% All rights researved.
%
% This source code or any portion thereof must not be reproduced or used
% in any manner whatsoever.
%
% Create: Dec 13, 2021
%

%% 0. Preamble ***********************************************************
close all;
clear;
clc;

%% 1. Set parameters *****************************************************
% 1-1) model parameters
beta  = 0.94;
sig   = 2.00;
a0    = 3.00;
b     = 2.00;
w_min = 1.00;
w_max = 3.00;
p     = 0.50;

% 1-2) simulation paramters
nGrids   = 1000      ; % # of grids for w
max_iter = 500       ; % maximum # of iterations
tol      = 10^(-20)  ; % convergence criteria
a_min    = 0.00      ;
a_max    = 5.00      ;

% 1-3) Graphics option
a2plot   = [0, 1, 2, 3]; % values of a, for which V will be illustrated
colors   = {'red', 'blue', 'green', 'black'};

% 1-4) Prepare grids for w and a, and calculate P(grid(i)):
%      automatically computed
w_vec = linspace(w_min, w_max, nGrids); % Grids for w
Pw    = 1 /nGrids; % prob. that each grid is drawn

a_vec = linspace(a_min, a_max, nGrids)'; % Grids for a

% 1-5) calculate r based on beta: automatically computed
r = 1/beta - 1;

%% 2. Run simluation *****************************************************
% initial guess of value function (V_0: b is expired, V_b: b > 0)
V_0  = zeros(length(a_vec), nGrids);
V_b  = zeros(length(a_vec), nGrids);

i    = 0; % iteration counter
diff = 1; % improvement from the previous guess

% 2-1) value of accepting offer
%     trick: In MATLAB, "horizontal vector plus vertical vector"
%            yields a matrix whose i-j entry is sum of j-th element
%            of vertical vector and i-th element horizontal one.
VA = (1 + r)/r * ((r * a_vec + w_vec).^(1 - sig))/ (1 - sig);

aA = repmat(a_vec, 1, size(VA, 2));  % best a(t+1) when employed


while (i < max_iter) & (diff > tol)
    % ******************* [The case b_t = 0] *************************
    % 2-2) value of rejecting offer
    %      --- horizontal: a(t+1), vertical : a(t)
    c_mat   = (1 + r)* a_vec + 0 - a_vec';
    VR_0_cand = (c_mat).^(1 - sig) / (1 - sig);
    VR_0_cand = VR_0_cand + 1/(1 + r) * Pw * sum(V_0, 2)';

    % if consumption is negative, give penalty
    VR_0_cand(c_mat < 0) = - Inf;

    % taking max w.r.t. the horizontal direction (a(t+1))
    [VR_0, argVR_0] = max(VR_0_cand, [], 2);
    VR_0    = repmat(VR_0   , 1, size(VA, 2)); % adjust the sizes of VA & VR

    aR_0 = a_vec(argVR_0); 
    aR_0 = repmat(aR_0, 1, size(VA, 2)); % best a(t+1) when unemployed

    % 2-3) taking max between accept and reject
    V_0_new = max(VA, VR_0);

    % ******************* [The case b_t = b] *************************
    % 2-4) value of rejecting offer
    %      --- horizontal: a(t+1), vertical : a(t)
    c_mat   = (1 + r)* a_vec + b - a_vec';
    VR_b_cand = (c_mat).^(1 - sig) / (1 - sig);
    VR_b_cand = VR_b_cand + 1/(1 + r) * ...
              (p * Pw * sum(V_0, 2)' + (1 - p) * Pw * sum(V_b, 2)');


    % if consumption is negative, give penalty
    VR_b_cand(c_mat < 0) = - Inf;

    % taking max w.r.t. the horizontal direction (a(t+1))
    [VR_b, argVR_b] = max(VR_b_cand, [], 2);
    VR_b    = repmat(VR_b   , 1, size(VA, 2)); % adjust the sizes of VA & VR

    aR_b = a_vec(argVR_b); 
    aR_b = repmat(aR_b, 1, size(VA, 2)); % best a(t+1) when unemployed

    % 2-5) taking max between accept and reject
    V_b_new = max(VA, VR_b);

    % ***************************************************************
    % 2-6) ending i-th iteration  
    diff_0 = max(abs(V_0_new - V_0), [], 'all'); % improvement in V_0
    diff_b = max(abs(V_b_new - V_b), [], 'all'); % improvement in V_b
    diff   = max(diff_0, diff_b);

    V_0  = V_0_new; % Use the obtained V_0 as the next guess
    V_b  = V_b_new; % Use the obtained V_b as the next guess
    
    i    = i + 1; % Move forward the iteration counter

    if  mod(i, 100) == 0
        fprintf('\n Iteration: %d...', i);
    end
end

if  diff < tol
    fprintf('\n ********************************************** \n')
    fprintf(' Iteration stopped due to convergent criteria')
    fprintf('\n ********************************************** \n')
else
    fprintf('\n ********************************************** \n')
    fprintf(' Iteration stopped due to maxiter')
    fprintf('\n ********************************************** \n')
end

% 2-5) Find the reservation wage
w_r_0_id  = sum((V_0 == V_0(: , 1)), 2);
w_r_0_vec = w_vec(w_r_0_id);

w_r_b_id  = sum((V_b == V_b(: , 1)), 2);
w_r_b_vec = w_vec(w_r_b_id);

% 2-6) Find the policy function of a(t+1)
optimal_a_t1_0 = aA;
optimal_a_t1_0(V_0 == VR_0) = aR_0(V_0 == VR_0);

optimal_a_t1_b = aA;
optimal_a_t1_b(V_b == VR_b) = aR_b(V_b == VR_b);

%% 3. Graphics ***********************************************************

% ******************* [The case b_t = 0] *************************
f_0 = figure;

% 3-1) Value function
a2plot_id = sum(a_vec <= a2plot);

subplot(3, 1, 1)
hold on
for j = 1 : length(a2plot_id)
    plot(w_vec, V_0(a2plot_id(j), :), 'Color', colors{j},...
                'Linewidth', 1.5);
end
title('Value function: b_t = 0')
xlabel('w')
ylabel('V(a_t, 0, w)')
j = 1;
code = 'legend({';
while j < length(a2plot_id) + 1
    code = [code, sprintf(' ''a_t = %d'' ',a2plot(j))];
    if j < length(a2plot_id)
        code = [code, ','];
    else
        code = [code, '}, ''Location'', ''southeast'');'];
    end
    j = j + 1;
end
eval(code);
hold off

% 3-2) Policy function
subplot(3, 1, 2)
hold on
for j = 1 : length(a2plot_id)
    plot(w_vec, optimal_a_t1_0(a2plot_id(j), :), 'Color', colors{j},...
        'Linewidth', 1.5);
end
title('Policy function: b_t = 0')
xlabel('w')
ylabel('a_{t+1}')
j = 1;
code = 'legend({';
while j < length(a2plot_id) + 1
    code = [code, sprintf(' ''a_t = %d'' ',a2plot(j))];
    if j < length(a2plot_id)
        code = [code, ','];
    else
        code = [code, '}, ''Location'', ''northwest'');'];
    end
    j = j + 1;
end
eval(code);
hold off


% 3-3) reservation wage
subplot(3, 1, 3)
hold on
plot(a_vec, w_r_0_vec, 'Color', 'red', 'Linewidth', 1.5);
title('Reservation wage: b_t = 0')
xlabel('a_t')
ylabel('reservation wage')
xlim([0,3])
hold off

print(f_0, 'Result_Q6_1','-dpng') %'-dpng':PNG, '-dpdf':PDF


% ******************* [The case b_t = b] *************************

f_b = figure;

% 3-1) Value function

subplot(3, 1, 1)
hold on
for j = 1 : length(a2plot_id)
    plot(w_vec, V_b(a2plot_id(j), :), 'Color', colors{j},...
                'Linewidth', 1.5);
end
title('Value function: b_t = b')
xlabel('w')
ylabel('V(a_t, b, w)')
j = 1;
code = 'legend({';
while j < length(a2plot_id) + 1
    code = [code, sprintf(' ''a_t = %d'' ',a2plot(j))];
    if j < length(a2plot_id)
        code = [code, ','];
    else
        code = [code, '}, ''Location'', ''southeast'');'];
    end
    j = j + 1;
end
eval(code);
hold off

% 3-2) Policy function
subplot(3, 1, 2)
hold on
for j = 1 : length(a2plot_id)
    plot(w_vec, optimal_a_t1_b(a2plot_id(j), :), 'Color', colors{j},...
        'Linewidth', 1.5);
end
title('Policy function: b_t = b')
xlabel('w')
ylabel('a_{t+1}')
j = 1;
code = 'legend({';
while j < length(a2plot_id) + 1
    code = [code, sprintf(' ''a_t = %d'' ',a2plot(j))];
    if j < length(a2plot_id)
        code = [code, ','];
    else
        code = [code, '}, ''Location'', ''northwest'');'];
    end
    j = j + 1;
end
eval(code);
hold off


% 3-3) reservation wage
subplot(3, 1, 3)
hold on
plot(a_vec, w_r_b_vec, 'Color', 'red', 'Linewidth', 1.5);
title('Reservation wage: b_t = b')
xlabel('a_t')
ylabel('reservation wage')
xlim([0,3])
hold off

print(f_b, 'Result_Q6_2','-dpng') %'-dpng':PNG, '-dpdf':PDF