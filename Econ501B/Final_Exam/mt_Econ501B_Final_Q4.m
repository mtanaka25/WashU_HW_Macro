%% mt_Econ501B_Final_Q4.m
%  is the MATLAB program for Question 4 in the final exam of Economics
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
b     = 0.70;
w_min = 1.00;
w_max = 3.00;

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
V    = zeros(length(a_vec), nGrids); % initial guess of value function.
i    = 0; % iteration counter
diff = 1; % improvement from the previous guess

% 2-1) value of accepting offer
%     trick: In MATLAB, "horizontal vector plus vertical vector"
%            yields a matrix whose i-j entry is sum of j-th element
%            of vertical vector and i-th element horizontal one.
VA = (1 + r)/r * ((r * a_vec + w_vec).^(1 - sig))/ (1 - sig);

aA = repmat(a_vec, 1, size(VA, 2));  % best a(t+1) when employed


while (i < max_iter) & (diff > tol)
    % 2-2) value of rejecting offer
    %      --- horizontal: a(t+1), vertical : a(t)
    c_mat   = (1 + r)* a_vec + b - a_vec';
    VR_cand = (c_mat).^(1 - sig) / (1 - sig);
    VR_cand = VR_cand + 1/(1 + r) * Pw * sum(V, 2)';


    % if consumption is negative, give penalty
    VR_cand(c_mat < 0) = - Inf;

    % taking max w.r.t. the horizontal direction (a(t+1))
    [VR, argVR] = max(VR_cand, [], 2);
    VR    = repmat(VR   , 1, size(VA, 2)); % adjust the sizes of VA & VR

    aR = a_vec(argVR); 
    aR = repmat(aR, 1, size(VA, 2)); % best a(t+1) when unemployed

    % 2-3) taking max between accept and reject
    V_new = max(VA, VR);

    % 2-4) ending i-th iteration  
    diff = max(abs(V_new - V), [], 'all'); % improvement in V
    V    = V_new; % Use the obtained V as the next guess
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
w_r_id  = sum((V == V(: , 1)), 2);
w_r_vec = w_vec(w_r_id);

% 2-6) Find the policy function of a(t+1)
optimal_a_t1 = aA;
optimal_a_t1(V == VR) = aR(V == VR);

%% 3. Graphics ***********************************************************
f = figure;

% 3-1) Value function
a2plot_id = sum(a_vec <= a2plot);

subplot(3, 1, 1)
hold on
for j = 1 : length(a2plot_id)
    plot(w_vec, V(a2plot_id(j), :), 'Color', colors{j},...
                'Linewidth', 1.5);
end
title('Value function')
xlabel('w')
ylabel('V(a_t, w)')
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
    plot(w_vec, optimal_a_t1(a2plot_id(j), :), 'Color', colors{j},...
        'Linewidth', 1.5);
end
title('Policy function')
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
plot(a_vec, w_r_vec, 'Color', 'red', 'Linewidth', 1.5);
title('Reservation wage')
xlabel('a_t')
ylabel('reservation wage')
xlim([0,3])
hold off

print(f, 'Result_Q4','-dpng') %'-dpng':PNG, '-dpdf':PDF