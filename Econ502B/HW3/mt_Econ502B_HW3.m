%% mt_Econ502B_HW3.m
% the MATLAB program to solve Question 3 of HW3 of Econ 502B (
% Washington University in St. Louis, 2022SP).
%
% ------------------------------------------------------------------------
% Copyright (c) 2022 Masaki Tanaka (Washington University in St Louis)
% All rights researved.
%
% This source code or any portion thereof must not be reproduced or used
% in any manner whatsoever.
%
% Create: Apr 26, 2022
%

%% 0. Housekeeping *******************************************************
close all;
clear;
clc;


%% 1. Set Parameters *****************************************************

% 1-1) Structural parameters

theta  = 0.5;   % concavity of production function
lmbd   = 0.1;   % prob. of the shift from sL to sH
sL     = 1  ;   % productivity: Low
sH     = 2  ;   % productivity: High
cf     = 0.1;   % fixed cost  
ce     = 0.1;   % entry cost
r      = 0.02;  % interest rate
g_pre  = 0.02;  % population growth rate: pre-shock
g_post = 0.012; % population growth rate: post-shock

% 1-2) Setting for the simulation
Tmax = 50; % # of periods to be simulated
colors   = {'red', 'blue', 'green', 'black'};

% 1-3) Define auxiliary parameters (for readability)
sgm    = 1/(1 - theta); 
sL_sgm = sL^(sgm);
sH_sgm = sH^(sgm);


%% 2. Comupute Stationary Demography of Firms ****************************
% 2-1) Pre-shock stationary distribution =================================

% Steady-state price: pre-shock
z_pre = 1/theta * (...
    (r + lmbd)*( 1/(1+r)*ce + 1/r*cf )...
    /...
    ( ( sL_sgm + lmbd/r*sH_sgm)*(1/theta - 1) )...
    )^(1 - theta);

% Steady-state avarage firm size: pre-shock
e_pre = (theta*z_pre)^sgm *...
        (lmbd * sH_sgm + g_pre * sL_sgm)/(lmbd + g_pre) +...
         cf + g_pre/(1+g_pre)*ce;

% Stationary distibution of firms:: pre-shock
M_pre    = 1/e_pre;
M_sH_pre = lmbd  / (lmbd + g_pre) * M_pre;
M_sL_pre = g_pre / (lmbd + g_pre) * M_pre;
m_pre    = g_pre / (  1  + g_pre) * M_pre;


% 2-2) Post-shock stationary distribution ================================

% Steady-state price: post-shock
z_post = z_pre; 

% Steady-state avarage firm size: post-shock
e_post = (theta*z_post)^sgm *...
         (lmbd * sH_sgm + g_post * sL_sgm)/(lmbd + g_post) +...
         cf + g_post/(1+g_post)*ce;

% Stationary distibution of firms:: post-shock
M_post    = 1/e_post;
M_sH_post = lmbd   / (lmbd + g_post) * M_post;
M_sL_post = g_post / (lmbd + g_post) * M_post;
m_post    = g_post / (  1  + g_post) * M_post;


% 2-3) Compare stationary demography ====================================
% This part is not manadatory for the assigment. Just for reference.

data2plt = [M_pre  , M_sH_pre , M_sL_pre , m_pre;
            M_post , M_sH_post, M_sL_post, m_post]';

fig1 = figure;
hold on
x = categorical({'M', 'M(s^H)', 'M(s^L)', 'm'});
x = reordercats(x, {'M', 'M(s^H)', 'M(s^L)', 'm'});
bar(x', data2plt);
title('Stationary values (per worker)');
legend('Pre-shock', 'Post-shock');
hold off
print(fig1, 'SS_Dist','-dpng') 

%% 3. Solve Transitional Dynamics ***************************************

% Prepare vectors to store the results
M_vec      = zeros(Tmax, 1); % # of firms
M_sH_vec   = zeros(Tmax, 1); % # of high-tech
M_sL_vec   = zeros(Tmax, 1); % # of low-tech
m_vec      = zeros(Tmax, 1); % # of entrants
m_rate_vec = zeros(Tmax, 1); % entry rates
e_vec      = zeros(Tmax, 1); % average firm size
Y_vec      = zeros(Tmax, 1); % real GDP per worker

% Price does not change on the transition
z_vec      = z_post * ones(Tmax, 1); 

% Set initial values (the economy starts from the pre-shock s.s.)
M_0    = M_pre;
M_sL_0 = M_sL_pre;
M_sH_0 = M_sH_pre;

%define an auxiliary parameter (which often appears)
theta_z_sgm = (theta * z_post)^sgm;


for t = 1 : Tmax
    % [notation]
    % "_1" means a value at t
    % "_0" means a value at t-1 

    % Calculating m_1 ----------------------------------------------------
    % labor assigned for paying incumbents' fixed cost:
    contribA = cf * M_0/(1 + g_post);

    % labor assigned for production of firms who was low-tech in t-1.
    contribB = theta_z_sgm * ...
        ( lmbd*sH_sgm + (1 - lmbd)*sL_sgm ) * M_sL_0/(1 + g_post);

    % labor assigned for production of firms who was high-tech in t-1.
    contribC = theta_z_sgm * ...
        sH_sgm * M_sH_0/(1 + g_post);


    % # of entrants implied by labor market clearing.
    m_1 = (1 - contribA - contribB - contribC) ...
        / ((theta * z_post * sL)^sgm + cf + ce);


    % Based on m_1, calculate the firm demography at t -------------------
    
    % low-tech
    M_sL_1 = (1 - lmbd) * M_sL_0 /(1 + g_post) + m_1;

    % high-tech
    M_sH_1 = M_sH_0 /(1 + g_post) + lmbd * M_sL_0 /(1 + g_post);

    % total # of firms
    M_1 = M_sH_1 + M_sL_1;


    % Compute the rest of stats ------------------------------------------
    % entry rate
    m_rate_1 = m_1 / M_0;
    
    % average firm size
    e_1 = 1 / M_1;
   
    % real GDP per worker
    % --- production by one high-tech firm
    y_sH_1 = sH * (theta * z_post * sH)^(theta/(1-theta));
    % --- production by one low-tech firm
    y_sL_1 = sL * (theta * z_post * sL)^(theta/(1-theta));
    % --- aggregation
    Y_1 = y_sH_1 * M_sH_1 + y_sL_1 * M_sL_1;

    % Store values -------------------------------------------------------
    M_vec(t)      = M_1;
    M_sL_vec(t)   = M_sL_1;
    M_sH_vec(t)   = M_sH_1;
    m_vec(t)      = m_1;
    m_rate_vec(t) = m_rate_1;
    e_vec(t)      = e_1;
    Y_vec(t)      = Y_1;

    % Set the initial values for the next loop ---------------------------
    M_0    = M_1;
    M_sL_0 = M_sL_1;
    M_sH_0 = M_sH_1;
end


%% 4. Graphics  **********************************************************

% page 1 =================================================================
fig2 = figure;
subplot(3, 1, 1)
hold on
plot(M_vec, 'Color', colors{1},'Linewidth', 1.5);
title('# of firms per worker')
ylabel('M_t/N_t')
xlabel('elapsed periods')

subplot(3, 1, 2)
plot(M_sH_vec, 'Color', colors{1},'Linewidth', 1.5);
title('# of high-tech firms per worker')
ylabel('M_t(s^H)/N_t')
xlabel('elapsed periods')

subplot(3, 1, 3)
plot(M_sL_vec, 'Color', colors{1},'Linewidth', 1.5);
title('# of low-tech firms per worker')
ylabel('M_t(s^L)/N_t')
xlabel('elapsed periods')
hold off

print(fig2, 'Result1','-dpng') 

% page 2 =================================================================

fig3 = figure;
subplot(3, 1, 1)
hold on
plot(m_vec, 'Color', colors{1},'Linewidth', 1.5);
title('# of entries per worker')
ylabel('m_t/N_t')
xlabel('elapsed periods')

subplot(3, 1, 2)
plot(m_rate_vec, 'Color', colors{1},'Linewidth', 1.5);
title('entry rate')
ylabel('m_t/M_{t-1}')
xlabel('elapsed periods')

subplot(3, 1, 3)
plot(e_vec, 'Color', colors{1},'Linewidth', 1.5);
title('average firm size')
ylabel('e_t')
xlabel('elapsed periods')
hold off

print(fig3, 'Result2','-dpng') 

% page 3 =================================================================

fig4 = figure;
subplot(3, 1, 1)
hold on
plot(Y_vec, 'Color', colors{1},'Linewidth', 1.5);
title('real GDP per worker')
ylabel('Y_t/N_t')
xlabel('elapsed periods')

subplot(3, 1, 2)
plot(z_vec, 'Color', colors{1},'Linewidth', 1.5);
title('price')
ylabel('z_t')
xlabel('elapsed periods')

print(fig4, 'Result3','-dpng') 


%% Postamble **************************************************************
fprintf('*********************\n')
fprintf('      finished       \n')
fprintf('*********************\n')