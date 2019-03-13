%code created by Loïc Marrec


clear all;
close all;


% Parameters

Nit = 1e4;                          % Number of replicates 
N_gen = 1e2;                        % Maximum number of generations of hosts for each replicate (assume no extinction if it is reached)
N_host = 1e2;                       % Number of contaminated individuals at the current generation of hosts beyond which we assume that no extinction happens
N_bs = 1e2;                         % Number of transmitted bacteria (bottleneck size and cluster size)
lambda_i = 2;                       % Mean number of infected host from an immune host
lambda_n = 2;                       % Mean number of infected host from a naive host
q_i = 0.55;                         % Probability of treating a immune host
q_n = 0.55;                         % Probability of treating a naive host
fs = 1;                             % Fitness of sensitive (S) bacteria
fr = .9;                            % Fitness of resistant (R) bacteria
gs = 0;                             % Death rate of sensitive (S) bacteria
gr = 0;                             % Death rate of resistant (R) bacteria
mu1 = 1e-4;                         % Probability of mutation S -> R
n_div_lim = N_bs*1024;              % Growth time in number of divisions
frac = 1;                           % Fraction of immune hosts in the population
n = 1;                              % Initial number of resistant bacteria


% Construct the table of results - one result per replicate

endTest_list = NaN(1, Nit);

% In endTest_list:
% 0 will mean no extinction because the number of generations of hosts reached N_gen
% 1 will mean extinction
% 2 will mean no extinction because the number of contaminated hosts at the current generation of hosts exceeded N_host


% Generate the data - call the function BP_two_types_stoch

[endTest_list(1, :)] = BP_two_types_stoch(Nit, N_gen, N_host, N_bs, lambda_i, lambda_n, q_i, q_n, frac, n, fs, fr, gs, gr, mu1, n_div_lim);
