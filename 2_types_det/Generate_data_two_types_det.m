%code created by Loïc Marrec


clear all;
close all;


% Parameters

Nit = 1e4;                          % Number of replicates 
N_gen = 1e2;                        % Maximum number of generations of hosts for each replicate (assume no extinction if it is reached)
N_host = 1e2;                       % Number of contaminated individuals at the current generation of hosts beyond which we assume that no extinction happens
N_bs = 1e2;                         % Number of transmitted bacteria (bottleneck size and cluster size)
lambda_i = 2;                       % Mean number of infected hosts from an immune host
lambda_n = 2;                       % Mean number of infected hosts from a naive host
q_i = .55;                          % Probability of treatment for an immune host
q_n = .55;                          % Probability of treatment for a naive host
fs = 1;                             % Fitness of sensitive (S) bacteria
fr = .9;                            % Fitness of resistant (R) bacteria
mu1 = 1e-4/2/log(2);                % Probability of mutation S -> R
tau = 7;                            % Incubation time (such that G=tau/log(2) with G number of generations of bacteria within a host) 
frac = 1;                           % Fraction of immune hosts in the population
n = 1;                              % Initial number of resistant bacteria


% Construct the table of results - one result per replicate

endTest_list = NaN(1, Nit);	

% In endTest_list:
% 0 will mean no extinction because the number of generations of hosts reached N_gen
% 1 will mean extinction
% 2 will mean no extinction because the number of contaminated hosts at the current generation of hosts exceeded N_host


% Generate the data - call the function BP_two_types_det

[endTest_list(1, :)] = BP_two_types_det(Nit, N_gen, N_host, N_bs, lambda_i, lambda_n, q_i, q_n, frac, n, fs, fr, mu1, tau);                                        
