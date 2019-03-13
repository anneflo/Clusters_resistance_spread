%code created by Loïc Marrec and Anne-Florence Bitbol

function [endTest_list] = BP_two_types_stoch(Nit, N_gen, N_host, N_bs, lambda_i, lambda_n, q_i, q_n, frac, n, fs, fr, gs, gr, mu1, n_div_lim)

    endTest_list = NaN(Nit, 1);

    for Index = 1 : Nit

        ResI = zeros(1, N_bs+1);       % Table for immune individuals: ResI(i) is the number of contaminated immune individuals infected with i-1 resistant bacteria at the current generation of hosts
        ResN = zeros(1, N_bs+1);       % Table for naive individuals: ResN(i) is the number of contaminated naive individuals infected with i-1 resistant bacteria at the current generation of hosts
        endTest = 0;

        % Initialization: is the first individual naive or immune?

        gen = 0;                   % Initialization of the generation index
        i_or_n = rand;             % Is the first individual immune or naive?

        if i_or_n > frac           % The first individual is naive

            ResN(n+1) = 1;    % Let's update the table of naive individuals

        else                       % The first individual is immune

            ResI(n+1) = 1;    % Let's update the table of immune individuals

        end

        while gen <= N_gen-1 && endTest == 0

            gen = gen+1;    % We prepare the new generation "gen" of hosts

            ResI_new = zeros(1, N_bs+1);   % Table for immune individuals of the generation "gen" of hosts
            ResN_new = zeros(1, N_bs+1);   % Table for naive individuals of the generation "gen" of hosts

            % How many infected naive and infected immune individuals in the generation "gen-1" of hosts?  
            N_inf_i = sum(ResI);
            N_inf_n = sum(ResN);

            % N_R is the number of R bacteria
            % N_ind is the number of individuals who have N_R R bacteria
            [~, N_R_i, N_ind_i] = find(ResI);
            [~, N_R_n, N_ind_n] = find(ResN);
            N_R_i = N_R_i-1;
            N_R_n = N_R_n-1;

            % Let us start with the naive contaminating individuals
            % Go through the loop if and only if there is at least one naive infected individual in the generation "gen-1" of hosts
            if N_inf_n ~= 0

                for i = 1 : length(N_ind_n)

                    for j = 1 : N_ind_n(i)

                        t_or_nt = rand;     % Is the individual treated?

                        if t_or_nt <= q_n       % Yes, the individual is treated

                            if not(N_R_n(i) == 0)      % He/she has at least one R bacterium

                                k_n = poissrnd(lambda_n);       % How many individuals is he/she going to contaminate?

                                if k_n ~= 0     % He/she contaminates at least one individual
                                    
                                    if k_n > N_bs       % An individual cannot contaminate more than N_bs individuals 

                                        k_n = N_bs;

                                    end

                                    for q = 1 : k_n

                                        i_or_n = rand;      % Is the contaminated individual is naive or immune?

                                        if i_or_n > frac % The contaminated individual is naive

                                            ResN_new(N_bs+1) = ResN_new(N_bs+1)+1; % Let's update the table of naive individuals (only resistant bacteria here since treatment happened)

                                        else % The contaminated individual is immune

                                            ResI_new(N_bs+1) = ResI_new(N_bs+1)+1; % Let's update the table of immune individuals (only resistant bacteria here since treatment happened)

                                        end

                                    end

                                end

                            end

                        else        % No, the individual is not treated

                            [~, ~, ~, r_n, ~, ~] = Gillespie_two_types_stoch(fs, gs, N_bs-N_R_n(i), fr, gr, N_R_n(i), mu1, n_div_lim); % Compute the proportion of R bacteria in the naive individual

                            if r_n < 0 % If numerical accuracy problem

                                r_n = 0;

                            end

                            N_R_growth_n = binornd(N_bs^2, r_n); % Compute the number of R bacteria after time tau among the N_bs^2 bacteria in the naive individual

                            k_n = poissrnd(lambda_n);       % How many individuals is he going to contaminate?

                            if k_n ~= 0     % He/she contaminates at least one individual
                                
                                if k_n > N_bs % An individual cannot contaminate more than N_bs individuals

                                    k_n = N_bs;

                                end

                                p = randperm(N_bs^2, N_bs*k_n); % Prepare k_n packets of N_bs bacteria from the N_bs^2 bacteria
                                start = 1;

                                for q = 1 : k_n

                                    i_or_n = rand;      % Is the contaminated individual is naive or immune?
                                    N_R_trans_n = size(p(p(start : start+N_bs-1) <= N_R_growth_n), 2); % How many R bacteria does the contaminated individual receive?

                                    if i_or_n > frac % The contaminated individual is naive

                                        ResN_new(N_R_trans_n+1) = ResN_new(N_R_trans_n+1)+1; % Let's update the table of naive individuals

                                    else % The contaminated individual is immune

                                        ResI_new(N_R_trans_n+1) = ResI_new(N_R_trans_n+1)+1; % Let's update the table of immune individuals

                                    end

                                    start = start+N_bs;

                                end

                            end

                        end 

                    end

                end

            end

            % Let us consider the immune contaminating ind
            % Go through the loop if and only if there is at least one immune infected individual in the generation "gen-1" of hosts
            if N_inf_i ~= 0

                for i = 1 : length(N_ind_i)

                    for j = 1 : N_ind_i(i)

                        t_or_nt = rand; % Is the individual treated?

                        if t_or_nt <= q_i       % Yes, the individual is treated

                            if not(N_R_i(i) == 0)          % He/she has at least one R bacterium

                                k_i = poissrnd(lambda_i); % How many individuals is he/she going to contaminate?

                                if k_i ~= 0 % He/she contaminates at least one individual
                                    
                                    if k_i > N_bs % An individual cannot contaminate more than N_bs individuals

                                        k_i = N_bs;

                                    end

                                    for q = 1 : k_i

                                        i_or_n = rand;      % Is the contaminated ind naive or immune?

                                        if i_or_n > frac        % The contaminated individual is naive

                                            ResN_new(N_bs+1) = ResN_new(N_bs+1)+1; % Let's update the table of naive individuals

                                        else % The contaminated individual is immune

                                            ResI_new(N_bs+1) = ResI_new(N_bs+1)+1; % Let's update the table of immune individuals

                                        end

                                    end

                                end

                            end

                        else        % No, the individual is not treated

                            [~, ~, ~, r_i, ~, ~] = Gillespie_two_types_stoch(fs, gs, N_bs-N_R_i(i), fr, gr, N_R_i(i), mu1, n_div_lim); % Compute the proportion of R bacteria after time tau in the immune individual

                             if r_i < 0 % If numerical accuracy problem

                                 r_i = 0;

                             end

                             N_R_growth_i = binornd(N_bs, r_i); % Compute the number of R bacteria after time tau among the N_bs^2 bacteria in the immune individual                                            

                             k_i = poissrnd(lambda_i); % How many individuals is he/she going to contaminate?

                             if k_i ~= 0 % He/she contaminates at least one individual
                                 
                                 if k_i > N_bs % An individual cannot contaminate more than N_bs individuals

                                    k_i = N_bs;

                                 end

                                 p = randperm(N_bs, k_i); % Prepare k_i clonal packets

                                 for q = 1 : k_i

                                     if p(q) <= N_R_growth_i % It is a packet of R bacteria

                                         N_R_trans_i = N_bs; 

                                     else % It is a packet of S bacteria

                                         N_R_trans_i = 0;

                                     end

                                     i_or_n = rand;      % Is the contaminated individual naive or immune?

                                     if i_or_n > frac        % The contaminated individual is naive

                                         ResN_new(N_R_trans_i+1) = ResN_new(N_R_trans_i+1)+1;

                                     else % The contaminated individual is immune
 
                                         ResI_new(N_R_trans_i+1) = ResI_new(N_R_trans_i+1)+1;

                                     end

                                 end

                             end

                        end

                    end

                end

            end

            ResI = ResI_new; % Let's update the table for immune individuals
            ResN = ResN_new; % Let's update the table for naive individuals

            % Check if there is extinction
            if  sum(ResI)+sum(ResN) == 0

                endTest = 1;

            end

            % Assume that the epidemic won't go extinct if the number of contaminated individuals at the generation "gen" of hosts is larger than N_host
            if  sum(ResI)+sum(ResN) > N_host

                endTest = 2;

            end

        end

        endTest_list(Index, 1) = endTest;

    end
    
end
