%code created by Loïc Marrec

function [XA_f, XB_f, Prop_A, Prop_B, t, n] = Gillespie_two_types_exp(fA, gA, XA_i, fB, gB, XB_i, mu1, n_div_lim)

        % Initialization

        XA = XA_i;              % Initialization of the number of A individuals
        XB = XB_i;              % Initialization of the number of B individuals
        t = 0;                  % Initialization of time
        cumul = zeros(5, 1);    % To build the sampling tower
        n = 0;                  % Initialization of the number of iterations
        n_div = 0;              % Initialization of the number of divisions

        while XA ~= 0 && n_div <= n_div_lim	% The process stops if the B type fixes or if the number of divisions is reached
            
            % Compute the transition rates

            T_A_rep_no_mut = fA*(1-mu1)*XA;	% Reproduction without mutation rate of the A type   
            T_A_rep_mut = fA*mu1*XA;		% Reproduction with mutation rate of the A type  
            T_A_death = gA*XA;			% Death rate of the A type			
            T_B_rep = fB*XB;			% Reproduction rate of the B type			
            T_B_death = gB*XB;			% Death rate of the B type			
            T = T_A_rep_no_mut+T_A_rep_mut+T_A_death+T_B_rep+T_B_death;		% Sum of all the transition rates

            % Increase the time t by the randomly generated time tau which is exponentially distributed with mean 1/T

            r1 = rand;
            tau = 1/T*log(1/r1);
            t = t+tau;

            % Build a sampling tower

            ir2 = 1;
            r2 = rand;
            cumul(1) = T_A_rep_no_mut;
            cumul(2) = T_A_rep_no_mut+T_A_rep_mut;
            cumul(3) = T_A_rep_no_mut+T_A_rep_mut+T_A_death;
            cumul(4) = T_A_rep_no_mut+T_A_rep_mut+T_A_death+T_B_rep;
            cumul(5) = T;

            % Determine which event occurs and update the number of the different types of bacteria

            while cumul(ir2) < r2*T

                 ir2 = ir2+1;

             end

             if ir2 == 1	% Reproduction of one A bacterium without mutation

                 XA = XA+1;
                 n_div = n_div+1;

             elseif ir2 == 2	% Reproduction of one A bacterium with mutation

                 XB = XB+1;
                 n_div = n_div+1;

             elseif ir2 == 3	% Death of one A bacterium

                 XA = XA-1;

             elseif ir2 == 4	% Reproduction of one B bacterium 

                 XB = XB+1;
                 n_div = n_div+1;

             elseif ir2 == 5	% Death of one B bacterium

                 XB = XB-1;

             end
             
             n = n+1;	% Update the number of iterations
            
        end
        
        XA_f = XA;	% The final number of A bacteria
        XB_f = XB;	% The final number of B bacteria
        Prop_A = XA/(XA+XB);	% The final proportion of A bacteria
        Prop_B = XB/(XA+XB);	% The final proportion of B bacteria

end
