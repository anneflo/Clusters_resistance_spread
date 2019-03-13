%code created by Loïc Marrec

function [Prop] = Compute_prop_for_BP_det(fs, fr, mu1, n, N, tau)

    T = [fs*(1-mu1) 0; fs*mu1 fr]; % Matrix involved in the system of differential equations describing the within-host growth

    [V, D] = eig(T); % Compute the eigenvectors and the eigenvalues of the matrix T

    s0 = 1-n/N;	% Initial proportion of S bacteria
    r0 = n/N; % Initial proportion of R bacteria

    A = [V(1, 1)-s0*(V(1, 1)+V(2, 1)) V(1, 2)-s0*(V(1, 2)+V(2, 2));...
         V(2, 1)-r0*(V(1, 1)+V(2, 1)) V(2, 2)-r0*(V(1, 2)+V(2, 2))];

    B = null(A);

    alpha = B(1, 1);
    beta = B(2, 1);
    
    Prop = 1-(V(1, 1)*alpha*exp(D(1, 1)*tau)+V(1, 2)*beta*exp(D(2, 2)*tau))/...
           (alpha*(V(1, 1)+V(2, 1))*exp(D(1, 1)*tau)+beta*(V(1, 2)+V(2, 2))*exp(D(2, 2)*tau)); % Propotion of R bacteria
    
end


