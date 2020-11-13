function [S,A,L,P] = solve_ode(t,x0,phivec,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type)
% Function that solves the ODEs describing S, A, L and P dynamics
% Input:
% t         time vector
% x0        Input vector with initial values for S, A, vec(L), vec(P).
% phivec    Vector with phi-values of the different strains
% r         Growth rate bacteria
% K         Scaled carrying capacity of bacteria
% B         Effective burst size
% alpha     Reactivation rate lysogens
% deltaP    Decay rate phages
% a         Scaled adsorption rate phages to cells
% cL        Production rate arbitrium by lysogens
% deltaA    Decay rate arbitrium
% u         Scaled uptake and breakdown rate arbitrium by cells
% mu        Mutation rate (prob of mutation in strain j)
% mut_type  Mutation type ("step" or "rand")

% Number of strains
ns = length(phivec);

% Check initial condition vector
if length(x0) ~= 2*ns+2
    disp('Error solve_ode: initial vector should have length 2ns + 2')
    return
end

% Check mut_type input
if (mut_type ~= "step" && mut_type ~= "rand")
    disp('Error solve_ode: mut_type should be step or rand')
    return
end

% Define mutation matrix
if ns == 1
    M = 1;
else
    if mut_type == 'step'       % Only include mutations to one strain up or down
        m_low = 0.5*mu*ones(1,ns-1);
        m_up = 0.5*mu*ones(1,ns-1);
        m_diag = [1-0.5*mu, (1-mu)*ones(1,ns-2), 1-0.5*mu];
        M = diag(m_diag) + diag(m_low,-1) + diag(m_up,1);
    else                        % Include mutations between all strains
        M = mu/(ns-1)*ones(ns);
        M = M - diag(diag(M)) + (1-mu)*eye(ns);
    end
end
       

% Solve the ODE
options1 = odeset('NonNegative',1:length(x0));  % Enforce positive solutions only
[t,x] = ode45(@RHS,t,x0,options1);

S = (x(:,1))';
A = (x(:,2))';
L = (x(:,3:(ns+2)))';
P = (x(:,(ns+3):end))';

    function fx = RHS(t,x)      % Function solved by ode45
       fx = zeros(2*ns+2,1);    % Vector to describe the outcome
       
       % Total number of cells
       N = x(1) + sum(x(3:(ns+2)));
       
       % Susceptible dynamics
       fx(1) = r*x(1)*(1-N/K) - B*a*x(1)*sum(x((ns+3):end));
       
       % Arbitrium dynamics
       fx(2) = B*a*x(1)*sum(x((ns+3):end)) + cL*sum(x(3:(ns+2))) - x(2)*(deltaA + u*N);
       
       % Lysogen and phage dynamics
       
       % Construction of full matrix Q that captures all effects of
       % lysogens and phages on the lysogen and phage dynamics
       Q = zeros(2*ns,2*ns);
       
       % Lysogens
       Q(1:ns,1:ns) = (r*(1-(x(1) + sum(x(3:(ns+2))))/K)-alpha)*eye(ns);       % Replication and reactivation of lysogens
       Q(1:ns,(ns+1):2*ns) = B*a*x(1)*diag(phivec);   % Influx from infections
        
       % Phages
       Q((ns+1):2*ns,1:ns) = alpha*M;                  % Influx from reactivation of lysogens
       Q((ns+1):2*ns,(ns+1):2*ns) = M * B*a*x(1)*diag(ones(1,ns) - phivec) - (deltaP + a*(x(1) + sum(x(3:(ns+2)))))*eye(ns); % Influx from infections, efflux from degradation and adsorption
       
       fx(3:end) = Q*x(3:end);
       
    end
end