function [S,A,L,P] = solve_ode(t,x0,phimaxvec,thresvec,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type)
% Function that solves the ODEs describing S, A, L and P dynamics
% Input:
% t         time vector
% x0        Input vector with initial values for S, A, vec(L), vec(P).
% phimaxvec Vector with phimax-values of the different strains
% thresvec  Vector with threshold-values of the different strains
% nphi      Number of different phi-values tested (= jump in strain vector to find the next thres value) 
% r         Growth rate bacteria
% K         Scaled carrying capacity of bacteria
% B         Effective burst size
% alpha     Reactivation rate lysogens
% deltaP    Decay rate phages
% a         Scaled adsorption rate phages to cells
% cL        Scaled production rate arbitrium by lysogens
% deltaA    Decay rate arbitrium
% u         Scaled uptake and breakdown rate arbitrium by cells
% mu        Mutation rate (prob of mutation in strain j)
% mut_type  Mutation type ("step" or "rand")

% Number of strains
ns = length(phimaxvec);

% Number of different threshold values
nthres = ns/nphi;

% Check threshold vector
if length(thresvec) ~= ns
    disp('Error solve_ode: thresvec should have length ns')
    return
end

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
    % Stepwise mutations: only small increase/decrease in phi_max and thres
    if mut_type == 'step'
        if nphi == 1
            Mphi = 1;
        else
            % Construct Mphi, matrix describing mutations in phi
            m_low = 0.5*mu*ones(1,nphi-1);
            m_up = 0.5*mu*ones(1,nphi-1);
            m_diag = [1-0.5*mu,(1-mu)*ones(1,nphi-2), 1-0.5*mu];
            Mphi = diag(m_diag) + diag(m_low,-1) + diag(m_up,1);
        end
        % Construct M, matrix describing all mutations
        if nthres == 1
            M = Mphi;
        else
            M = zeros(ns);
            M(1:nphi,1:nphi) = (1-0.5*mu) * Mphi;
            M((nphi+1):(2*nphi),1:nphi) = 0.5*mu * Mphi;
            for j=2:(nthres-1)
                M((nphi*(j-2)+1):(nphi*(j-1)),(nphi*(j-1)+1):(nphi*j)) = 0.5*mu * Mphi;
                M((nphi*(j-1)+1):(nphi*j),(nphi*(j-1)+1):(nphi*j)) = (1-mu) * Mphi;
                M((nphi*j+1):(nphi*(j+1)),(nphi*(j-1)+1):(nphi*j)) = 0.5*mu * Mphi;
            end
            M((nphi*(nthres-1)+1):(nphi*nthres),(nphi*(nthres-1)+1):(nphi*nthres)) = (1-0.5*mu) * Mphi;
            M((nphi*(nthres-2)+1):(nphi*(nthres-1)),(nphi*(nthres-1)+1):(nphi*nthres)) = 0.5*mu * Mphi;
        end
    % Random mutations
    else
        M = mu/(ns-1)*ones(ns);
        M = M - diag(diag(M)) + (1-mu)*eye(ns);
    end
end

% Solve the ODE
options1 = odeset('NonNegative',1:length(x0)); % Ensure solutions are positive
[t,x] = ode45(@RHS,t,x0,options1);

S = (x(:,1))';
A = (x(:,2))';
L = (x(:,3:(ns+2)))';
P = (x(:,(ns+3):end))';

    % Function that is solved by ODE-solver
    function fx = RHS(t,x)
       % Define output vector
       fx = zeros(2*ns+2,1);
       
       % Total number of cells
       N = x(1) + sum(x(3:(ns+2)));
       % Dynamics of susceptibles
       fx(1) = r*x(1)*(1-N/K) - B*a*x(1)*sum(x((ns+3):end));
       % Dynamics of arbitrium
       fx(2) = B*a*x(1)*sum(x((ns+3):end)) + cL*sum(x(3:(ns+2))) - x(2)*(deltaA + u*N);
       
       % Dynamics of lysogens and phages
       
       % Construction of full matrix Q that captures all effects of
       % lysogens and phages on the lysogen and phage dynamics
       Q = zeros(2*ns,2*ns);
       
       % phivec: contains current lysogenization propensities based on
       % current arbitrium concentration (stored in x(2))
       phivec = phi(phimaxvec,thresvec,x(2));
       
       % Lysogens
       Q(1:ns,1:ns) = (r*(1-(x(1) + sum(x(3:(ns+2))))/K)-alpha)*eye(ns);       % Replication and reactivation
       Q(1:ns,(ns+1):2*ns) = B*a*x(1)*diag(phivec);   % Influx from infections
        
       % Phages
       Q((ns+1):2*ns,1:ns) = alpha*M;                  % Influx from reactivation lysogens
       Q((ns+1):2*ns,(ns+1):2*ns) = M * B*a*x(1)*diag(ones(1,ns) - phivec) - (deltaP + a*(x(1) + sum(x(3:(ns+2)))))*eye(ns); % Influx from infections, efflux from degradation and adsorption
       
       fx(3:end) = Q*x(3:end);
       
    end
end