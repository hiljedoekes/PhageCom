function z = one_transfer_results(t,x0,phimaxvec,thresvec,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type)
% Function that returns the results of simulating one transfer episode

[S,A,L,P] = solve_ode(t,x0,phimaxvec,thresvec,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
z = [S(end), A(end), L(:,end)', P(:,end)'];

end