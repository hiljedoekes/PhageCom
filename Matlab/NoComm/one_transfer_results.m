function z = one_transfer_results(t,x0,phivec,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type)
% Return the results of simulating a single transfer episode

[S,A,L,P] = solve_ode(t,x0,phivec,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
z = [S(end), A(end), L(:,end)', P(:,end)'];

end