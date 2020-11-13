function pL = phi(phimax,thres,conc)
% This function calculates the lysogeny propensity pL as function of the
% arbitrium concentration
% Lysogenization propensity phi is phimax if conc > threshold, 0 otherwise
% Returns a vector of length phimax, containing phi's for all strains

ns = length(phimax);
pL = zeros(1,ns);

for i = 1:ns
    if conc > thres(i)
        pL(i) = phimax(i);
    end
end

end