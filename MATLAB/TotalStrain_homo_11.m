function out = TotalStrain_homo_11(P1, P2, P3, Constants)
    out = Constants.Us_11; % Thin film
%    out = 0; % Strain free
%    out = vol_avg(Constants.Q11 .* P1.^2 + Constants.Q12 .* (P2.^2 + P3.^2),Constants); % Stress free
end