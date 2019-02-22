function out = TotalStrain_homo_22(P1, P2, P3, Constants)
    out = Constants.Us_22; % thin film
%     out = 0; % strain free
%     out = vol_avg(Constants.Q11 * P2.^2 + Constants.Q12 .* (P1.^2 + P3.^2),Constants); % Stress free 
end