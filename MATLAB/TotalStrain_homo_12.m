function out = TotalStrain_homo_12(P1, P2, P3, Constants)
    out = Constants.Us_12; % thin film
%    out = 0; % strain free
%    out = Constants.Q44 .* vol_avg(P1 .* P2,Constants); % Stress Free
end