function out = TotalStrain_homo_23(P1, P2, P3, Constants)
%    out = 0; % Strain free
%    out = Constants.Q44 .* vol_avg(P2 .* P3,Constants); % Stress free
    out = Constants.Q44 .* vol_avg(P2 .* P3,Constants);
end