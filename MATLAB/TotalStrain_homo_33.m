function out = TotalStrain_homo_33(P1, P2, P3, Constants)
%     out = -( Constants.C12 * TotalStrain_homo_11(P1,P2,P3,Constants) + ...
%         Constants.C12 * TotalStrain_homo_22(P1,P2,P3,Constants) ) / Constants.C44; % Thin film

%    out = 0; % Strain free

%    out = vol_avg(Constants.Q11 * P3.^2 + Constants.Q12 .* (P1.^2 + P2.^2),Constants); % Stress free
    
    e_33 = vol_avg(Constants.Q11 * P3.^2 + Constants.Q12 .* (P1.^2 + P2.^2),Constants);
    e_11 = vol_avg(Constants.Q11 * P1.^2 + Constants.Q12 .* (P3.^2 + P2.^2),Constants);
    e_22 = vol_avg(Constants.Q11 * P2.^2 + Constants.Q12 .* (P3.^2 + P1.^2),Constants);
    out = (Constants.C12 * e_11 + Constants.C12 * e_22 + Constants.C44 * e_33 -(Constants.C12 * TotalStrain_homo_11(P1,P2,P3,Constants) + ...
         Constants.C12 * TotalStrain_homo_22(P1,P2,P3,Constants)) ) / Constants.C44;
    
end