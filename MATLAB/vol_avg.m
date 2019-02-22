function out = vol_avg( f, Constants )

    out = trapz(trapz(trapz(f))) /  sum(sum(sum(Constants.in_film)));

end