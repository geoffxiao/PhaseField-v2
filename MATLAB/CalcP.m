% Calculate P within the film
P1_film = P1(:,:,interface_index+1:film_index-1);
P2_film = P2(:,:,interface_index+1:film_index-1);
P3_film = P3(:,:,interface_index+1:film_index-1);

P1_val = vol_avg(abs(P1_film),Constants)
P2_val = vol_avg(abs(P2_film),Constants)
P3_val = vol_avg(abs(P3_film),Constants)

P_val = sqrt(P1_val^2+P2_val^2+P3_val^2)