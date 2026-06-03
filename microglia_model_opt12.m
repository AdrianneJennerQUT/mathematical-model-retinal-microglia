function dydt = microglia_model_opt12(t, y, p)
% State unpacking
O = y(1);
I = y(2);
G = y(3);
R = y(4);
B = y(5);

hill_G = B / (B + p.eta_G);
hill_O = B / (B + p.eta_O);
flux_G = p.k_G * G * (1 - G / p.G0);

dO = -p.k_OI * O * hill_G + p.k_IO * I * hill_O;

dI =  p.k_OI * O * hill_G ...
    - p.k_IO * I * hill_O ...
    - p.k_IG * I * hill_G ...
    - flux_G;

dG =  p.k_IG * I * hill_G ...
    + flux_G;

dR = -p.k_R * R;

dB =  p.F * p.k_R * (p.Vol_P / p.Vol_B) * R ...
    - p.k_B * B;

dydt = [dO; dI; dG; dR; dB];
end
