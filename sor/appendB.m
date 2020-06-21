function Vrf = appendB(theta_1, theta_2, zeta_B, zeta_D, eta_B, eta_D, A, N)
V_1 = exp(-1j*theta_1);
V_2 = exp(-1j*theta_2);
inv_A = inv(A);
f_1 = N*trace(inv_A) - N*(zeta_B + 2*real(conj(V_1) * eta_B))/(1 + zeta_D + 2*real(conj(V_1) * eta_D));
f_2 = N*trace(inv_A) - N*(zeta_B + 2*real(conj(V_2) * eta_B))/(1 + zeta_D + 2*real(conj(V_2) * eta_D));
if f_1 >= f_2
    Vrf = V_2;
else
    Vrf = V_1;
end
end