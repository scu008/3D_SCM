function [zeta_B, zeta_D, eta_B, eta_D] = appendA(N, i, j, B, D, Vrf)
zeta_B = 0;
zeta_D = 0;
eta_B = 0;
eta_D = 0;
for m = 1 : N
    for n = 1 : N
        if (m ~= i && n ~= i)
            zeta_B = zeta_B + conj(Vrf(m,j))*B(m,n)*Vrf(n,j);
            zeta_D = zeta_D + conj(Vrf(m,j))*D(m,n)*Vrf(n,j);
        end
    end
    if (m ~= i)
        eta_B = eta_B + B(i,m)*Vrf(m,j);
        eta_D = eta_D + D(i,m)*Vrf(m,j);
    end
    
end
zeta_B = B(i,i) + 2*real(zeta_B);
zeta_D = D(i,i) + 2*real(zeta_D);
end