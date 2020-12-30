function S_norm = plot_normalized_PSD(rho)

f = -0.5:0.001:0.5;
den=0;

for i = 2:length(rho)
    den = den + rho(i)*exp(-1j*2*pi*(i-1)*f);
end
S = 1./abs(1-den).^2;
S_norm = S/max(S);

figure
plot(f,10*log10(S_norm))
grid on
title('Normalized PSD')

end

