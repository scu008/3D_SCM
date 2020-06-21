h = Rayleigh_channel(10);

A = fft(h, 30)*sqrt(30);
A = log(1+ abs(A).^2);
A = A'
plot(A);