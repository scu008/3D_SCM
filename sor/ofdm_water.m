function lambda = ofdm_water(H, pcon, SNR, tol)
if nargin < 4, tol = 1e-5;
end
% Nsubcarrier = 128;
% tol = 1e-5;
% pcon = 10^(-10/10) * Nsubcarrier;
% h = Rayleigh_channel(7);
% H = fft(h,Nsubcarrier)*sqrt(Nsubcarrier);
% SNR = 0;

H = H * sqrt(length(H));
SNR_ = 10^(-SNR/10);
vec = SNR_./(abs(H).^2);

N = length(vec);

%first step of water filling
wline=min(vec)+pcon/N; %initial waterline
ptot=sum(max(wline-vec,0)); %total power for current waterline
 
%gradual water filling
while abs(pcon-ptot)>tol
    wline=wline+(pcon-ptot)/N;
    ptot=sum(max(wline-vec,0));
end
lambda = max(wline-vec,0);
% g = [vec'  max(wline-vec,0)'];
% bar(g,'stack');
% legend ('Noise Level','Power Level','')
% ylabel('Noise & Power Level','fontsize',12)
% xlabel('Number of Channels (N)','fontsize',12)
% title('Power Allocation for Waterfilling Alogorithm','fontsize',12)
% axis([0 N 0 mean(sum(g,2))])