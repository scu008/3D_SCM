function y = dopp_filter(fd, Nfft, spectrum)
% fd = Maximum Doppler frequency
% Nfft= Number of frequency domain points

% frequency spacing
df = 2*fd / (Nfft-1); 
f = -fd:df:fd;
y = zeros(1, length(f) );

% 기본 스펙트럼을 Jake model로 설정
if nargin < 3, spectrum = 'j'; end


if spectrum == 'j'
    % Jake model
    y = 1.5 ./ ( pi*fd*sqrt( 1 - (f/fd).^2 ) );
    y(1) = 2 * y(2) - y(3);
    y(end) = 2 * y(end-1) - y(end-2);
    
elseif spectrum == 's'
    % SUI model
    y = 1 - 1.72 * (f ./ fd).^2 + 0.785 * (f ./ fd).^4;
    
elseif spectrum == 'f'
    % Flat model
    y(:) = 1/(2*fd);
    
elseif spectrum == 'g'
    % Gaussian
    sigma = 0.5;
    y = (1/sqrt(2*pi*sigma^2) * exp( -1 * f.^2 / (2*sigma^2) ) );
    
end

% Output spectrum
y = ifftshift(y);