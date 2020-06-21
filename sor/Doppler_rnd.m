function [coeff] = Doppler_rnd(len, tap_len, fs, dopp_f, mode)
% 시변 채널을 위한 랜덤 시퀀스 출력

% len - 송신 심볼 길이
% tap_len - 다중 경로 수
% fs - sampling frequency
% dopp_f - 1 by n or 2 by n 도플러 주파수 행렬 (첫 행의 엘리먼트는 탭의 도플러 주파수, 둘 째 행의 엘리먼트는 탭의 도플러 스펙트럼 )
% dopp_f의 경우 도플러 주파수만 존재하는 벡터로 사용 가능 ( default spectrum - Jake model)
% dopp_f에 스칼라 값을 넣을 시 모든 탭에 동일한 도플러 효과 적용
% mode - 1: NLOS, mode - 2: LOS
% 작동을 위해선 추가로 dopp_filter.m 필요


% 매개변수에 따른 초기화
if nargin < 5, mode = 1; end

% Rician 설정
K = zeros(tap_len,1);
if mode == 2, K(1) = 10^(30/10); end
k1 = sqrt(K./(K+1));
k2 = sqrt(1./(K+1));

% 도플러 적용 여부 결정
if sum( dopp_f(1,:) ) == 0
    % 가우시안 랜덤 벡터 생성
    rnd = sqrt(0.5) * (randn(tap_len,1) + 1j*randn(tap_len,1));
    rnd = k1 + k2 .* rnd;
    rnd = rnd * ones(1,len);
else
    % 벡터 생성
    rnd = zeros(tap_len, len);
    
    for i = 1:tap_len
        
        % index 조정
        j = i;
        [f_row, f_len] = size(dopp_f);
        if j > f_len, j = f_len; end
        
        % 도플러 주파수와 스펙트럼을 설정
        df = dopp_f(1,j) + 0;
        if f_row == 2
            spectrum = dopp_f(2,j);
        else
            spectrum = 'j';
        end
        
        % 벡터 길이 정의
        Nfft = 2^max( 3, nextpow2(2*df/fs*len) );
        Nifft = ceil( Nfft*fs/( 2*df ) );
        
        % 가우시안 벡터 생성
        rnd_re = randn(1,Nfft) * sqrt(0.5);
        rnd_im = randn(1,Nfft) * sqrt(0.5);
        
        % 도플러 필터링 수행
        doppler_coeff = dopp_filter(df,Nfft, spectrum);
        f_rnd_re = fft(rnd_re) .* doppler_coeff;
        f_rnd_im = fft(rnd_im) .* doppler_coeff;
        rnd_re = ifft( [f_rnd_re(1:Nfft/2) zeros(1,Nifft-Nfft) f_rnd_re(Nfft/2+1:Nfft)] );
        rnd_re = rnd_re * sqrt( 0.5/mean( abs(rnd_re).^2 ));
        rnd_im = ifft( [f_rnd_im(1:Nfft/2) zeros(1,Nifft-Nfft) f_rnd_im(Nfft/2+1:Nfft)] );
        rnd_im = rnd_im * sqrt( 0.5/mean( abs(rnd_im).^2 ));
        
        
        % 채널 계수 발생
        if k1(i) ~= 0 && mode == 2
            % Rician
            t = 0 : 1/fs : (len-1)/fs;
            rnd(i,:) = k1(i)*exp( 1j .* ( 2*pi .* df .* cos(rand()*2*pi) .* t ) )  +  k2(i)*( rnd_re(1:len) - 1j*rnd_im(1:len) );
        else
            % Rayleigh
            rnd(i,:) = rnd_re(1:len) - 1j*rnd_im(1:len);
        end
        
    end
end

% 랜덤 시퀀스 출력
coeff = rnd;

end