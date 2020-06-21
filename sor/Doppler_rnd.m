function [coeff] = Doppler_rnd(len, tap_len, fs, dopp_f, mode)
% �ú� ä���� ���� ���� ������ ���

% len - �۽� �ɺ� ����
% tap_len - ���� ��� ��
% fs - sampling frequency
% dopp_f - 1 by n or 2 by n ���÷� ���ļ� ��� (ù ���� ������Ʈ�� ���� ���÷� ���ļ�, �� ° ���� ������Ʈ�� ���� ���÷� ����Ʈ�� )
% dopp_f�� ��� ���÷� ���ļ��� �����ϴ� ���ͷ� ��� ���� ( default spectrum - Jake model)
% dopp_f�� ��Į�� ���� ���� �� ��� �ǿ� ������ ���÷� ȿ�� ����
% mode - 1: NLOS, mode - 2: LOS
% �۵��� ���ؼ� �߰��� dopp_filter.m �ʿ�


% �Ű������� ���� �ʱ�ȭ
if nargin < 5, mode = 1; end

% Rician ����
K = zeros(tap_len,1);
if mode == 2, K(1) = 10^(30/10); end
k1 = sqrt(K./(K+1));
k2 = sqrt(1./(K+1));

% ���÷� ���� ���� ����
if sum( dopp_f(1,:) ) == 0
    % ����þ� ���� ���� ����
    rnd = sqrt(0.5) * (randn(tap_len,1) + 1j*randn(tap_len,1));
    rnd = k1 + k2 .* rnd;
    rnd = rnd * ones(1,len);
else
    % ���� ����
    rnd = zeros(tap_len, len);
    
    for i = 1:tap_len
        
        % index ����
        j = i;
        [f_row, f_len] = size(dopp_f);
        if j > f_len, j = f_len; end
        
        % ���÷� ���ļ��� ����Ʈ���� ����
        df = dopp_f(1,j) + 0;
        if f_row == 2
            spectrum = dopp_f(2,j);
        else
            spectrum = 'j';
        end
        
        % ���� ���� ����
        Nfft = 2^max( 3, nextpow2(2*df/fs*len) );
        Nifft = ceil( Nfft*fs/( 2*df ) );
        
        % ����þ� ���� ����
        rnd_re = randn(1,Nfft) * sqrt(0.5);
        rnd_im = randn(1,Nfft) * sqrt(0.5);
        
        % ���÷� ���͸� ����
        doppler_coeff = dopp_filter(df,Nfft, spectrum);
        f_rnd_re = fft(rnd_re) .* doppler_coeff;
        f_rnd_im = fft(rnd_im) .* doppler_coeff;
        rnd_re = ifft( [f_rnd_re(1:Nfft/2) zeros(1,Nifft-Nfft) f_rnd_re(Nfft/2+1:Nfft)] );
        rnd_re = rnd_re * sqrt( 0.5/mean( abs(rnd_re).^2 ));
        rnd_im = ifft( [f_rnd_im(1:Nfft/2) zeros(1,Nifft-Nfft) f_rnd_im(Nfft/2+1:Nfft)] );
        rnd_im = rnd_im * sqrt( 0.5/mean( abs(rnd_im).^2 ));
        
        
        % ä�� ��� �߻�
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

% ���� ������ ���
coeff = rnd;

end