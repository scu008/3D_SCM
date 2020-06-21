function [ s_hat, W ] = MMSE_detect( Y, H, No )
    % MMSE detection을 수행하는 함수

    % 채널 행렬의 사이즈를 구한다.
    [fft_len, Mr, Mt] = size(H);
    
    % 반복 행렬 선언
    s_hat = zeros(Mt,fft_len);
    t_Y = [Y; zeros(Mt,fft_len)];
    W = zeros(fft_len,Mt,Mr);
    
    % 모든 서브케리어에 대해 수행
    for k = 1 : fft_len
        
        % 서브케리어 하나의 채널 행렬을 저장
        t_H(:,:) = H(k,:,:);
        e_H = [t_H; sqrt(No) * eye(Mt)];
        
        % MMSE 역행렬 계산
        G = inv( e_H' * e_H ) * e_H';
        
        % detect한 심볼 저장
        s_hat(:,k) = G * t_Y(:,k);
        W(k,:,:) = G;
        
    end

end

