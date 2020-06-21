function [ s_hat ] = ML_detect( Y, H, mod_type )
    % ML detection을 수행하는 함수
    
    % 채널 행렬의 사이즈를 구한다.
    [fft_len, Mr, Mt] = size(H);
    all = 2^(Mt * mod_type);
    
    % 심볼 pool 생성
    sym_pool = base_mod( de2bi( 0 : all - 1, mod_type * Mt), mod_type).';
    
    % 반복 행렬 선언
    s_hat = zeros(Mt,fft_len);
    
    
    % 모든 서브케리어에 대해 수행
    for k = 1 : fft_len
        
        % 서브케리어 하나의 채널 행렬을 저장
        t_H(:,:) = H(k,:,:);
        
        % 거리 계산
        t_dis = repmat( Y(:,k), 1, all ) - (t_H * sym_pool);
        t_dis = sum( real(t_dis).^2 + imag(t_dis).^2 );
        [~, idx] = min(t_dis);
        
        % detect한 심볼 저장
        s_hat(:,k) = sym_pool(:,idx);
        
    end
end
