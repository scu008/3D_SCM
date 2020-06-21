function y = FD_fading(sym, coeff)
% 3차원 MIMO 채널을 적용하는 함수

% sym - 송신 심볼 벡터 ( (N_h*N_v) by 1 ) = 수평 수직 안테나 수의 곱 by 1
% coeff - FD_channel의 계수
% y - 수신 안테나에 따른 심볼 벡터 또는 행렬

% 매개변수에 따른 초기화
%=========================================================
[tap_len, sym_len, Nrx, Ntx] = size(coeff);

% 최종 수신 신호 계산
%==================================================================
% 송신 신호에 채널 계수 적용
temp_y = zeros(tap_len, sym_len, Nrx);
for i= 1:Nrx
    for j= 1:Ntx        
        % 수신 안테나마다 채널과 송신 신호 내적
        temp_y(:,:,i) = temp_y(:,:,i) + coeff(:,:,i,j) .* ( ones(tap_len,1) * sym(j,:) ) ;
    end
end

% 다중 경로 적용
y = zeros(Nrx,sym_len + tap_len -1);
for m = 1:Nrx
    % 다중 경로 배열 생성
    flat_m = zeros(tap_len, sym_len + tap_len -1);
    
    for i = 1:tap_len
        flat_m(i, i:i + sym_len - 1) = temp_y(i,:,m);
    end
    
    % 다중 경로 중첩
    y(m,:) = sum( flat_m, 1 );
end


end