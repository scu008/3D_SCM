function y = FD_fading(sym, coeff)
% 3���� MIMO ä���� �����ϴ� �Լ�

% sym - �۽� �ɺ� ���� ( (N_h*N_v) by 1 ) = ���� ���� ���׳� ���� �� by 1
% coeff - FD_channel�� ���
% y - ���� ���׳��� ���� �ɺ� ���� �Ǵ� ���

% �Ű������� ���� �ʱ�ȭ
%=========================================================
[tap_len, sym_len, Nrx, Ntx] = size(coeff);

% ���� ���� ��ȣ ���
%==================================================================
% �۽� ��ȣ�� ä�� ��� ����
temp_y = zeros(tap_len, sym_len, Nrx);
for i= 1:Nrx
    for j= 1:Ntx        
        % ���� ���׳����� ä�ΰ� �۽� ��ȣ ����
        temp_y(:,:,i) = temp_y(:,:,i) + coeff(:,:,i,j) .* ( ones(tap_len,1) * sym(j,:) ) ;
    end
end

% ���� ��� ����
y = zeros(Nrx,sym_len + tap_len -1);
for m = 1:Nrx
    % ���� ��� �迭 ����
    flat_m = zeros(tap_len, sym_len + tap_len -1);
    
    for i = 1:tap_len
        flat_m(i, i:i + sym_len - 1) = temp_y(i,:,m);
    end
    
    % ���� ��� ��ø
    y(m,:) = sum( flat_m, 1 );
end


end