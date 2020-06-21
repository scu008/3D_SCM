function [Output] = interleaver(Input, mode)

% mode : 1 ==> Interleaver   ,  mode : 2  ==>  Deinterleaver


N_tint = 10;
if mode == 1
    for k = 1:(length(Input)/(N_tint*10))
        for m = 1:N_tint*10
            Output((k-1)*(N_tint*10)+m) = Input( floor((m-1)/N_tint) + 10 * rem(m-1,N_tint) +1 +(k-1)*(N_tint*10) );
        end
        for m = floor( (length(Input)/(N_tint*10)) )*N_tint*10 + 1 :length(Input)
            Output(m) = Input(m);
        end
    end
else


    % Tone De-interleaver
    for k = 1:length(Input)/(N_tint*10)
        for m = 1:N_tint*10
            Output( floor((m-1)/N_tint) + 10 * rem(m-1,N_tint) +1 +(k-1)*(N_tint*10) ) = Input((k-1)*(N_tint*10)+m);
        end
        for m = floor( (length(Input)/(N_tint*10)) )*N_tint*10 + 1 :length(Input)
            Output(m) = Input(m);
        end
        
    end

end
