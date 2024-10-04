%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
% Created on   : July 2023
% Description  : QAM Slicer
%-----------------------------------------------------------------------------%

function o_data_v = my_slicer(i_data_v, M)
    
    %--------------------------%
    %           ERRORS
    %--------------------------%
    
    if mod(log2(M), 2) ~= 0 && M~=2
        error('The modulation levels must be M=2^k with even k (or k=1)')
    end
    
    %--------------------------%
    %   CONSTANTS & VARIABLES
    %--------------------------%
    
    sqrt_M = sqrt(M);
    
    %--------------------------%
    %          PROCESS
    %--------------------------%

    if M==2     % PAM-2
        
        % Output
        o_data_v = real(i_data_v);
        o_data_v(real(i_data_v)<0) = -1;
        o_data_v(real(i_data_v)>0) = 1;
        
    else        % QAM-M
        
        % Real part
        real_v = round((real(i_data_v) + (sqrt_M-1)) ./2 );
        real_v(real_v <= -1) = 0;
        real_v(real_v > (sqrt_M - 1)) = sqrt_M - 1;
        real_v = 2*real_v - sqrt_M + 1;
        
        % Imag part
        imag_v = round((imag(i_data_v) + (sqrt_M-1)) ./2 );
        imag_v(imag_v <= -1) = 0;
        imag_v(imag_v > (sqrt_M - 1)) = sqrt_M - 1;
        imag_v = 2*imag_v - sqrt_M + 1;
        
        % Output
        o_data_v = real_v + 1j * imag_v;

    end
    
end
