% mPSK-16QAM modulation
% MATLAB R2018a

clc;
clear all;
close all;

m = 2; % PSK modulation 


m_PSK = [2 4 8 16];
modulation = 'PSK'; %The type of the modulation 'PSK' or '16-QAM'

%execute for all PSK modulations (2,4,8,16)

% for z = 1:length(m_PSK)
% m = m_PSK(z);   


Tb = 1;   
Eb = 1;
N = log2(m)*64; %number of bit in the input
bits = randi([0 1],1,N); %random 0s and 1s with equal probability 
Ts = Tb*log2(m); %period of symbols
Rs = 1/Ts; %Symbol Rate
OS = 1; %OVERSAMPLING rate
SNR_db = [0:2:12]; %SNR in db



%variance and SNRlinear for each SNRdb
for n=1:length(SNR_db)
    SNR(n) = 10.^(SNR_db(n)/10);
    variance(n) = 1/SNR(n);
    N0(n) = 2*variance(n);
end

%contellation look up tables
%lut(n) contains the complex number
%and lut_b(n) contains lut(n)'s 
%binary encoding
if modulation == 'PSK'
    if m == 2
        lut = [1 -1];
        lut_b = [1 ; 0];
    elseif m == 4
        lut = [1+1i -1+1i -1-1i 1-1i];
        dec = [0:1:m-1];
        lut_b = de2bi(dec,'left-msb');
    else
        %find the phase of the modulation
        %then create polar coordinates
        %convert them to cartesian and 
        %finally create mPSK constelation 
        phase = 2*pi/m;
        theta = [0:phase:2*pi];
        r = sqrt(Eb);
        [x1,y1] = pol2cart(theta,r);
        lut = [x1 + y1*1i];
        lut = lut(1:end-1);
        dec = [0:1:m-1];
        lut_b = de2bi(dec,'left-msb');

    end
 elseif modulation == 'QAM'
    m = 16;
    Ts = Tb*log2(m); 
    Rs = 1/Ts; 
    dec = [0:1:15];
    lut(1) = [-3*sqrt(Eb) + 3*sqrt(Eb)*1i]; 
    lut(2) = [-1*sqrt(Eb) + 3*sqrt(Eb)*1i];
    lut(3) = [-1*sqrt(Eb) + 1*sqrt(Eb)*1i];
    lut(4) = [-3*sqrt(Eb) + 1*sqrt(Eb)*1i];
    lut(5) = [1*sqrt(Eb) + 3*sqrt(Eb)*1i];
    lut(6) = [3*sqrt(Eb) + 3*sqrt(Eb)*1i];
    lut(7) = [3*sqrt(Eb) + 1*sqrt(Eb)*1i];
    lut(8) = [1*sqrt(Eb) + 1*sqrt(Eb)*1i];
    lut(9) = [1*sqrt(Eb) - 1*sqrt(Eb)*1i];
    lut(10) = [3*sqrt(Eb) - 1*sqrt(Eb)*1i];
    lut(11) = [3*sqrt(Eb) - 3*sqrt(Eb)*1i];
    lut(12) = [1*sqrt(Eb) - 3*sqrt(Eb)*1i];
    lut(13) = [-3*sqrt(Eb) - 1*sqrt(Eb)*1i];
    lut(14) = [-1*sqrt(Eb) - 1*sqrt(Eb)*1i];
    lut(15) = [-1*sqrt(Eb) - 3*sqrt(Eb)*1i];
    lut(16) = [-3*sqrt(Eb) - 3*sqrt(Eb)*1i];
    lut_b = de2bi(dec,'left-msb');
end


    
% %     change input depending on the
% %     number of bits by adding 0s to
% %     the end of input string
% %     example:we have 8PSK and
% %     input=01-> input=010
% %     R = rem(length(slot),log2(m));
% %     input = slot;
% %     if R~=0
% %         for n=R: log2(m)-1
% %             input = [input 0];
% %         end
% %     end

    %%%%%%%%%% SYMBOL ENCODER %%%%%%%%%%

    %find the group of 0s and 1s depending on the modulation
    symbols_b = [];
    for n=0: length(bits)/log2(m)-1
        symb = bits([n*log2(m)+1:(n+1)*log2(m)]); 
        symbols_b = [symbols_b ;symb];   
    end


    %create symbols in complex form
    X = [];
    for n =1: size(symbols_b,1)
        j = find(ismember(lut_b,symbols_b(n,:),'rows'));
        X = [X lut(j)];
    end

    
    
    %%%%%%%%%% IDFT %%%%%%%%%%
    X1 = fftshift(real(X));
    X2 = fftshift(imag(X));
    x1 = ifft(X1);
    x2 = ifft(X2);
    x = x1 + x2*1i;
    
        
    
    %%%%%%%%%% Tx FILTER %%%%%%%%%%

    Tx = ones(1,OS);%filter creation

    s1 = [];
    s2 = [];
    %convolution of bits(encoded) with Ts filter
    for(n=1:1:length(x))
        b = conv(real(x(n)),Tx);
        c = conv(imag(x(n)),Tx);
        s1 = [s1 b]; % real part samples
        s2 = [s2 c]; % imaginary part samples
    end
    
    x_n = s1 + s2*1i;%After Tx convolution

    %%%%%%%%%% Add cyclic prefix %%%%%%%%%%
    cp = 16;
    k = 0;
    for n=length(x_n)+1:1:length(x_n)+cp
        k = k+1;
        x_n(n) = x_n(k);
    end
    

    %%%%%%%%%% LTI %%%%%%%%%%
    h = [0.9 + 0.9*1i 0.6 + 0.6*1i 0.3 + 0.3*1i];   
    y0 = cconv(x_n,h);
    v = zeros(1,61);
    h = [h v];
    
    
    %%%%%%%%%% AWGN %%%%%%%%%%
   

    w = [];
    % %values of AWGN and after add it to the signal
    for n=1: length(variance)
        w(n,:) = sqrt(variance(n))*(randn(length(y0),1) + randn(length(y0),1)*1i);
        y(n,:) = y0 + w(n,:); 
    end
    
    %%remove cyclic prefix %%
    for n=1:length(variance)
        y1(n,1:64) = y(n,1:64);
    end
    
    %%%%%%%%%% DFT %%%%%%%%%%
    for n=1:length(variance)
        Y1(n,:) = fftshift(real(y1(n,:)));
        Y2(n,:) = fftshift(imag(y1(n,:)));
        Y1(n,:) = fft(Y1(n,:));
        Y2(n,:) = fft(Y2(n,:));
        Y(n,:) = Y1(n,:)+Y2(n,:)*1i;
    end
    
    %%%%%%%%%% Equalizer/Channel inversion %%%%%%%%%%
    H = fft(h);
    
    for n=1:length(variance)
        for k=1:length(H)
            Y0(n,k) = Y(n,k)/H(k);
        end
    end
    
    %%%%%%%%% Decision Device %%%%%%%%%%

    d_sampled = Y0;
    d = [];

    %for all the values of variance create d[n] the output of Decision device
    for n=1: length(variance)

        output = [];
        %Go through all the sampled values in the first for loop
        %then loop in the look up tables find the nearest complex symbol
        %and give the binary value depending on the complex number of the
        %constellation
        %Maximum Likelihood Detector
        for k=1:length(d_sampled(n,:))

            dis_min = 1000;

            for l=1: length(lut)
                distance = norm(lut(l) - d_sampled(n,k));

                if distance< dis_min
                    dis_min = distance;
                    dis(n,k) = distance;
                    pos(n,k) = l;
                    j = l;
                end

            end
            output = [output lut_b(j,:)];

        end
        d(n,:) = output;
    end

    

%count of the wrong bits at the output of the system
    error_counts = [];
    for n=1:length(variance)
        c=0;
        for k=1:length(bits)
            if  bits(k)~=d(n,k) %%xor(bits(k),d(n,k))
                c=c+1;
            end
        end
        error_counts(n) = c;
        %BER_vector(z,n) = error_counts(n)/N;
        BER_test(n) = error_counts(n)/N;
    end
  
    
% end   %big for end
BER_theory = qfunc(sqrt(2*SNR));

figure(1);
sem1 = semilogy(SNR_db, BER_test, 'bo-');
% hold on;
% sem2 = semilogy(SNR_db,BER_vector(1,:),'rx-');
% hold on;
% sem3 = semilogy(SNR_db, BER_vector(2,:), 'g*-');
% hold on;
% sem4 = semilogy(SNR_db, BER_vector(3,:), 'k>-');
% hold on;
% sem5 = semilogy(SNR_db, BER_vector(4,:), 'm^-');
% grid on;
% legend([sem1(1),sem2(1),sem3(1),sem4(1),sem5(1)],'BPSK theory','BPSK - uncoded - AWGN','QPSK - uncoded - AWGN','8-PSK - uncoded - AWGN','16-PSK - uncoded - AWGN')
title('BER Results');
xlabel('Eb/No (dB)');
ylabel('BER');
axis([0 10 10^-6 1]);     
