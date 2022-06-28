%mPSK modulation

clc;
clear all;
close all;

Rb = 1000; %bit rate
Tb = 1/Rb; %period of bits
Eb = 1; %Power of bit
m = 2; % PSK modulation 

T_flat_fading = 1; %Period of flat fading channel
%We declare the type of channel(s)
%'MRC' for Maximum Ratio Combining
%'MIMO' for Multiple Input-Multiple Output
Flat_fading = 'MIMO';

m_PSK = [2 4 8 16];
modulation = 'QAM'; %The type of the modulation 'PSK' or '16-QAM'

%execute for all PSK modulations (2,4,8,16)
% for(PSK=1:1:length(m_PSK))
%     m = m_PSK(PSK)
    
    
N = 2000; %number of bit in the input
num_of_slots = ceil(N/Rb);
bits = randi([0 1],1,N); %random 0s and 1s with equal probability 
Ts = Tb*log2(m); %period of symbols
Rs = 1/Ts; %Symbol Rate
OS = 1; %OVERSAMPLING rate

Ps = N*Eb; %total power of all bits-symbols
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


%We start a big for loop at which its iteration will create the
%serial-to-paraller handware. 1st will be the slot x1 (1-1000 bits) 
%and 2nd the slot x2(1001-2000 bits) etc
y = [];
for CH=1 : num_of_slots
    slot = bits((CH-1)*Rb+1:CH*Rb);
    
    %change input depending on the
    %number of bits by adding 0s to
    %the end of input string
    %example:we have 8PSK and
    %input=01-> input=010
    R = rem(length(slot),log2(m));
    input = slot;
    if R~=0
        for n=R: log2(m)-1
            input = [input 0];
        end
    end

    %%%%%%%%%% SYMBOL ENCODER %%%%%%%%%%

    %find the group of 0s and 1s depending on the modulation
    symbols_b = [];
    for n=0: length(input)/log2(m)-1
        symb = input([n*log2(m)+1:(n+1)*log2(m)]); 
        symbols_b = [symbols_b ;symb];   
    end


    %create symbols in complex form
    x = [];
    for n =1: size(symbols_b,1)
        j = find(ismember(lut_b,symbols_b(n,:),'rows'));
        x = [x lut(j)];
    end


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
    s = s1 + s2*1i;%After Tx convolution

    
%     % %plot of Real part after the Tx filter
%     t = [Ts/length(Tx):Ts/length(Tx):length(x)*Ts];
%     figure(1);
%     stem(t,s1,'b');
%     xlabel('Ts');
%     ylabel('S_I[n]');
%     title('Real part of symbol');
% 
%     %  %plot of Imaginary part after Tx filter
%     t = [Ts/length(Tx):Ts/length(Tx):length(x)*Ts];
%     figure(2);
%     stem(t,s2,'r');
%     xlabel('Ts');
%     ylabel('S_Q[n]');
%     title('Imaginary part of symbol');

    %%%%%%%%%% FLAT FADING %%%%%%%%%%
    
    
    if  strcmp(Flat_fading,'MRC')==1
        h1 = sqrt(1/sqrt(2))*(randn(1) + randn(1)*1i);
        h2 = sqrt(1/sqrt(2))*(randn(1) + randn(1)*1i);
        hn = [h1 ; h2];
        s1 = hn(1)*s;
        s2 = hn(2)*s; 
    elseif strcmp(Flat_fading,'MIMO')==1
        for  i = 1 : num_of_slots
            h1 = sqrt(1/sqrt(2))*(randn(1) + randn(1)*1i);
            h2 = sqrt(1/sqrt(2))*(randn(1) + randn(1)*1i);
            H(1:2,CH) = [h1 ; h2];
            s1 = H(:,CH)*s;
        end
    end
    
    
 
    
    
    
    %%%%%%%%%% AWGN %%%%%%%%%%
    %s_h = h*s;
    
       
    %%plot of constellation with flat fading and without AWGN 
    figure(11);
    plot(real(s1),imag(s1),'*r'); 
    hold on;
%     plot(real(s2),imag(s2),'*g'); 
%     hold on;
    plot(real(lut),imag(lut),'*b'); 
    hold on;
    xlabel('Real');
    ylabel('Imaginary');
    title('Constellation with flat fading and without AWGN');
    grid on;

    w = [];
    % %values of AWGN and after add it to the signal
    % %first channel
   
    if strcmp(Flat_fading,'MRC')==1
        for n=1: length(variance)
            w(n,:) = sqrt(variance(n))*(randn(length(s1),1) + randn(length(s1),1)*1i);
            y1(n,:) = s1 + w(n,:); 
            w(n,:) = sqrt(variance(n))*(randn(length(s2),1) + randn(length(s2),1)*1i);
            y2(n,:) = s2 + w(n,:); 
        end
        
        for n=1 :length(variance)
            figure(n+2);
            plot(real(y1(n,:)),imag(y1(n,:)),'*r'); 
            hold on;
            plot(real(y2(n,:)),imag(y2(n,:)),'*y'); 
            hold on;
            plot(real(lut),imag(lut),'*b'); 
            xlabel('Real');
            ylabel('Imaginary');
            title(['Constellation for SNRdb = ',num2str(SNR_db(n)),' and variance = ',num2str(variance(n)),' with Flat Fading']);
            grid on;
        end
    
    elseif strcmp(Flat_fading,'MIMO')==1
        for n=1: length(variance)
            w(n,:) = sqrt(variance(n))*(randn(length(s1),1) + randn(length(s1),1)*1i);
            y1(n,:) = s1(1,:) + w(n,:); 
        end
        
        for n=1 :length(variance)
            figure(n+2);
            plot(real(y1(n,:)),imag(y1(n,:)),'*r'); 
            hold on;
            plot(real(lut),imag(lut),'*b'); 
            xlabel('Real');
            ylabel('Imaginary');
            title(['Constellation for SNRdb = ',num2str(SNR_db(n)),' and variance = ',num2str(variance(n)),' with Flat Fading']);
            grid on;
        end
    end

  %%plot of constellation after add AWGN


    %%%%%%%%%% Equalizer/Channel inversion %%%%%%%%%%
   
    c = 0;%constant of B2 question

    %MRC
    if strcmp(Flat_fading,'MRC')==1
        for n=1:length(variance)
             mrc = (conj(hn).*(hn)/norm(hn')^2).*[y1(n,:); y2(n,:)];
             y0(n,:) = mrc(1,:)+ mrc(2,:);
             y(n,(CH-1)*Rs+1:CH*Rs) = y0(n,:);
        end
            figure(14);
            plot(real(y(7,:)),imag(y(7,:)),'*r'); 
            hold on;
            plot(real(lut),imag(lut),'*b'); 
            xlabel('Real');
            ylabel('Imaginary');
            title(['Constellation for SNRdb = ',num2str(SNR_db(7)),' and variance = ',num2str(variance(7)),' after the Equalizer']);
            grid on;
    end
    
    if strcmp(Flat_fading,'MIMO')==1
        for i = 1 :length(variance)
            y0(n,(CH-1)*Rs+1:CH*Rs) = y1(n,:);
        end
    end
    

    
end %end of slot loop
    
%%%MIMO%%%
if  strcmp(Flat_fading,'MIMO')==1
    for n =1: length(variance)
        y11 = y0(n,1:Rs); 
        y12 = y0(n,Rs+1:2*Rs);
        A = inv(conj(H)*(H'))*conj(H).*(H');
        y_mimo1 = A(1,1)*y11(1,:)+A(1,2)*y11(1,:);
        y_mimo2 = A(2,1)*y12(1,:)+A(2,2)*y12(1,:);
        y(n,:)=[y_mimo1 y_mimo2];
    end
    figure(14);
    plot(real(y(7,:)),imag(y(7,:)),'*r'); 
    hold on;
    plot(real(lut),imag(lut),'*b'); 
    xlabel('Real');
    ylabel('Imaginary');
    title(['Constellation for SNRdb = ',num2str(SNR_db(7)),' and variance = ',num2str(variance(7)),' after the Equalizer']);
    grid on;
end
%%%%%%%%%% Rx FILTER %%%%%%%%%%
    Rx = Tx;
    z1 = [];
    z2 = [];

    %  %convlution of y[n] = s[n]+w[n] with Rx
    for(n=1:1:length(variance))
       for(k=1:OS:length(s))
            b = conv(real(y(n,k:k+OS-1)),Rx);
            c = conv(imag(y(n,k:k+OS-1)),Rx);
            z1 = [z1,b]; % real part samples
            z2 = [z2,c]; % imaginary part samples
       end
    end
    z11 = [];
    z21 = [];
    z11 = reshape(z1(:),length(z1)/length(variance), []).';
    z21 = reshape(z2(:),length(z1)/length(variance), []).';

%     for(k=1:1:length(variance))
%         % %plot of Real part after the Rx filter
% 
%         figure(k+100);
%         stem(z11(k,:),'b');
%         xlabel('Ts');
%         ylabel('Z_I[n]');
%         title(['Real part of symbol after Rx filter for variance = ',num2str(variance(k))]);
% 
%         % %plot of Imaginary part after Rx filter
% 
%         figure(k+200);
%         stem(z21(k,:),'r');
%         xlabel('Ts');
%         ylabel('Z_Q[n]');
%         title(['Imaginary part of symbol after Rx filter for variance = ',num2str(variance(k))]);
%     end

    %%%%%%%%%% DECISION DEVISE %%%%%%%%%%
    z = [];
    z = z11 + z21*i; %z=the output of Rx filter for the different variances/SNRdb
    d1 = [];
    d2 = [];
    d11 = [];
    d21 = [];
    %sample the min/max values for each symbol
    for n=1: length(variance)
        
        d1 = reshape(z11(n,:),2*OS-1,[]).';
        d2 = reshape(z21(n,:),2*OS-1, []).';
        
        for k=1: size(d1,1)
            %for each symbol we take the min or max value
            
            %we check the if the values after the convolution
            %if they are positive or negative and we sample 
            %the min or max
            
            if d1(k,1)>0
                d11(n,k) = max(d1(k,:));
            else
                d11(n,k) = min(d1(k,:));
            end

            if d2(k,1)>0
                d21(n,k) = max(d2(k,:));
            else
                d21(n,k) = min(d2(k,:));
            end        
        end
    end
    
    d_sampled = [];
    d_sampled = d11 + d21*i;
    d = [];

    %for all the values of variance create d[n] the output of Decision devise
    for n=1: length(variance)

        output = [];
        %Go through all the sampled values in the first for loop
        %then loop in the look up tables find the nearest complex symbol
        %and give the binary value depending on the complex number of the
        %constellation
        %Maximum Likelihood Detector
        for k=1:length(d_sampled(n,:))

            dis_min = 100;

            for l=1: length(lut)-1
                distance = norm(lut(l) - d_sampled(n,k));

                if distance < dis_min
                    dis_min = distance;
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
        for k=1:length(input)
            if xor(input(k),d(n,k)) %input(k)~=d(n,k))
                c=c+1;
            end
        end
        error_counts(n,:) = c;
    end
    %BER_vector(PSK,:) = error_counts/N;
    
    BER_test = error_counts/N;


BER_theory = qfunc(sqrt(2*SNR));

figure(1);
sem1 = semilogy(SNR_db, BER_theory, 'bo-');
hold on;
sem2 = semilogy(SNR_db,BER_test,'rx-');
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
% 
% 
%     % %plot of Real part after the Tx filter
%     t = [Ts/length(Tx):Ts/length(Tx):length(x)*Ts];
%     figure(2);
%     stem(t,s1,'b');
%     xlabel('Ts');
%     ylabel('S_I[n]');
%     title('Real part of symbol after Tx filter');
% 
%     %  %plot of Imaginary part after Tx filter
%     t = [Ts/length(Tx):Ts/length(Tx):length(x)*Ts];
%     figure(3);
%     stem(t,s2,'r');
%     xlabel('Ts');
%     ylabel('S_Q[n]');
%     title('Imaginary part of symbol after Tx filter');
% 
% 
%     for(k=1:1:length(variance))
%         % %plot of Real part after the Rx filter
%         figure(k+3);
%         stem(z11(k,:),'b');
%         xlabel('Ts');
%         ylabel('Z_I[n]');
%         title(['Real part of symbol after Rx filter for variance = ',num2str(variance(k))]);
% 
%         % %plot of Imaginary part after Rx filter
%         figure(k+10);
%         stem(z21(k,:),'r');
%         xlabel('Ts');
%         ylabel('Z_Q[n]');
%         title(['Imaginary part of symbol after Rx filter for variance = ',num2str(variance(k))]);
%     end
%     
%         %%plot of constellation after add AWGN
%     for(n=1:1:length(variance))
%         figure(n+20);
%         plot(real(y(n,:)),imag(y(n,:)),'*r'); 
%         hold on;
%         plot(real(lut),imag(lut),'*b'); 
%         xlabel('Real');
%         ylabel('Imaginary');
%         title(['Constellation for SNRdb = ',num2str(SNR_db(n)),' and variance = ',num2str(variance(n))]);
%         grid on;
%     end
