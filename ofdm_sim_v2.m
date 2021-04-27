clc
clear
close

M=16;
N=600; 
N1 = N/10; % 10 is the number of packet // N1 is the size of each packet
cp_size=9;%round(N1/4);
a=1;
F=dftmtx(N1);

index = [-3*a+j*3*a;-a+3*j*a;a+3*a*j;3*a+3*j*a;-3*a+j*a;-a+j*a;a+j*a;3*a+a*j;-3*a-a*j;-a-j*a;a-j*a;3*a-j*a;-3*a-3*a*j;-a-3*a*j;a-a*3*j;3*a-3*a*j];
ofdm_symb=6;

for i=1:ofdm_symb
x(i,:)=randi([0 M-1],1,N1);
x_mod(i,:) =qammod(x(i,:),M);
end

% IFFT
for i=1:ofdm_symb
x_mod_ifft(i,:)=ifft(x_mod(i,:));
end

% adding cyclic prefix
for i=2:ofdm_symb
x_mod_ifft_cp(i,:)=[x_mod_ifft(i-1,end-cp_size + 1:end) x_mod_ifft(i,:)];
end
x_mod_ifft_cp(1,:)=[zeros(1,cp_size) x_mod_ifft(1,:)];
pilot=ones(1,length(x_mod_ifft_cp));
x_mod_ifft_cp1=cat(1,pilot,x_mod_ifft_cp);

% for i=1:ofdm_symb
% r(i,:)=C*x_mod_ifft_cp(i,:)';
% end
v = 1:10;
%fft
% for i=1:ofdm_symb
% r_F(i,:)=fft(r(i,:));
% end
% % observation ( compare results)
%    diag_matrix1=F*C(:,4:end)*inv(F); 
%    diag_matrix2=(fft(r_F*x_mod'));
npilot=1;
for i = 1:ofdm_symb+npilot
    r(i,:)= conv(x_mod_ifft_cp1(i,:),v);
end
rr = r(:,cp_size + 1:length(x_mod_ifft_cp1));
rr_fft = fft(rr);

for j=1:ofdm_symb+npilot
    for i=1:length(rr)
    dif=rr(j,i)-compare;
    end
end
% for i = 1:ofdm_symb+npilot
% compare(i,:) = ifft(fft(x_mod_ifft(i,:)).*fft(v,length(x_mod_ifft)));
% end
%equality = compare - rr;






