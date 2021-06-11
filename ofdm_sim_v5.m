clc
clear
close

% parameters

QAM_order=16;                %16 QAM
M=QAM_order;
N_Subcarrier=48;             % Data subcarrier
OFDM_symbol_size=128;        % Number of OFDM symbol 
channel_size=10;
Channel = randi([0 M-1],1,channel_size);
Cyclic_Prefix=channel_size*2;
N_pilot=64;  %standard N_pilot=4
modulation_coefficient=1;


%initialization
N1=N_Subcarrier; 
ofdm_symb=OFDM_symbol_size;
cp_size=Cyclic_Prefix;
a=modulation_coefficient;
F=dftmtx(N1);
index = [-3*a+j*3*a;-a+3*j*a;a+3*a*j;3*a+3*j*a;-3*a+j*a;-a+j*a;a+j*a;3*a+a*j;-3*a-a*j;-a-j*a;a-j*a;3*a-j*a;-3*a-3*a*j;-a-3*a*j;a-a*3*j;3*a-3*a*j];
npilot=N_pilot;
 


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



[sx_x_mod_ifft_cp sy_x_mod_ifft_cp]=size(x_mod_ifft_cp);
for i = 1:ofdm_symb
    r(i,:)= conv(x_mod_ifft_cp(i,:),Channel);
end
rr = r(:,cp_size + 1:sy_x_mod_ifft_cp);
rr_fft = fft(rr);


[sx_x_mod_ifft sy_x_mod_ifft]=size(x_mod_ifft);
for i = 1:ofdm_symb
compare(i,:) = ifft((fft(x_mod_ifft(i,:))).*fft(Channel,(sy_x_mod_ifft))); % in time domain
end
equality = compare - rr;
equality(:,1:length(Channel-1))=[];

%%

pilot_symbol=randi([0 M-1],N_pilot,N1) + j*randi([0 M-1],N_pilot,N1);
x_mod_ifft_pilot=cat(1,pilot_symbol,x_mod);

% IFFT
for i=1:ofdm_symb+N_pilot
x_mod_ifft_pilot(i,:)=ifft(x_mod_ifft_pilot(i,:));
end

% adding cyclic prefix
for i=2:ofdm_symb+N_pilot
x_mod_ifft_cp_pilot(i,:)=[x_mod_ifft_pilot(i-1,end-cp_size + 1:end) x_mod_ifft_pilot(i,:)];
end
x_mod_ifft_cp_pilot(1,:)=[zeros(1,cp_size) x_mod_ifft_pilot(1,:)];



[sx_x_mod_ifft_cp_pilot sy_x_mod_ifft_cp_pilot]=size(x_mod_ifft_cp_pilot);

%adding channel 
for i = 1:ofdm_symb+N_pilot
    r(i,:)= conv(x_mod_ifft_cp_pilot(i,:),Channel);
end

% in receiver
rr = r(:,cp_size + 1:sy_x_mod_ifft_cp_pilot);
rr_fft = fft(rr);

%estimation
     for j=1:N_pilot
         %H(i,:,j)=(rr_fft(i,:).*conj(pilot_symbol(j)))./(abs(pilot_symbol(j))).^2;
         H(j,:)=(rr_fft(j,:).*conj(pilot_symbol(j)))/(abs(pilot_symbol(j))).^2;
     end
     
%averaging over all pilots
     for j=1:N_Subcarrier
     H_est(j)=sum(H(:,j))/N_pilot;
     end
 %end


[sx_x_mod_ifft sy_x_mod_ifft]=size(x_mod_ifft_pilot);
for i = 1:ofdm_symb
compare_est(i,:) = ifft(fft(x_mod_ifft_pilot(i,:)).*fft(H_est,(sy_x_mod_ifft)));
end
equality_est = compare_est - rr(N_pilot+1:end,:);

%% figures1 : compare 

figure
plot(abs(equality(1,:)))
title('proof equality , Diag(.)==fft(h_zp)  ')
%% figures2 : channel estimation
Channel_zp=[Channel zeros(1,length(H_est)-length(Channel))];
plot(abs(Channel_zp),'LineWidth',2)
hold on
plot(abs(H_est),'LineWidth',2)

grid on
grid minor
sgtitle('channel estimation')
legend('real channel','estimated channel')



