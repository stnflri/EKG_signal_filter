clear

load('C:\Users\Iustin\Desktop\psb proiect\iaf1_afwm.mat');
X=val(6,:);
Fs=1000; %fr esantionare
Ts=1/Fs; %pas
L=10000; %esantioane
t=(0:(L-1))*Ts; %Vectorul timp
figure()
subplot(311);
plot(t,X)
xlabel('timp[s]'); ylabel('amplitudine[mV]');title('Semnal in timp')
grid on;
y=X+sin(2*pi*50*t)*100;
Y=fft(y);
amp=abs(Y);
la=length(amp);
amp1=amp(1:(la/2));
%amp2=zeros(1,la/2);
%amp2(1,500)=amp(1,500)*25;
%amp3=amp2+amp1;
f1=0:0.1:(1000-0.1);
lb=length(f1);
f=f1(1:(lb/2));
subplot(312)
plot(t,y)
title('Semnal cu zgomot in timp')
grid on
subplot(313)
plot(f,amp1)
xlabel('frecventa[Hz]'); ylabel('amplitudine[mV]');title('Semnal in frecventa cu 50 Hz zgomot')
grid on;


%Filtru ideal

zer = zeros(1, length(amp));
zer (1:490) = 1;
zer (510:9490) = 1;
zer(9510:10000)=1;
semnal_filtru_ideal = amp.*zer;
semnal_final_ideal = ifft(Y.*zer);
figure()
subplot(211)
plot(t(1:length(semnal_final_ideal)), semnal_final_ideal)
title('Semnal filtrat cu filtru ideal'), xlabel('timp[s]'), ylabel('amplitudine')
grid on
subplot(212)
plot(f1, semnal_filtru_ideal)
title('semnal in frecventa dupa aplicarea filtrului ideal'), xlabel('frecventa[Hz]'), ylabel('amplitudine')
grid on

%Filtru fir1

firi = fir1 (40, 0.05, 'low');
figure()
freqz(firi)
semnal_firi = filter (firi, 0.27, y)/2;
firi_filtrat = abs ( fft ( semnal_firi ) );
amp_firi = firi_filtrat (1:(length(firi_filtrat)/2));
amp_firi_final = zeros (1, 10000);
amp_firi_final (1:5000) = amp_firi;
j=5000;
for i = length(amp_firi):-1:1
    amp_firi_final(j) = amp_firi(i);
    j=j+1;
end
figure()
subplot(211)
plot(f, amp_firi_final(1:(length(amp_firi_final)/2)))
title('Semnal fir1 in frecventa'); xlabel('frecventa[Hz]'); ylabel('amplitudine')
grid on
subplot(212)
plot(t, semnal_firi)
title('Semnal fir1 in timp'); xlabel('timp[s]'); ylabel('amplitudine')
grid on

%Filtru fir2

firii = fir2(45,[0 0.08 0.08 1],[1 1 0 0]);

figure()
freqz(firii)

semnal_firii = filter(firii, 0.17, y)/4;
firii_filtrat = abs(fft(semnal_firii));
amp_firii = firii_filtrat (1:(length(firii_filtrat)/2));
amp_firii_final = zeros (1, 10000);
amp_firii_final (1:5000) = amp_firii;
j=5000;
for i = length(amp_firii):-1:1
    amp_firii_final(j) = amp_firii(i);
    j=j+1;
end
figure()
subplot(211)
plot(f, amp_firii_final(1:(length(amp_firii_final)/2)))
title('Semnal fir2 in frecventa'); xlabel('frecventa[Hz]'); ylabel('amplitudine');
grid on
amp_firi_t = ifft (amp_firii_final);
subplot(212)
plot(t, semnal_firii)
title('Semnal fir 2 in timp'); xlabel('ftimp[s]'); ylabel('amplitudine');
grid on

%Filtru eliptic

ft=40;
fop=50;
Wp=ft/(Fs/2);
Ws=fop/(Fs/2);
Rp = 3; %riplu in banda de trecere (db)
Rs = 40; %riplu in banda de oprire (db)
figure;
subplot(211)
[N,Wp]=ellipord(Wp,Ws,Rp,Rs);

[B,F]=ellip(N,Rp,Rs,Wp);

subplot(212)
freqz(B,F);

elip_filtrat = filter(B,F,y)*10;
amp_elip = abs(fft(elip_filtrat));
figure;
subplot(211)
plot(f, amp_elip(1:(length(amp_elip)/2)))
title('Semnal filtru eliptic in frecventa'); 
grid on
subplot(212)
plot(t, elip_filtrat);
title('Semnal filtru eliptic in timp');
grid on

%Filtru Cebasev 1

figure;
subplot(211)
[N,Wn]=cheb1ord(Wp,Ws,Rp,Rs,'s');
[B,D]=cheby1(N,Rp,Wp);
subplot(212)
freqz(B,D);

csvi_filtrat = filter(B,D,y);
amp_cvsi = abs(fft(csvi_filtrat));
figure;
subplot(211)
plot(f, amp_cvsi(1:5000))
title('Semnal filtru cebasev1 in frecventa'); 
grid on
subplot(212)
plot(t, csvi_filtrat)
title('Semnal filtru cebasev1 in timp');
grid on

%Filtru Cebasev 2

figure;
subplot(211)
[N,Ws]=cheb2ord(Wp,Ws,Rp,Rs,'s');
[B,E]=cheby2(N,Rs,Ws);
subplot(212)
freqz(B,E);

cvsii_filtrat = filter(B,E,y);
amp_cvsii = abs(fft(cvsii_filtrat));
figure;
subplot(211)
plot(f, amp_cvsii(1:5000))
title('Semnal Cebasev 2 in frecventa')
grid on
subplot(212)
plot(t, cvsii_filtrat)
title('Semnal Cebasev 2 in timp')
grid on

%Filtru butterworth

[B,C] = butter(16,45/500);
btt = filter(B,C,y);
amp_btt = abs(fft(btt));

figure;
freqz(btt);

figure;
subplot(211)
plot(f, amp_btt(1:5000))
title('Semnal butterworth in frecventa')
grid on
subplot(212)
plot(t, btt)
title('Semnal butterworth in timp')
grid on



diferenta_ideal = mean(abs(semnal_final_ideal - X))
diferenta_fir1 = mean(abs(semnal_firi - X))
diferenta_fir2 = mean(abs(semnal_firii - X))
diferenta_eliptic = mean(abs(elip_filtrat - X))
diferenta_ceb1 = mean(abs(csvi_filtrat - X))
diferenta_ceb2 = mean(abs(cvsii_filtrat - X))
diferenta_butter = mean(abs(btt - X))