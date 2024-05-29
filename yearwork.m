clear all
close all
load("ponte_mkr");
ndgl=82;
ntot=90;
MFF=M(1:ndgl, 1:ndgl);
KFF=K(1:ndgl, 1:ndgl);
CFF=R(1:ndgl, 1:ndgl);
MFC=M(1:ndgl, ndgl+1:ntot);
KFC=K(1:ndgl, ndgl+1:ntot);
CFC=R(1:ndgl, ndgl+1:ntot);
MCF=M(ndgl+1:ntot, 1:ndgl);
KCF=K(ndgl+1:ntot, 1:ndgl);
CCF=R(ndgl+1:ntot, 1:ndgl);
MCC=M(ndgl+1:ntot, ndgl+1:ntot);
KCC=K(ndgl+1:ntot, ndgl+1:ntot);
CCC=R(ndgl+1:ntot, ndgl+1:ntot);
%------------------------------------
% point 1/2 - natural frequencies and modes of vibration
[eigenvectors eigenvalues]=eig(MFF\KFF);
freq=sqrt(diag(eigenvalues))/2/pi;
[freq_sort,idx]=sort(freq);
freq_sort;
eigenvectors=eigenvectors(1:end,idx);
%------------------------------------
% point 3 - damped natural frequencies
alpha=0.8;
beta=3.0e-5;
for k=1:6;
wi=2*pi*freq_sort(k);
wdi=sqrt(1-2*[alpha/(2*wi)+beta*wi*0.5]^2)*wi;
damped_freq(k)=wdi/2/pi;
h(k)=alpha/2/wdi+0.5*beta*wdi;
end
1
%------------------------------------
% point - 4 frequency response
i=sqrt(-1);
vett_f=0:0.01:24;
b=idb(9,2);
a=idb(7,2);
vect_F0=zeros(ndgl,1);
vect_F0(b)=1;
for k=1:length(vett_f)
ome=vett_f(k)*2*pi;
A=-ome^2*MFF+i*ome*CFF+KFF;
x=A\vect_F0;
y7=x(a);
y9=x(b);
acc_nodo7=-ome^2*y7;
acc_nodo9=-ome^2*y9;
mod1(k)=abs(acc_nodo7);
fas1(k)=angle(acc_nodo7);
mod2(k)=abs(acc_nodo9);
fas2(k)=angle(acc_nodo9);
end
figure(1)
subplot 211;plot(vett_f,mod1);grid;
hold
xlabel('Freq. [Hz]');
subplot 212;plot(vett_f,fas1);grid
hold
xlabel('Freq. [Hz]');
figure(2)
subplot 211;plot(vett_f,mod2);grid;
hold
xlabel('Freq. [Hz]');
subplot 212;plot(vett_f,fas2);grid
hold
xlabel('Freq. [Hz]');
phi_3=eigenvectors(:,1:3);
M_3=phi_3'*MFF*phi_3;
K_3=phi_3'*KFF*phi_3;
R_3=phi_3'*CFF*phi_3;
[eigenvectors eigenvalues]=eig(M_3\K_3);
freq=sqrt(diag(eigenvalues))/2/pi;
[freq_sort,idx]=sort(freq);
freq_sort;
eigenvectors=eigenvectors(1:end,idx);
2
%------------------------------------
% point 5 - frequency response related to first 3 modes
i=sqrt(-1);
vett_f=0:0.01:24;
vect_F1_q=phi_3'*vect_F0;
for k=1:length(vett_f)
ome=vett_f(k)*2*pi;
A_3=-ome^2*M_3+i*ome*R_3+K_3;
q=A_3\vect_F1_q;
x=phi_3*q;
y7=x(a);
y9=x(b);
acc_nodo7=-ome^2*y7;
acc_nodo9=-ome^2*y9;
mod3(k)=abs(acc_nodo7);
fas3(k)=angle(acc_nodo7);
mod4(k)=abs(acc_nodo9);
fas4(k)=angle(acc_nodo9);
end
figure(1)
subplot 211;plot(vett_f,mod3);grid on;
hold
title('Acceleration node A');
xlabel('Freq. [Hz]');
subplot 212;plot(vett_f,fas3);grid on;
xlabel('Freq. [Hz]');
figure(2)
subplot 211;plot(vett_f,mod4);grid on;
hold
title('Acceleration node B');
xlabel('Freq. [Hz]');
subplot 212;plot(vett_f,fas4);grid on;
xlabel('Freq. [Hz]');
3
%------------------------------------
% point 6 - frequency response of bending moment at point C
i=sqrt(-1);
vett_f=0:0.01:24;
xa=idb(7,1);
ya=idb(7,2);
ta=idb(7,3);
xb=idb(9,1);
yb=idb(9,2);
tb=idb(9,3);
xc=idb(8,1);
vect_F0=zeros(ndgl,1);
vect_F0(b)=1;
for k=1:length(vett_f)
ome=vett_f(k)*2*pi;
A=-ome^2*MFF+i*ome*CFF+KFF;
x=A\vect_F0;
x7=x(xa);
y7=x(ya);
t7=x(ta);
x9=x(xb);
y9=x(yb);
t9=x(tb);
eps=x(xc);
Lk=5;
fvpp=[0, 12*eps/Lk-6/Lk, 6*eps-4, 0, -12/Lk*eps+6/
Lk, 6*eps-2];
xk=[x7 y7 t7 x9 y9 t9];
wpp=fvpp*xk';
E=2.06e11;
I=2.313e-4;
EI=E*I;
M_c=EI*wpp;
mod5(k)=abs(M_c);
fas5(k)=angle(M_c);
end
figure(5)
subplot 211;plot(vett_f,mod5);grid
title('M_c/F');
xlabel('Freq. [Hz]');
subplot 212;plot(vett_f,fas5);grid
xlabel('Freq. [Hz]');
4
%------------------------------------
% point 7 - bridge response due to a seismic motion of the ground
% Load seismic displacement data from a text file
data = load('seismic_displ.txt');
% Extract time and displacement values
time = data(:, 1);
y1 = data(:, 2);
% Compute the discrete Fourier transform (DFT) of y with respect
to time
Fs = 1 / (time(2) - time(1)); % Sampling frequency
N = length(y1); % Length of the signal
frequencies = linspace(0, Fs/2, N/2 + 1); % Frequency axis up to
Nyquist frequency
% Compute the DFT using FFT (Fast Fourier Transform)
Y1 = fft(y1);
amplitude_spectrum1 = abs(Y1(1:N/2 + 1));
% Initialize the combined spectrum array
amplitude_spectrumtot = zeros(1, 2 * length(amplitude_spectrum1));
% Copy the first half of amplitude_spectrumtot from
amplitude_spectrum1
amplitude_spectrumtot(1:length(amplitude_spectrum1)) =
flip(amplitude_spectrum1);
% Copy the second half of amplitude_spectrumtot from
amplitude_spectrum
amplitude_spectrumtot(length(amplitude_spectrum1) + 1:end) =
amplitude_spectrum1;
% Initialize the combined freq array
freq_spectrumtot = zeros(1, 2 * length(amplitude_spectrum1));
% Copy the first half of amplitude_spectrumtot from
amplitude_spectrum1
freq_spectrumtot(1:length(frequencies)) = -flip(frequencies);
% Copy the second half of amplitude_spectrumtot from
amplitude_spectrum
freq_spectrumtot(length(frequencies) + 1:end) = frequencies;
% Display the amplitude spectrum
figure;
plot(freq_spectrumtot, amplitude_spectrumtot);
5
title('Combined Amplitude Spectrum O_{11}/O_{12}');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
hold
y2 = data(:, 3);
% Compute the discrete Fourier transform (DFT) of y with respect
to time
Fs = 1 / (time(2) - time(1)); % Sampling frequency
N = length(y2); % Length of the signal
frequencies = linspace(0, Fs/2, N/2 + 1); % Frequency axis up to
Nyquist frequency
% Compute the DFT using FFT (Fast Fourier Transform)
Y2 = fft(y2);
amplitude_spectrum2 = abs(Y2(1:N/2 + 1));
% Initialize the combined spectrum array
amplitude_spectrumtot2 = zeros(1, 2 *
length(amplitude_spectrum2));
% Copy the first half of amplitude_spectrumtot from
amplitude_spectrum1
amplitude_spectrumtot2(1:length(amplitude_spectrum2)) =
flip(amplitude_spectrum2);
% Copy the second half of amplitude_spectrumtot from
amplitude_spectrum
amplitude_spectrumtot2(length(amplitude_spectrum2) + 1:end) =
amplitude_spectrum2;
% Display the amplitude spectrum
plot(freq_spectrumtot, amplitude_spectrumtot2);
title('Combined Amplitude Spectrum O_{21}/O_{22}');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
% frequency response
i=sqrt(-1);
vett_f=0:0.01:24;
vect_F0=zeros(ntot-ndgl,1);
for k=1:length(frequencies)
6
ome=frequencies(k)*2*pi;
vect_F0(2)=Y1(k);
vect_F0(6)=Y1(k);
vect_F0(4)=Y2(k);
vect_F0(8)=Y2(k);
vect_F=-(-ome^2*MFC+i*ome*CFC+KFC)*vect_F0;
A=-ome^2*MFF+i*ome*CFF+KFF;
x=A\vect_F;
y7=x(a);
acc_nodo7=-ome^2*y7;
displ(k)=abs(y7);
acc(k)=abs(acc_nodo7);
ya(k,1)=y7;
dd_ya(k,1)=acc_nodo7;
end
% Display the displacement amplitude spectrum
figure;
plot(frequencies, displ);
title(["vertical displacement point A due to O_{12}/O_{12} andO_{21}/O_{22}"]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
% Display the acceleration amplitude spectrum
figure;
plot(frequencies, acc);
hold
title(["blu = vertical acceleration point A due to O_{12}/O_{12} and O_{21}/O_{22}"]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
7
% plot ya in time domain
ya_full = cat (1, ya, zeros (N/2+1-length (ya), 1)) ; %concatenare
zeri dopo ya per avere stesso numero di campioni dell' input
ya_tot = cat (1, ya_full, flip (conj (ya_full (2:end-1)))) ; %per
essere fft devo avere vettore speculare
% inverse DFT
xj = ifft (ya_tot); %dovrebbe tornare numeri reali ma torna numeri
complessi
figure;
plot (time,xj)
xlabel ('time (s]');
ylabel ('displacement [m]');
title('Vertical displacement due to earthquake excitation');
% plot dd_ya in time domain
dd_ya_full = cat (1, dd_ya, zeros (N/2+1-length (dd_ya), 1)) ;
%concatenare zeri dopo ya per avere stesso numero di campioni
dell' input
dd_ya_tot = cat (1, dd_ya_full, flip (conj (dd_ya_full
(2:end-1)))) ; %per essere fft devo avere vettore speculare
% inverse DFT
xj = ifft (dd_ya_tot); %dovrebbe tornare numeri reali ma torna
numeri complessi
figure;
plot (time,xj)
xlabel ('time (s]');
ylabel ('acceleration [m/s^2]');
title('Vertical acceleration due to earthquake excitation');
8
%------------------------------------
% point 8 - FRF of the vertical acceleration of point B (modified
system)
load("ponte_dampe_mkr");
i=sqrt(-1);
vett_f=0:0.01:24;
b=idb(9,2);
a=idb(7,2);
vect_F0=zeros(ndgl,1);
vect_F0(b)=1;
for k=1:length(vett_f)
ome=vett_f(k)*2*pi;
A=-ome^2*MFF+i*ome*CFF+KFF;
x=A\vect_F0;
y9=x(b);
acc_nodo9=-ome^2*y9;
mod6(k)=abs(acc_nodo9);
fas6(k)=angle(acc_nodo9);
end
figure
subplot 211;plot(vett_f,mod6);grid;
hold
xlabel('Freq. [Hz]');
subplot 212;plot(vett_f,fas6);grid
hold
xlabel('Freq. [Hz]');
subplot 211;plot(vett_f,mod2);grid on;
hold
title('Acceleration node B');
xlabel('Freq. [Hz]');
subplot 212;plot(vett_f,fas2);grid on;
xlabel('Freq. [Hz]');