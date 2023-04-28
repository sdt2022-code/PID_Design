%% Open-Loop Transfer Function
% Parameters
K = 1.4*10^-4; 
z1 = 1.5*10^-3; z2 = 3.6*10^-5; 
p1 = 2.87*10^-2; p2 = 7.75*10^-3; p3 = 2.85*10^-4; p4 = 8*10^-5; 
Td = 8.5; 
% Transfer Function 
s = tf('s'); 
Lag = exp(-Td*s); 
H1 = K*((s+z1)*(s+z2))/((s+p1)*(s+p2)*(s+p3)*(s+p4)); 
G1 = Lag*H1; 
figure(1); % Exact Bode Plot
bode(G1);
title("Exact Bode Plot"); 
grid on; 


%% System ID Using Frequency Response Estimation - Results
% This part of the code extracts the results of the Frequency Response
% Estimation done in Simulink and plots them in Matlab. The code won't work
% if the simulator is not opened, so this section is commented. 

% freq = H_exp.Frequency; % Retrieving Test Input Frequency Range
% resp = H_exp.ResponseData; % Retrieving Output Response Data as an Array of Complex Numbers
% resp = permute(resp,[3 1 2]); % Converting Output Data Into 1D Vector
% mag = 20*log10(abs(resp)); phase = angle(resp)*180/pi; 
% figure(2); 
% subplot(2,1,1); 
% semilogx(freq,mag);
% title("Experimental Bode Magnitude Plot of G(s)"); 
% xlabel("Frequency (rad/s)"); 
% ylabel("Magnitude (dB)"); 
% xlim([1e-06 1]); 
% ylim([-80 20]); 
% grid on; 
% % axis([1e-06 1 -80 80]);
% subplot(2,1,2); 
% semilogx(freq,phase); 
% title("Experimental Bode Phase Plot of G(s)"); 
% xlabel("Frequency (rad/s)"); 
% ylabel("Phase(degrees)"); 
% grid on; 

%% Validation of Identified System 
K_id = 1.28*10^-4; 
z1_id = 3*10^-5; z2_id = 1.5*10^-3;
p1_id = 8*10^-5; p2_id = 2.4*10^-4; p3_id = 7*10^-3; p4_id = 2.6*10^-2; 
Td_id = 8.51; 

% Identified Transfer Function 
Lag_id = exp(-Td*s); 
H1_id = K_id*((s+z1_id)*(s+z2_id))/((s+p1_id)*(s+p2_id)*(s+p3_id)*(s+p4_id)); 
G1_id = Lag_id*H1_id; 
figure(3);
subplot(2,2,[1 3]); 
bode(G1); 
title('Exact Transfer Function'); 
grid on; 
subplot(2,2,[2 4]); 
bode(G1_id,'r'); 
title('Identified Transfer Function'); 
grid on; 

%% Pade Approximation 
Lag_pade = pade(Lag,1); 
G1_pade = Lag_pade*H1; 
figure(4); 
bode(G1_pade);
title("Pade Approximation Bode Plot"); 
grid on; 

%% System ID of Pade Approximation Manual Approach - Results 
% Extract Data from Excel File 
data = xlsread('pade_frequency_response.xlsx');
freq = data(:,1); 
mag = data(:,2); 
phase = data(:,3); 
figure(5);  
subplot(2,1,1); 
semilogx(freq,mag,'*-');
title("Experimental Bode Magnitude Plot of Gpade(s)"); 
xlabel("Frequency (rad/s)"); 
ylabel("Magnitude (dB)"); 
grid on; 
subplot(2,1,2); 
semilogx(freq,phase,'*-'); 
title("Experimental Bode Phase Plot of Gpade(s)"); 
xlabel("Frequency (rad/s)"); 
ylabel("Phase(degrees)"); 
axis([1e-06 1 -400 100]);
grid on; 

%% Higher Order Pade Approximations
Lag_pade2 = pade(Lag,2); 
G1_pade2 = K*Lag_pade2*((s+z1)*(s+z2))/((s+p1)*(s+p2)*(s+p3)*(s+p4)); 
Lag_pade3 = pade(Lag,3);
G1_pade3 = K*Lag_pade3*((s+z1)*(s+z2))/((s+p1)*(s+p2)*(s+p3)*(s+p4));
% Comparing Bode Plots
figure(6);
subplot(2,2,1);
bode(G1); 
title('Exact');
axis([1e-06 10 -900 180]);
subplot(2,2,2);
bode(G1_pade); 
title('1st Order');
axis([1e-06 10 -540 540]);
subplot(2,2,3);
bode(G1_pade2);
title('2nd Order');
axis([1e-06 10 -540 540]);
subplot(2,2,4);
bode(G1_pade3);
title('3rd Order');
axis([1e-06 10 -180 900]);
% Comparing Root Locus Plots
figure(7);
subplot(2,2,1);
rlocus(G1_pade);
title('1st Order'); 
subplot(2,2,2);
rlocus(G1_pade2);
title('2nd Order'); 
subplot(2,2,3);
rlocus(G1_pade3); 
title('3rd Order'); 


%% Linear Approximation of G
[num, den] = linmod('G_original');
G1_linear = tf(num,den); 
figure(8); 
subplot(2,2,1);
bode(G1); 
title('Exact');
subplot(2,2,2);
bode(G1_linear); 
title('Linmod Approximation');
subplot(2,2,3);
bode(G1_pade); 
title('1st Order Pade Approximation');
subplot(2,2,4);
bode(G1_pade2); 
title('2nd Order Pade Approximation');



