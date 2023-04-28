%% Open-Loop Transfer Function
% Parameters
K = 1.7*10^-1; 
z1 = 5*10^-2; z2 = 8*10^-2; 
p1 = 10*10^-2; p2 = 3*10^-2; p3 = 4*10^-3; p4 = 0.5*10^-3; 
Td = 8.5; 

% Transfer Function 
s = tf('s'); 
Lag = exp(-Td*s);  
H2 = K*((s+z1)*(s+z2))/((s+p1)*(s+p2)*(s+p3)*(s+p4)); 
G2 = Lag*H2; 

% Pade Approximation of G2 
Lag_pade = pade(Lag,1);
G2_pade = Lag_pade*H2; 

% Open-Loop Static Position Error: Kp = lim(G2(s)) as s->0, 
Kp = evalfr(G2_pade,0); % The function evalfr(G2,s) evaluates G2 at the specified complex number s

%% Design Requirements 
ts = 1200; % 1200s or 20 minutes 
Mp = 3; 
ess = 0.02; % 0.02 or 2% for a unit step response 

% Corresponding Dominant Closed-Loop Parameters 
[zeta, wn] = SecondOrder(ts,Mp); 

% Corresponding Static Position Error; ess = 1/(1+Kp)
Kp_des = 1/ess - 1; 

% Corresponding Set Point 
SP = -zeta*wn + j*(wn*sqrt(1-zeta^2)); 

%% Angle of Deficiency 
angle_G2_pade = angle(evalfr(G2_pade,SP))*180/pi; % Evaluating the angle of G2_pade at SP 
D = 180 - angle_G2_pade; 

%% Lag-Lead Compensator Design: Kc1*(s+zc1)/(s+pc1)*(s+zc2)/(s+pc2)
% Lead Compensator Zero: found using bisection method  
theta = angle(SP)*180/pi; 
phi_zero = 180-theta/2-D/2;
zc1 = -real(SP) + imag(SP)*tand(phi_zero-90); 

% Lead Compensator Pole: found using bisection method  
phi_pole = phi_zero + abs(D);
pc1 = -real(SP) + imag(SP)*tand(phi_pole-90); 
Gamma = pc1/zc1;

% Lead Compensator 
Gc_lead = (s+zc1)/(s+pc1);

% Compensator Gain: Found Analytically
mag_G2_pade = abs(evalfr(G2_pade,SP)); % Evaluating the magnitude of G at SP. 
mag_Gc_lead = abs(evalfr(Gc_lead,SP)); % Evaluating the magnitude of Gc_lead at SP. 
Kc1 = 1/(mag_G2_pade*mag_Gc_lead); % Computing the gain from the magnitude condition 


% Lag Compensator 
Beta = (Kp_des*Gamma)/(Kc1*Kp); % Found from static position error requirement

% We must now specify the zero and pole so that we satisfy the magnitude and angle 
% constraints. We seek the largest possible values that satisfy the constraints, 
% this is because we do not want values that are too small to be implemented by 
% physical components. This procedure is implemented in the subsequent for loop. 
% The range of values for zc2 was specified by tune through trial and
% error, then a step size of 0.00005 and tolerances for the magnitude and
% angle constraints were specified and the the program was run.

for zc2 = 0.001:-0.00001:0.0001
    pc2 = zc2/Beta; 
    Gc_lag = (s+zc2)/(s+pc2); 
    mag_Gc_lag = abs(evalfr(Gc_lag,SP));
    angle_Gc_lag = angle(evalfr(Gc_lag,SP))*180/pi;
    if (abs(1-mag_Gc_lag)< 0.1 && angle_Gc_lag > -4)
        break;
    end
end

% Closed-Loop Compensated System (using exact transfer function G2)
G_lag_lead = Kc1*Gc_lead*Gc_lag; 
sys_lag_lead = feedback(G_lag_lead*G2,1); 

% Step Response 
figure(1); 
step(sys_lag_lead);
title("Step Response of Lag-Lead Compensated System"); 

% Step Response Perfromance 
info1 = stepinfo(sys_lag_lead) ;
y1 = step(sys_lag_lead,100000);
[ts1, Mp1] = deal(info1.SettlingTime,info1.Overshoot);
ess1 = abs(1-y1(end));
fprintf('Lag-Lead Compensated System Performance: \n'); 
fprintf('Settling Time: %.2f (Required: %.2f)\n', ts1,ts);
fprintf('Percent Overshoot: %.2f (Required: %.2f)\n', Mp1,Mp);
fprintf('Steady-State Error: %.4f (Required: %.4f) \n\n', ess1,ess);

%% PID Compensator Deisgn: Kc2*(s+zc3)*(s+zc4)/s; 
% We select the first compensator zero zc3 such that it eliminates the
% plant pole p4.
zc3 = p4; 

% Now we find zc4 using the angle of deficiency computed before: 
alpha = D + angle(SP)*180/pi - angle(SP+zc3)*180/pi;
zc4 = -real(SP) + imag(SP)/tand(alpha); 
G_pid = ((s+zc3)*(s+zc4))/s;  

% Controller Gain : found analytically 
mag_G_PID = abs(evalfr(G_pid,SP));
Kc2 = 1/(mag_G_PID*mag_G2_pade); 

% Extract PID gains kp ki and kd for later tuning 
[num, den] = tfdata(Kc2*G_pid, 'v'); 
kprop = num(2); 
ki = num(3);
kd = num(1); 
G_PID = kprop + kd*s +ki/s; 

% We found that the system exhibits a slight overshoot beyond the required
% value. We solve this problem by slightly increasing the derivative gain
% kd. 

% Closed-Loop Compensated System 
sys_PID = feedback(G_PID*G2,1); 

% Step Response 
figure(2); 
step(sys_PID);
title("Step Response of PID Compensated System Before and After Tuning)"); 
hold on; 

% Step Response Perfromance 
info2 = stepinfo(sys_PID) ;
y2 = step(sys_PID);
[ts2, Mp2] = deal(info2.SettlingTime,info2.Overshoot);
ess2 = abs(1-y2(end));
fprintf('PID Compensated System Performance (Before Tuning): \n'); 
fprintf('Settling Time: %.2f (Required: %.2f)\n', ts2,ts);
fprintf('Percent Overshoot: %.2f (Required: %.2f)\n', Mp2,Mp);
fprintf('Steady-State Error: %.4f (Required: %.4f) \n\n', ess2,ess); 

% Tune PID parameters to achieve desired performance 
kprop = kprop; % kprop and ki unchanged.
ki = ki;
kd = 0.01375; % kd slightly increased. 

G_PID = (kd*s^2+kprop*s+ki)/s; 
sys_PID = feedback(G_PID*G2,1); 

% Re-evaluate Performance 
step(sys_PID);
legend('Before Tuning','After Tuning');

info2 = stepinfo(sys_PID) ;
y2 = step(sys_PID);
[ts2, Mp2] = deal(info2.SettlingTime,info2.Overshoot);
ess2 = abs(1-y2(end));
fprintf('PID Compensated System Performance (After Tuning): \n'); 
fprintf('Settling Time: %.2f (Required: %.2f)\n', ts2,ts);
fprintf('Percent Overshoot: %.2f (Required: %.2f)\n', Mp2,Mp);
fprintf('Steady-State Error: %.4f (Required: %.4f) \n\n', ess2,ess); 



%% Comparaison to an Equivalent Ideal Second-Order System
sys_ideal = wn^2/(s^2+2*zeta*wn*s+wn^2); 
figure(3); 
step(sys_ideal);  
hold on; 
step(sys_lag_lead);
step(sys_PID); 
title('Design Assessment Against Ideal System'); 
legend('Ideal System','Lag-Lead Compensator','PID Compensator'); 
    
%% Bode Plots and Stability Margins (Margins Displayed in Figures)
figure(4); 
subplot(2,2,[1 3]); 
margin(G_lag_lead*G2);
legend('Lag-Lead Compensated System');

subplot(2,2,[2 4]); 
margin(G_PID*G2);
legend('PID Compensated System');
