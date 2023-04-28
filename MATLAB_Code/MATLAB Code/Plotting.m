% This script is for plotting data extracted from the simulink scope. 
figure('DefaultAxesFontSize',18); 
plot(non_linear_resp);
title("Effect of Noise on Step Response", 'FontSize', 26);
xlabel("Time (s)",'FontSize', 22); 
ylabel("Amplitude",'FontSize', 22); 
legend('Lag-Lead','PID (without filter)','PID (with filter)','FontSize', 18); 
axis([0 7000 0 1.5]); 
