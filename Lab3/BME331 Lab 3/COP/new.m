perturbed = importdata('Pulling.txt');
perturbed_interpolated = interp1(0:1/200:120-1/200, perturbed, 0:1/960:120-1/960);

Mx3=perturbed_interpolated(:,4);
Fx3=perturbed_interpolated(:,1);
My3=perturbed_interpolated(:,5);
Fy3=perturbed_interpolated(:,2);
Fz3=perturbed_interpolated(:,3);

x3=-(My3 + Fx3*dz)./Fz3;
y3=(Mx3-Fy3*dz)./Fz3;
x3(isnan(x3))=0;
y3(isnan(y3))=0;

[b a] = butter(4, 5/480, 'low');
x3_filtered = filtfilt(b, a, x3);
x3_filtered = x3_filtered - mean(x3_filtered);
y3_filtered = filtfilt(b, a, y3);
y3_filtered = y3_filtered - mean(y3_filtered);

plot(y3_filtered)
ylabel('AP');
xlabel('Time (s)');
title('AP Data');
