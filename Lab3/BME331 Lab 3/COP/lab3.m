dz=0.043;
quietopen = importdata('Eye open standing.txt');
quietopen_interpolated = interp1(0:1/200:120-1/200, quietopen, 0:1/960:120-1/960);

Mx1=quietopen_interpolated(:,4);
Fx1=quietopen_interpolated(:,1);
My1=quietopen_interpolated(:,5);
Fy1=quietopen_interpolated(:,2);
Fz1=quietopen_interpolated(:,3);

x1=-(My1 + Fx1*dz)./Fz1;
y1=(Mx1-Fy1*dz)./Fz1;
x1(isnan(x1))=0;
y1(isnan(y1))=0;

[b a] = butter(4, 5/480, 'low'); %the first parameter is the filter order, and the second parameter is
%the cutoff frequency divided by half the sampling frequency.
x1_filtered = filtfilt(b, a, x1); %x and y are your original and filtered signals, respectively, to be
%renamed as appropriate
x1_filtered = x1_filtered - mean(x1_filtered); %make the signal have zero-mean

y1_filtered = filtfilt(b, a, y1);
y1_filtered = y1_filtered - mean(y1_filtered);


%%----------------------------------------------------------------------%%

quietclosed = importdata('Closed eye standing.txt');
quietclosed_interpolated = interp1(0:1/200:120-1/200, quietclosed, 0:1/960:120-1/960);
Mx2=quietclosed_interpolated(:,4);
Fx2=quietclosed_interpolated(:,1);
My2=quietclosed_interpolated(:,5);
Fy2=quietclosed_interpolated(:,2);
Fz2=quietclosed_interpolated(:,3);

x2=-(My2 + Fx2*dz)./Fz2;
y2=(Mx2-Fy2*dz)./Fz2;
x2(isnan(x2))=0;
y2(isnan(y2))=0;

[b a] = butter(4, 5/480, 'low'); %the first parameter is the filter order, and the second parameter is
%the cutoff frequency divided by half the sampling frequency.
x2_filtered = filtfilt(b, a, x2); %x and y are your original and filtered signals, respectively, to be
%renamed as appropriate
x2_filtered = x2_filtered - mean(x2_filtered); %make the signal have zero-mean
y2_filtered = filtfilt(b, a, y2);
y2_filtered = y2_filtered - mean(y2_filtered);

%%----------------------------------------------------------------------%%

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

%%----------------------------------------------------------------------%%
n= length(x1_filtered);
MDISTx1 = (1/n)*sum(abs(x1_filtered));
MDISTx2 = (1/n)*sum(abs(x2_filtered));
MDISTx3 = (1/n)*sum(abs(x3_filtered));
MDISTy1 = (1/n)*sum(abs(y1_filtered));
MDISTy2 = (1/n)*sum(abs(y2_filtered));
MDISTy3 = (1/n)*sum(abs(y3_filtered));

%%----------------------------------------------------------------------%%
RDISTx1 = sqrt( (1/n)* sum((x1_filtered).^2));
RDISTx2 = sqrt( (1/n)* sum((x2_filtered).^2));
RDISTx3 = sqrt( (1/n)* sum((x3_filtered).^2));
RDISTy1 = sqrt( (1/n)* sum((y1_filtered).^2));
RDISTy2 = sqrt( (1/n)* sum((y2_filtered).^2));
RDISTy3 = sqrt( (1/n)* sum((y3_filtered).^2));

%%----------------------------------------------------------------------%%
RANGEx1= max(x1_filtered) - min(x1_filtered);
RANGEx2= max(x2_filtered) - min(x2_filtered);
RANGEx3= max(x3_filtered) - min(x3_filtered);
RANGEy1= max(y1_filtered) - min(y1_filtered);
RANGEy2= max(y2_filtered) - min(y2_filtered);
RANGEy3= max(y3_filtered) - min(y3_filtered);

%%----------------------------------------------------------------------%%
MVELOx1=0;
 for i=1:n-1
MVELOx1 = (MVELOx1+abs(x1_filtered(i+1)-x1_filtered(i)));
 end
MVELOx1=MVELOx1/120;

 MVELOx2=0;
 for i=1:n-1
MVELOx2 = (MVELOx2+abs(x2_filtered(i+1)-x2_filtered(i)));
 end
 MVELOx2=MVELOx2/120;
 
 MVELOx3=0;
 for i=1:n-1
MVELOx3 = (MVELOx3+abs(x3_filtered(i+1)-x3_filtered(i)));
 end
 MVELOx3=MVELOx3/120;
 
MVELOy1=0;
 for i=1:n-1
MVELOy1 = (MVELOy1+abs(y1_filtered(i+1)-y1_filtered(i)));
 end
MVELOy1=MVELOy1/120;

 MVELOy2=0;
 for i=1:n-1
MVELOy2 = (MVELOy2+abs(y2_filtered(i+1)-y2_filtered(i)));
 end
MVELOy2=MVELOy2/120;

 MVELOy3=0;
 for i=1:n-1
MVELOy3 = (MVELOy3+abs(y3_filtered(i+1)-y3_filtered(i)));
 end
MVELOy3=MVELOy3/120;
%%----------------------------------------------------------------------%%

MVELOx1;
MVELOy1;
MVELOx2;
MVELOy2;
MVELOx3;
MVELOy3;

quietopen_emg = importdata('standingeyeopen.data');
%plot(quietopen_emg)
TA1=quietopen_emg(:,1);
SOL1=quietopen_emg(:,2);
[b a] = butter(4, 5/480, 'low'); %the first parameter is the filter order, and the second parameter is
%the cutoff frequency divided by half the sampling frequency.
TA1_filtered = filtfilt(b, a, TA1);%x and y are your original and filtered signals, respectively, to
%be renamed as appropriate
SOL1_filtered = filtfilt(b, a, SOL1);
%plot(TA1_filtered)
% hold on
%plot(SOL1_filtered)


%%----------------------------------------------------------------------%%
quietclosed_emg = importdata('standingeyeclose.data');
%plot(quietclosed_emg)
TA2=quietclosed_emg(:,1);
SOL2=quietclosed_emg(:,2);
[b a] = butter(4, 5/480, 'low');
TA2_filtered = filtfilt(b, a, TA2);
SOL2_filtered = filtfilt(b, a, SOL2);
%plot(TA2_filtered)
% hold on
plot(SOL2_filtered)
%title('Eyes Closed Standing Experiment EMG - SOL Data');

%%----------------------------------------------------------------------%%
freq = 960;
seg_length = 8192;
for k = 1:13
 a = (k-1)*(seg_length)+ 4801;
 b = (seg_length+4800) + (k-1)*(seg_length);
 ccf3 = xcorr(SOL2_filtered(a:b) - mean(SOL2_filtered(a:b)), y2_filtered(a:b) - mean(y2_filtered(a:b)),'coeff');
ccf_seg3(k,:) = ccf3;
end
avg_ccf3 = mean(ccf_seg3);
time = (-(seg_length-1):1:(seg_length-1))/freq;
CVAL = max(avg_ccf3);
index = find(avg_ccf3 ==max(avg_ccf3(:)));
time_shift = time(index);
%
plot(time,avg_ccf3);
xlabel('Time Shift (s)');
title('CCF Between AP and SOL in Eyes Closed Standing Experiment');

CVAL;
time_shift;

%%----------------------------------------------------------------------%%

freq = 960;
seg_length = 8192;
for k = 1:13
 a = (k-1)*(seg_length)+ 4801;
 b = (seg_length+4800) + (k-1)*(seg_length);
 ccf4 = xcorr(TA2_filtered(a:b) - mean(TA2_filtered(a:b)), y2_filtered(a:b) - mean(y2_filtered(a:b)),'coeff');
 ccf_seg4(k,:) = ccf4;
end
avg_ccf4 = mean(ccf_seg4);
time = (-(seg_length-1):1:(seg_length-1))/freq;
plot(time,avg_ccf4);
legend('TA - quiet closed');

CVAL = max(avg_ccf4);
index = find(avg_ccf4 ==max(avg_ccf4(:)));
time_shift = time(index);

plot(time,avg_ccf4);
xlabel('Time Shift (s)');
title('CCF Between AP and TA in Eyes Closed Standing Experiment');

%%----------------------------------------------------------------------%%

freq = 960;
seg_length = 8192;
for k = 1:13
 a = (k-1)*(seg_length)+ 4801;
 b = (seg_length+4800) + (k-1)*(seg_length);
 ccf5 = xcorr(SOL1_filtered(a:b) - mean(SOL1_filtered(a:b)), y1_filtered(a:b) -mean(y1_filtered(a:b)),'coeff');
 ccf_seg5(k,:) = ccf5;
end
avg_ccf5 = mean(ccf_seg5);
time = (-(seg_length-1):1:(seg_length-1))/freq;
plot(time,avg_ccf5);
legend('SOL - quiet open');
CVAL = max(avg_ccf5);
index = find(avg_ccf5 ==max(avg_ccf5(:)));

time_shift = time(index);

xlabel('Time Shift (s)');
title('CCF Between AP and SOL in Eyes Open Standing Experiment');

CVAL;
time_shift;

%%----------------------------------------------------------------------%%
freq = 960;
seg_length = 8192;
for k = 1:13
 a = (k-1)*(seg_length)+ 4801;
 b = (seg_length+4800) + (k-1)*(seg_length);
 ccf6 = xcorr(TA1_filtered(a:b) - mean(TA1_filtered(a:b)), y1_filtered(a:b) -mean(y1_filtered(a:b)),'coeff');
 ccf_seg6(k,:) = ccf6;
 end
avg_ccf6 = mean(ccf_seg6);
time = (-(seg_length-1):1:(seg_length-1))/freq;

plot(time,avg_ccf6);
CVAL = max(avg_ccf6);
index = find(avg_ccf6 ==max(avg_ccf6(:)));
time_shift = time(index);
xlabel('Time Shift (s)');
title('CCF Between AP and TA in Eyes Open Standing Experiment');

CVAL;
time_shift;

%%----------------------------------------------------------------------%%

perturbed_emg = importdata('standingpulling.data');
TA3=perturbed_emg(:,1);
SOL3=perturbed_emg(:,2);
[b a] = butter(4, 5/480, 'low'); %the first parameter is the filter order, and the second parameter is
%the cutoff frequency divided by half the sampling frequency.
TA3_filtered = filtfilt(b, a, TA3);%x and y are your original and filtered signals, respectively, to
%be renamed as appropriate
SOL3_filtered = filtfilt(b, a, SOL3);
plot(TA3_filtered)
legend('TA - Perturbed')
max(TA3_filtered)
hold on

perturbed_emg = importdata('standingpulling.data');
TA3=perturbed_emg(:,1);
SOL3=perturbed_emg(:,2);
[b a] = butter(4, 5/480, 'low'); %the first parameter is the filter order, and the second parameter is
%the cutoff frequency divided by half the sampling frequency.
TA3_filtered = filtfilt(b, a, TA3);%x and y are your original and filtered signals, respectively, to
%be renamed as appropriate
SOL3_filtered = filtfilt(b, a, SOL3);
plot(SOL3_filtered)
legend('SOL - Perturbed')
%max(TA3_filtered)
hold off

hold on
plot(y3_filtered)
legend('TA', 'SOL', 'AP')
hold off






