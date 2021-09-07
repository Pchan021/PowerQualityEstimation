%POWER QUALITY MEASURE

C_values=RC_load1;
V_values=RCV_load2;

%loading data
%taking phaseA values
Curr=C_values(:,2);%column1 is time column2 is phaseA values column3 is phaseB...
Volt=V_values(:,2);
%data set up
number_of_readings = length(Curr);
stop_time = C_values(number_of_readings,1);%time stamp of the last reading
%situation setup
n=100; % measuring up 'n' number of harmonics
time_step = stop_time/number_of_readings;
time_clock = 0:time_step:stop_time;
freq=50;

%Processing setup
%generating sine and cosine waves 
for i=1:n %we are generating waves upto the nth harmonic
    for j= 1:number_of_readings
        %%column matrix with harmonics(i) along x  axis and values(j) in y axis
        sine_wave(j,i) = sin(i * freq*2*pi* time_clock(j));% of the form sin(k * 2pi/L *x) 
        cos_wave(j,i) = cos(i * freq*2*pi* time_clock(j));
    end;
end;

%feature Extraction process
sine_part_I=0; cos_part_I=0;
sine_part_V=0; cos_part_V=0;
for i=1:n
    for j=1:number_of_readings
        sine_part_I(j,i)= Curr(j)*(sine_wave(j,i));
        cos_part_I(j,i)= Curr(j)*(cos_wave(j,i));
        sine_part_V(j,i)= Volt(j)*(sine_wave(j,i));
        cos_part_V(j,i)= Volt(j)*(cos_wave(j,i));
        
    end
end

%developing the component waves
sum_sq_I=0;sum_sq_V=0;
for i=1:n
    a=0;b=0;c=0;d=0;
    for j=1:number_of_readings %The integration part
        a=a+ sine_part_I(j,i);
        b=b+ cos_part_I(j,i);
        c=c+ sine_part_V(j,i);
        d=d+ cos_part_V(j,i);
    end
    %%Fourier Coefficients Ak and Bk
    Bk_for_I= (2/number_of_readings)*a; %2/L is normalizing
    Ak_for_I= (2/number_of_readings)*b;
    Bk_for_V= (2/number_of_readings)*c;
    Ak_for_V= (2/number_of_readings)*d;
    
    I_mag(i)= sqrt(Bk_for_I^2 + Ak_for_I^2); %current magnitude of I'th harmonic
    I_arg(i)= atan(Ak_for_I/Bk_for_I) * 180/pi; % current angle of I'the harmonic
    V_mag(i)= sqrt(Bk_for_V^2 + Ak_for_V^2); %voltage magnitude of I'th harmomic
    V_arg(i)= atan(Ak_for_V/Bk_for_V) * 180/pi; %voltage angle of V'th harmonic
    
    sum_sq_I= sum_sq_I+I_mag(i)^2; % for rms calculation
    sum_sq_V= sum_sq_V+V_mag(i)^2;
end

%Calculation finding RMS, THD and PowerFactor
I_rms= sqrt(sum_sq_I /2); %rms of total current waveform
I_fundamental= I_mag(1)/sqrt(2); %rms value of fundamental component of I
THD_I= 100*(sqrt(I_rms^2 - I_fundamental^2)) /I_fundamental %Total Harmonic Distortion
CF_I= max(Curr)/I_rms %Crest Factor

V_rms= sqrt(sum_sq_V /2);
V_fundamental= V_mag(1)/sqrt(2);
THD_V= 100*(sqrt(V_rms^2 - V_fundamental^2)) /V_fundamental%Total harmonic Distortion
CF_V= max(Volt)/V_rms %Crest Factor

disp_pf= cosd(I_arg(1)-V_arg(1)); % pf between fundamental components
dist_pf= 1/(sqrt(1+(THD_I/100)^2)*sqrt(1+(THD_V/100)^2));
pf = disp_pf * dist_pf %true pf = distortion pf * displacement pf


%%Grephs
% FFT spectrum
t2=1:n;
figure(1)
bar(t2,I_mag)
xlabel('Harmonic Number')
ylabel('Magnitude')
title('Current Harmonic spectrum')

figure(2)
bar(t2,V_mag)
xlabel('Harmonic Number')
ylabel('Magnitude')
title('Voltage Harmonic spectrum')


%%----------------------GRADING-----------------------------%%


%%Power Quality Measures in per-unit
PQ(1)= THD_I /100;
PQ(2)= THD_V /100;
PQ(3)= sqrt(1- pf^2); %% Orthogonal Current Factor == sin theta


        %%weights
%%K1 & K2 for THD_I and THD_V
K(1)=1/3; K(2)=1/3;
%%K3 for OCF 
K(3)=1/3;

PQF=0; %%powerQUalityFactor
for i=1:3
   PQF = PQF + K(i)*(1 - PQ(i));
end
 
PQF=PQF/1