%% Question 1
clear; clc; close;
%h(t)=g(Tb-t) because receiver is a matched filter
%1-First Question :
%a)
%assume t=x
%First figure :
figure()

Tb=2
%will draw 2*Tb*fs points ... 200 points Tb*fs and Tb*fs points 1
Ts=1/200
fs = 1/Ts

%Tb/Ts=Tb*fs
t = linspace(0, Tb, Tb*fs);
g = [2* t(1:(length(t)/2)) /Tb, ones(1, Tb*fs/2)];

plot(t,g)
title('Input signal & Impulse response  ')
ylabel('g(t)')
legend({'g(t)'},'Location','southwest')
xlabel('time')
hold off
%---------------------------------
%Plotting h(t) only
figure()
h=fliplr(g)
plot(t,h) %note h(t)=k*g(t-Tb))---assume k = 1
title('Impulse response h(t) ')
ylabel('h(t)')
legend({'h(t)'},'Location','southwest')
xlabel('time')
hold off
%----------------------------
%b)
%If the input array is [-1 1 -1 -1 1] so g1 is [-g g -g -g g]
%Plotting g1(t): g1 is train of pulses 
figure()
t1=linspace(0,5*Tb,5*(Tb/Ts))
Input=[-1 1 -1 -1 1]
%empty matrix 
%g1=
%length(g1)=5*length(g)
k=1;

for i= 1 :length(g):5*length(g)
    g1(i:i+length(g)-1)=Input(k)*g
    k=k+1
end   
% adding noise
% uncomment in the presence of awgn
%g1 = 0.5* awgn(g1,5);

plot(t1,g1)
title('Input signal(Pulses) g1(t)')
ylabel('g1(t)')
legend({'g1(t)'},'Location','southwest')
xlabel('time')
hold off

%----------
%covolution 
%h(t)*g(t)*0.5=g0(t)=y(t) 
%L : length of time should = [ length(g1)+length(h)-1 ]-convolution 
L=(length(g1)+length(h)-1)
figure()
t2=linspace(0,(length(Input)+1)*Tb,L)
y=Ts*0.5*conv(h,g1)

%uncomment in the presence of awgn
%y=Ts*conv(h,g1)


plot(t2,y)
title('Convolution Output of h(t) ')
ylabel('y(t)')
legend({'y(t)'},'Location','southwest')
xlabel('time')
hold off

%Sampling :

figure()
%%%%%%%%%%%%%%%%%%
%Tb/Ts %L 
for i=Tb/Ts:Tb/Ts:L
    stem(t2(i),y(i))
    hold on
end   
title('sampled output of the optimum receiver, h1(t)')
ylabel('y1(t)')
legend({'y1(t)'},'Location','southwest')
xlabel('time')
hold off

size(g1)
size(t1)

figure()

t4 = linspace(0, Tb, Tb*fs);
h2 = [zeros(1, Tb*fs/2),ones(1, Tb*fs/2)];

plot(t4,h2) 
title('Impulse response h2(t) ')
ylabel('h2(t)')
legend({'h2(t)'},'Location','southwest')
xlabel('time')
hold off

%covolution 
%h2(t)*g(t)*0.5=g0(t)=y3(t) 
%L : length of time should = [ length(g1)+length(h2)-1 ]-convolution 
L=(length(g1)+length(h2)-1)
figure()
t5=linspace(0,(length(Input)+1)*Tb,L)
y3=Ts*0.5*conv(h2,g1)
plot(t5,y3)
title('Convolution Output of h2(t)')
ylabel('y2(t)')
legend({'y2(t)'},'Location','southwest')
xlabel('time')
hold off

%Sampling :

figure()
%%%%%%%%%%%%%%%%%%
%Tb/Ts %L 
for i=Tb/Ts:Tb/Ts:L
    stem(t5(i),y3(i))
    hold on
end   

title('sampled output of the receiver filter, h2(t)')
ylabel('y2(t)')
legend({'y2(t)'},'Location','southwest')
xlabel('time')
hold off

%% Question 3
clear; clc; close;
%below we ask the user for all the parameters needed before we proceed
fprintf('Please choose the pulse shape:\n');
fprintf('1. NRZ\n2. Ideal Nyquist Pulse\n3. Raised-Cosine\n4. Square Root Raised-Cosine\n');
pulse_type = input('Enter the number of the chosen type: ');
pam_level = input('Enter numebr of PAM levels(2 or 8): ');
duration = input('Enter 1 for eye digram of duration Ts or 2 for duration 2Ts: ');

%generate 300 random bits uniformly distributed(equiprobable bits)
no_bits = 300; %number of bits to generate
bits = randi([0, 1], 1, no_bits);
Tb = 1; %each pulse is 1 sec
fs = 1000; %sampling freq used to plot pulses

switch pulse_type %decide which pulse type to use based on the user's choice
    case 1 %polar
        encoded = PolarPAM(bits, pam_level);
        %repeat each value in the encoded array Ts*fs times to make it a pulse
        encoded_pulses = repelem(encoded, Tb*fs);
        %define the time axis for the desired time(either Ts or 2Ts)
        t = linspace(0, duration*Tb, duration*Tb*fs);
        %draw eye diagram by taking one/two pulse each iteration and imposing it on the plot
        for i = 1 : duration*Tb*fs : length(encoded_pulses)
            encoded_pulses(1, i) = 0;
            encoded_pulses(1, i + duration*Tb*fs - 1) = 0;
            plot(t, encoded_pulses(1, i : i + duration*Tb*fs - 1))
            hold on
        end
        grid on; xlim([-0.1 duration*Tb + 0.1]);
        title(sprintf('Eye Diagram of %d level Polar PAM with duaration %dTs', pam_level, duration))
        xlabel('time s')
        ylabel('Amplitude V')
    case 2 %Ideal nyquist
        %define the time axis on which we want to generate the sinc signal
        t_sinc = linspace(-5, 5, Tb*1001*10);
        %generate the sinc signal using the time axis defined above
        raised_cos = sinc(t_sinc);
        %apply polar pam to the input bits to transform the input into correct levels
        input = PolarPAM(bits, pam_level);
        %repeat each value in the encoded array Tb*1001 times to make it a pulse of same width as sinc
        input_pulses = repelem(input, Tb*1001);
        %we want to turn the pulses into just an impulse in the middle of the pulse
        %generate array of zeros same length as the input_pulses
        temp = zeros(1, length(input_pulses));
        %only put 1 in the middle of each pulse
        temp(501:1001:end) = 1;
        %finally to get the impulses multiply element by element in temp and input_pulses
        impulses = temp.*input_pulses;
        %convolve sinc signal and input impulses together to obtain pam modulation
        result = conv(raised_cos, impulses, 'same');
        t = linspace(0, duration*Tb, duration*Tb*1001);
        %draw eye diagram by taking one/two pulse each iteration and imposing it on the plot
        for i = 1 : duration*Tb*1001 : length(result)
            plot(t, result(1, i : i + duration*Tb*1001 - 1))
            hold on
        end
        grid on; xlim([-0.1 duration*Tb + 0.1]);
        title(sprintf('Eye Diagram of %d level Ideal Nyquist PAM with duaration %dTs', pam_level, duration))
        xlabel('time s')
        ylabel('Amplitude V')
    case 3 %raised cos
        %define the time axis on which we want to generate the sinc signal
        t_sinc = linspace(-5, 5, Tb*1001*10);
        %define the roll off factor
        alpha = 1;
        %generate the raised cos signal using the time axis defined above
        raised_cos = sinc(t_sinc).*cos(2*pi*alpha*t_sinc)./(1-16*(alpha^2)*(t_sinc.^2));
        %apply polar pam to the input bits to transform the input into correct levels
        input = PolarPAM(bits, pam_level);
        %repeat each value in the encoded array Tb*1001 times to make it a pulse of same width as sinc
        input_pulses = repelem(input, Tb*1001);
        %we want to turn the pulses into just an impulse in the middle of the pulse
        %generate array of zeros same length as the input_pulses
        temp = zeros(1, length(input_pulses));
        %only put 1 in the middle of each pulse
        temp(501:1001:end) = 1;
        %finally to get the impulses multiply element by element in temp and input_pulses
        impulses = temp.*input_pulses;
        %convolve sinc signal and input impulses together to obtain pam modulation
        result = conv(raised_cos, impulses, 'same');
        t = linspace(0, duration*Tb, duration*Tb*1001);
        %draw eye diagram by taking one/two pulse each iteration and imposing it on the plot
        for i = 1 : duration*Tb*1001 : length(result)
            plot(t, result(1, i : i + duration*Tb*1001 - 1))
            hold on
        end
        grid on; xlim([-0.1 duration*Tb + 0.1]);
        title(sprintf('Eye Diagram of %d level Raised Cos PAM with duaration %dTs', pam_level, duration))
        xlabel('time s')
        ylabel('Amplitude V')
    case 4 %sqrt raised cos
        %define the time axis on which we want to generate the sinc signal
        t_sinc = linspace(-5, 5, Tb*1001*10);
        %define the roll off factor
        alpha = 1; Tc = 0.2604167;
        %generate the raised cos signal using the time axis defined above
        num = cos(pi*t_sinc*(1+alpha)/Tc) + (Tc./(4*alpha*t_sinc)).*sin(pi*t_sinc*(1-alpha)/Tc);
        den = 1 - (4*alpha*t_sinc/Tc).^2;
        sqrt_raised_cos = (4*alpha/(pi*sqrt(Tc)))*num./den;
        %apply polar pam to the input bits to transform the input into correct levels
        input = PolarPAM(bits, pam_level);
        %repeat each value in the encoded array Tb*1001 times to make it a pulse of same width as sinc
        input_pulses = repelem(input, Tb*1001);
        %we want to turn the pulses into just an impulse in the middle of the pulse
        %generate array of zeros same length as the input_pulses
        temp = zeros(1, length(input_pulses));
        %only put 1 in the middle of each pulse
        temp(501:1001:end) = 1;
        %finally to get the impulses multiply element by element in temp and input_pulses
        impulses = temp.*input_pulses;
        %convolve sinc signal and input impulses together to obtain pam modulation
        result = conv(sqrt_raised_cos, impulses, 'same');
        t = linspace(0, duration*Tb, duration*Tb*1001);
        %draw eye diagram by taking one/two pulse each iteration and imposing it on the plot
        for i = 1 : duration*Tb*1001 : length(result)
            plot(t, result(1, i : i + duration*Tb*1001 - 1))
            hold on
        end
        grid on; xlim([-0.1 duration*Tb + 0.1]);
        title(sprintf('Eye Diagram of %d level Sqrt Raised Cos PAM with duaration %dTs', pam_level, duration))
        xlabel('time s')
        ylabel('Amplitude V')
end

%%
function encoded = PolarPAM(bits, pam_level)
if pam_level == 2
    %in pam2 encoded array has the same size as original array
    encoded = ones(size(bits)); %create an array filled with 1s same size as original array
    encoded(bits == 0) = -1; %only set indices where the bit = 0 to the level -1
else %pam level == 8
    %in pam2 encoded array has the (1/3) size as original array as log_2(8) = 3
    encoded = zeros(1, length(bits)/3);
    encoded_index = 1; %will be used as the index of the encoded array inside the loop
    for i = 1 : 3 : length(bits) %each iteration takes three consecutive bits to encode together
        bits_to_encode = bits(1, i : i+2); %extract the bits to be encoded
        num = bin2dec(num2str(bits_to_encode)); %convert the bits to decimal
        encoded(1, encoded_index) = 2*num-7; %the formula 2*num-7 finds the correct level of each num
        encoded_index = encoded_index + 1;
    end
end
end



