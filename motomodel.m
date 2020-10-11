%
% motomodel.m
% Marcus Wilson
% 19 September 2018
% p6321
%
% This code looks to take the output of layer 5 cortical neurons from 
% NFTSim and use them in a scheme of Li et al to construct a MEP.
%
close all; %clear

% input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 100;    %how many different motounits to consider
max_input = 900;  %Controls alpha in paper.  Alpha = (max_input/min_input)/M. maximum firing rate likely to be received as input from layer 5. Standard 900
min_input = 14;  % T_min in paper.  min firing rate per second to get a motounit to fire. Standard 10 or 14 to get no activity for 0 bg
min_firing = 8;    %q_min in paper. minimum firing rate of motoneuron. Standard 8
max_firing = 300;  %Qmmax in paper. maximum firing rate of motoneuron
lambda =  2e-3;  %the lambda parameter, (s), of moezzi for describing the MUAP function
time_delay = 10e-3;  %a 10 ms time delay (moezzi) in aligning the muaps
gain=1.0;   %a multiplying gain for the input.  
mV_scaling = 3.0;   %scaling factor to get the output correct (= M_0 / min_input in paper).  Standard 3 for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%calculate the thresholds
motoindex=[1:M];   %a list of integers
a=log(max_input/min_input)/M;   %the Li parameter 'a'
thresholds=min_input*exp(a*motoindex);   %a list of thresholds

figure(5)
plot(thresholds)
xlabel('motounit index')
ylabel('firing threshold (/s)')

%calculate the gradient of the output versus input line (depends on which
%Li method
% Assume that they all end at the max input, max firing point, but start at
% threshold, 8 Hz point
gradient_for_i = (max_firing - min_firing)*(max_input - thresholds).^(-1);  %increase in y over increase in x
%only thresholds is a vector so should get a vector out here.
%make sure last one doesn't have infinite gradient
gradient_for_i(end) = gradient_for_i(end-1); 

intercept_for_i = min_firing - (thresholds.*gradient_for_i);  %calculate intercept on y axis


%read in the firing rate data from marcus_reader_sept17.m
%load '../../../neurofield_sept17/neurofield-v.0.1.4/key_output_parameters.mat'

%loading from the neurofield file would be...
pwd
load 'key_output_parameters.mat'

phi=phi*gain;  %multiply up the input by a gain factor;

%Calculate the firing rate of each index
firing_rate_for_each = get_firing_rate(phi, thresholds, gradient_for_i, intercept_for_i); 

%Now find locations of firings for each
firing_times = get_firing_times(tt,firing_rate_for_each);

figure(1)
subplot(5,1,1)
plot(tt(1:40000)*1e3-20,phi(1:40000));
xlabel('shifted time (ms)')
ylabel('input to motoneurons (/s)')
subplot(5,1,2)
plot(tt(1:40000),firing_rate_for_each(1:40000,1:10))
xlabel('time (s)')
ylabel('firing rate for first 10 motounits')
subplot(5,1,3)
plot(tt(1:40000),firing_times(1:40000,1:10))

subplot(5,1,4)
plot(tt(1:40000),   M-sum( (firing_rate_for_each(1:40000,:)<0.1) , 2)   );   %find how many are firing (by counting up at each time the number that aren't firing)
xlabel('time (s)')
ylabel('number of units firing')

%Finally we need to convolve with our function. 

%Each unit comes in more strongly. 
%Take each firing sequence (0 and 1) and multiply it by the threshold for
%that unit

threshold_matrix = ones(length(tt),1)*thresholds;  %make a matrix where each row is identical
%but each column is the threshold for that motounit

weighted_firing = firing_times.*threshold_matrix;  %Now have the ith unit being 0 or its threshold. 

%sum them all up over all units (columns)
total_impact = mV_scaling*sum(weighted_firing,2); %6 May 2019. Multiply by a mV scale factor to get y-scale right.
subplot(5,1,5)
plot(tt(1:40000),total_impact(1:40000)); 
xlabel('time (s)')
ylabel('total motooutput strength (unnormalized)')

%Carry out convolution.  p6321
deltat=tt(2)-tt(1);
t_in_muap = [-5*lambda:deltat:5*lambda];
muap = -t_in_muap.*exp( -(t_in_muap/lambda).^2 );
figure(2)
plot(t_in_muap,muap); 
xlabel('time (s)')
ylabel('muap (unnormalized)')

motor_output = conv(total_impact,muap,'same');  %This will have a length longer than both total_impact and muap. We want the central bit
%take only middle bit for accurate timing


%shift this by 10 ms rightwards
motor_output_final = [ zeros(round(time_delay/deltat),1) ; motor_output(1:end - round(time_delay/deltat)) ];  %move by time_delay/deltat steps

figure(3)
plot(tt(1:100000),motor_output_final(1:100000))
xlabel('time (s)')
ylabel('motor output (unnormalized)')

%26 September 2018
%Do an 'average' output over meps 3 - 10...



figure(7)
%identify when we have a 1 and plot it
list=[];
for timeindex=[30000:40000]
    for motoindex=[1:100]
        if (firing_times(timeindex,motoindex)==1)
           list = [list;  (timeindex-30000)/10-20 motoindex 1];
        end
    end
end

if (size(list,1)>0);   %plot if there is something there.
     plot3(list(:,1),list(:,2),list(:,3),'k.')
     xlabel('time (ms)'); grid on
     ylabel('motounit')
     set(gca,'ylim', [0 100]); 
     set(gca,'xlim',[0 380])
     view(0,90)
end

%reshape the motor_output_final to allow an average
motor_shaped = reshape(motor_output_final,1/deltat,length(motor_output_final)*deltat); 
ave_mep = mean(motor_shaped(:,3:10), 2);
figure(6)
subplot(2,1,2)
%plot out the mep again in this window
plot(tt(1:4000)*1e3-20,ave_mep(1:4000),'k'); grid on
xlabel('time (ms)')
ylabel('motor output (unnormalized)')
title('(b)')
set(gca,'xlim',[0 380])
set(gca,'ylim',[-0.8 1.1])


figure(4)
plot(tt(1:10000)*1e3-20,ave_mep)

disp(['max in figure4 = ' num2str(max(ave_mep))])

%25 june 2019
%spit out the mep to a file
save 'mep_file.mat' ave_mep


figure(8) %This one for the paper. 31 may 2019
subplot(3,1,1)
%first plot the MEP
plot(tt(1:10000)*1e3-20,ave_mep,'k'); grid on;
set(gca,'xlim',[0 350]);
xlabel('time (ms)')
ylabel('MEP (mV)')
title('(a)')


%second plot the layer 5 rate.
%first do an average over several, like with the MEP itself
subplot(3,1,2)
phi_shaped = reshape(phi,1/deltat,length(phi)*deltat); 
ave_layer_5 = mean(phi_shaped(:,3:10), 2);
plot(tt(1:10000)*1e3-20,ave_layer_5(1:10000),'k');  grid on;
set(gca,'xlim',[0 350]);
xlabel('time (ms)')
ylabel('layer 5 pulse rate (/s)')
title('(b)')

%This one is the motounit firings.  Take the 3rd one (probably phases have
%been a bit more randomized by then.
subplot(3,1,3)
if (size(list,1)> 0)   %plot if something there
    plot3(list(:,1),list(:,2),list(:,3),'k.'); box on;
end
xlabel('time (ms)'); grid on
ylabel('motounit')
set(gca,'ylim', [0 100]); 
set(gca,'xlim',[0 350])
view(0,90)
title('(c)')