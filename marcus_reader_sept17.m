%
% neurofield_read_sept17.m
% Marcus Wilson 4 September 2017
% Read in the basic neurofield output from the example of version 0.1.4
%
close all; %clear
%
% Read with the nf_read file
%
%%%%%  INPUT here %%%%%%%%%%%%%%%%%%

node_select=[1];   %a list of nodes that we will plot out. Must be taken from 1 to the number of nodes that the code outputs
show_plasticity=true;  %do we show plasticity output? frue of false

fid=fopen('example.output');   %Put in name of file we wish to open
%
%11 July 2018  Parameters for raster plots and spike generation
plot_phi_or_Q=false;  %if true work with phi spikes, if false use Q
N=150;  %number of neurons
do_spike_reconstruction=false;  %do we reconstruct spikes?  Make note of which population we need
sigma=0.15;  %the standard deviation
do_raster=true;  %will only work if we also do spike reconstruction
time_needed_array = [1000 1050];  %in ms - what points do we plot in the raster plot
nearest_pulse_time=1000;  %in ms  - when is the closest pulse to the period we want to plot?
resynch_time=0.2;   %resynch every burst of 5Hz rTMS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tline='some_string_to_start';
row=0;

%first, we'll find the number of nodes

while( ~( (tline(1)=='N')  & (tline(2)=='o')  & (tline(3)=='d') ) );
tline=fgetl(fid); 
if (length(tline)<1)
    tline='empty_line';
end
disp(tline)
row=row+1;
end
%Now have reached the row where the nodes are. 
disp('*****')
nnodes_sim=str2num(tline(7:end));

%8 September 2017.  Look for where we have our write out time and nodes for
%output
while(~( (tline(1) == 'O') & (tline(2) == 'u') & (tline(3) == 't') ) ) %look for the string 'Out' for 'Output'
tline=fgetl(fid);  %read next line
if (length(tline)<3)
    tline='empty_line';
end
disp(tline);
row=row+1;   %how many rows have we read in?
%Now read the write time (how many samples are written out in results per
%second)
end
disp(['***found the output line*** it is ' tline])
pos_interval=strfind(tline,'Inter')    %find 'Inter' (start of 'interval')
%The value will start in position 9 following this
write_time=str2num(tline(pos_interval+9:end));   %turn to number

%Find whether we do all nodes or not.  Identify how many nodes we have
%outputed
all_find=strfind(tline,'All')   
if (isempty(all_find))
    %if we couldn't find the word 'All' then we are only outputing some
    %nodes
    %
    %Find where 'Node:' ends and 'Start:' starts. That's where our output
    %nodes are listed. 
    node_find=strfind(tline,'Node:');
    start_find=strfind(tline,'Start:');
    nodes_outputed=str2num( tline(node_find+5:start_find-1) )  %look in this portion of tline
    nnodes=length(nodes_outputed);   %now identifying the number of nodes
else
    %we have found the word 'All'. So we do all nodes.Keep nnodes as we
    %found it earlier
    nnodes=nnodes_sim;  %nnodes for output is the same as the simulation nodes
end


%Find the output data
tline='empty_line';
while(~(tline(1) == '=') )  %while we don't have '===' at start, read in a new line
tline=fgetl(fid);   
if (length(tline)<3)
    tline='empty_line';
end
disp(tline);
row=row+1;   %how many rows have we read in?
end

%Now we read in three empty lines
tline1=fgetl(fid);
tline2=fgetl(fid);
tline3=fgetl(fid);
row=row+3;  %add in another three to the row count

%Now we start on the data
result=dlmread('example.output','',row,0);  %read in file, with blank line as delimiter and skipping first 'row' rows and zero columns

%Now should have each node available. Formatted tt, Q for all nodes, V for
%all nodes, phi for all nodes



tt=result(:,1);  %seconds
Q=result(:,2:nnodes+1);  %per second
V=result(:,nnodes+2:2*nnodes+1); %volts
phi=result(:,2*nnodes+2:3*nnodes+1); %per second. The propagator

%Check whether we have long enough data to have printed out plasticity data
%  6 November 2017
if (size(result,2)>5)
    %correct length - carry on
    nu=result(:,3*nnodes+2:4*nnodes+1);
    nutilde=result(:,4*nnodes+2:5*nnodes+1);
    Ca=result(:,5*nnodes+2:6*nnodes+1);
else
    disp('result matrix not long enough to contain plasticity data')
    show_plasticity=false;
end

%  Q_square=zeros(length(tt),sqrt(nnodes),sqrt(nnodes));
%  for t=1:length(tt)
%  Qshape=reshape(Q(t,:),sqrt(nnodes),sqrt(nnodes));   %puts in a 2x2 matrix
%  Q_square(t,:,:)=Qshape(:,:); 
%  end   %Now Q_square is a 3d array with time, x, y as indices for the firing rate;

figure(1);
subplot(3,2,1)
plot(tt,Q(:,node_select))
xlabel('time (s)')
ylabel('firing rate (/s)')
title('output for selected nodes')
set(gca,'xlim',[0 10])

subplot(3,2,2)
plot(tt,1000*V(:,node_select))
xlabel('time (s)')
ylabel('voltage (mV)')

figure(1);
subplot(3,2,3)
plot(tt,phi(:,node_select))
xlabel('time (s)')
ylabel('axonal flux (/s)')

if (show_plasticity)

figure(1);
subplot(3,2,4)
plot(tt,nu(:,node_select))
xlabel('time (s)')
ylabel('synaptic coupling (V s)')

figure(1);
subplot(3,2,5)
plot(tt,nutilde(:,node_select))
xlabel('time (s)')
ylabel('target synaptic coupling (V s)')

figure(1);
subplot(3,2,6)
plot(tt,Ca(:,node_select))
xlabel('time (s)')
ylabel('Calcium concentration (M)')

end

% 8 September 2017
% Power desnity function calculation on last half of data
figure(2)
fs=1/write_time;
pxx_total=0;
for i=1:length(node_select)
    j=node_select(i);  %select the node
    Qsamps=Q(length(Q)/2+1:end,j);   %take the Q series at the i-th node
    Qsamps=Qsamps-mean(Qsamps);   %make it average zero
    [pxx,f]=pwelch(Qsamps,[],[],[],fs);   %do a pwelch power sample
    %now do a plot 
    loglog(f,pxx); hold on
    pxx_total=pxx_total+pxx;   %add to runing total;
end
pxx_average=pxx_total/length(node_select);  
%loglog(f,pxx_average,'k','linewidth', 2); grid on   %average of these
xlabel('frequency (Hz)')
ylabel('power in Q (s^{-2}/Hz)')

%2 July 2018
%integrate for spike trains
summed_up=NaN*Q;   
sum_so_far=0;
for i=1:length(Q);
    sum_so_far = sum_so_far + Q(i);
    summed_up(i)=sum_so_far;
end
integrated_up=summed_up*(tt(2)-tt(1));

figure(3)
subplot(3,1,1)
plot(tt,integrated_up);
xlabel('time (s)')
ylabel('firings')

identify_when = floor(integrated_up(2:end)) - floor(integrated_up(1:end-1));
subplot(3,1,2)
plot(tt(1:end-1),identify_when)
sum(identify_when)

%poisson version
Q_delta_t = Q*(tt(2)-tt(1));   %this is probability of spike for a small time
%random_nums = unifrnd(0,1,length(Q),1);
random_nums=0.5*ones(length(Q),1);
%if these are less than Q_delta_t we get a spike
spike_here = (random_nums < Q_delta_t);
subplot(3,1,3)
plot(tt,spike_here)
sum(spike_here)

deltat=tt(2)-tt(1);
if (length(Q)*deltat > 2.01)
     disp(['Maximum Q in first two second is ' num2str(max(Q(1:round(2/deltat)))) ' s-1'])

one_sec_point = round(1/deltat);
Q_baseline=Q(one_sec_point) %at one second
Q_half = 0.5*(max(Q(1:2*one_sec_point))+Q_baseline)
%find when we get to half amplitude
i=one_sec_point;
while ( (Q(i) < Q_half)  & (i < length(Q)) )
    i=i+1;
end
time_delay = (i-one_sec_point)*(tt(2)-tt(1));
disp(['time to half amplitude is ' num2str(1000*time_delay) ' ms'])
end
  
figure(4)
%17 August 2018
%p6242 Try making an actual MEP out of all of this
lambda=0.0230; %s-1
mepsize=exp(lambda*Q)-1;   %proportional MEP
plot(tt,mepsize(:,node_select))
xlabel('time (s)')
ylabel('MEP size (arbitrary units)')
set(gca,'xlim',[0 10])
title('A MEP done in Moezzi style. Only relevant if we had population 7 outputed')



if (do_spike_reconstruction)

%11 July 2018
%Construct synthetic spike trains for lots of cells



mean_value=1.0;
mult_facs=normrnd(mean_value,sigma,1,N);  %N random numbers
if (plot_phi_or_Q)
     phi_multiple = phi*mult_facs;  %150 columns
else
     phi_multiple=Q*mult_facs;   %its actually Q but we'll call it phi here
end
raster_record=zeros(length(phi),N);  %record 1 and 0 for plotting here
for i=1:N %go over all columns
    int_value=0;
    Nspike=0;  %how many spikes so far. 
    for j=1:length(phi);   
        %12 July 2018. Need to resynch every so many seconds
        if ( mod(j,resynch_time/deltat)==0)  %if j hits a multiple of resynch_time/deltat - the resynch number of steps
           %set the integral back to zero
            int_value=0;
        end
        int_new=int_value+phi_multiple(j,i)*deltat;  %look at jth row (time) and ith column (which cell)
        %is this new value going through an integer?
        if ((floor(int_new) - floor(int_value)) == 1)
            %we have gone through an integer
            Nspike=Nspike+1;  %add one to the total spike count
            spike_time(Nspike,i) = j*deltat;   %enter the time in seconds.  This array will keep changing in size. Not the best coding.
            raster_record(j,i)=1;
        end
        int_value=int_new;  %update the integral
    end  %next time point
end  %next cell


if (do_raster)
%Do raster plot
figure(4)
points_needed = (time_needed_array(2)-time_needed_array(1))/(1000*deltat); %number of points needed to plot
start_time=time_needed_array(1)/1000;
start_point=start_time/deltat;
X=[1:points_needed]'*deltat + start_time - nearest_pulse_time/1000;  %the time values, in seconds, calibrated from the TMS pulse
Y=[1:N]';
surf(1000*X,Y,raster_record(start_point+1:points_needed+start_point,:)','linestyle','none');
xlabel('time (ms)')
if (plot_phi_or_Q)
   ylabel('axon number')
   title('raster plot of spikes arriving at motoneurons')
else
    ylabel('layer 5 cell number')
    title('raster plot of spikes generated by layer 5 neurons')
end
view(0,90)
end

%save our output
save 'spike_time_file.mat' spike_time

save 'key_output_parameters.mat' tt Q V phi

end    
    
%Last of all, save our tt Q V phi anyway
save 'key_output_parameters.mat' tt Q V phi


figure(6)
%plot some output pics for the paper
subplot(2,1,1)
plot([1:4000]*1e-1-20,phi(40001:44000),'k'); grid on
xlabel('time (ms)')
ylabel('layer 5 pulse rate (/s)')
title('(a)')
set(gca,'xlim',[0 380])

