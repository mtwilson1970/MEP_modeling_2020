%This function generates the times at which a firing occurs by simple time
%integration. 
function spike_times = get_firing_times(tt, firing_rate)

%firing rate has rows corresponding to time points and columns
%corresponding to motounits

%Integrate firing rate
deltat = tt(2)-tt(1);  %this is the time step
integrated = zeros(1,size(firing_rate,2)); 
%       %Marcus Wilson 18 May 2020. Put in a random phase. p6943
%       integrated=rand(1,size(firing_rate,2));   %this is random numbers 0 to 1;
%       %
int_matrix = zeros(length(tt),size(firing_rate,2)); 
for time_step = 1:length(tt)   %go over all time steps
   integrated = integrated + deltat*firing_rate(time_step,:);   %integrate in time up to the required time. Integrated is the result up to this time
   %set phase to zero if firing_rate is actually 'zero' (less than 0.1)
   %hold to the integer it has got to. 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %18 May 2020. Snip out this code to give random phases.
  for j=1:size(firing_rate,2)
      if (firing_rate(time_step,j) < 0.1)   %if we have zero firing we go back to last integer
          integrated(j)=floor(integrated(j));
      end
  end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %integrated = integrated + (firing_rate(time_step,:) < 0.1).*(integrated - floor(integrated))
   
   int_matrix(time_step,:) = integrated(1,:);   
end

%Now work out when we go over an integer
num_firings_to_date = [floor(int_matrix) ;  NaN*zeros(1,size(firing_rate,2))];   %just take the integer value and a row of NaNs at the end;
num_firings_to_date_displace = [NaN*zeros(1,size(firing_rate,2)) ;  floor(int_matrix)]; %put in a row at the front

firing_has_happened = num_firings_to_date - num_firings_to_date_displace; 
%make first a zero not a NaN
firing_has_happened(1,:)=zeros(1,size(firing_rate,2));  

%knock off the last row
spike_times=firing_has_happened(1:end-1,:);   %a 1 when a firing happens, 0 otherwise

end