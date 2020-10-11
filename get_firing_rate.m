% This function works out what firing rates each of the motounits will
% have, given an input
function Q = get_firing_rate(phi, thresholds, gradient_for_each, intercept_for_each)
%first, work out how many units are firing
%Compare phi to threshold

Q = zeros(length(phi),length(thresholds)); 
for time_point = 1:length(phi)  %go through the length of phi in time

    which_are_firing = (phi(time_point)>thresholds);  %this gives us a list of 1 s and then 0 s

    rate_of_firing = phi(time_point)*gradient_for_each + intercept_for_each; %gives the firing rate for each one (this one can go below 8 Hz)
    
    Q(time_point,:) = which_are_firing.*rate_of_firing; %give the sequences for each
end
%Now have Q for each time point for each neuron
end