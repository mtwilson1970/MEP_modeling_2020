% 12 October 2020. Making code available

%
% Marcus Wilson 25 June 2019
% 
% This code reads in a templated configuration file, modifies it, and
% writes out a conf file. Then it runs NFTsim on it. 

close all; clear

%template_to_read_from='configs/meps/marcus_template_june2019.conf';  %the file for amp curves
template_to_read_from='marcus_template_june2019.conf';  %p7208. Paired pulse curve


do_we_save=false;  %do we save this run?
if (~do_we_save)
    disp('***WARNING*** The output will not be saved')
end

output_file_base = 'MEP_CSP_against_amplitude_';
extra_identifying_text = 'standard_retest';   %extra text to identify which case we have.
output_filename = [output_file_base extra_identifying_text '.mat'];


ISIcurve=true; ampcurve=false;  gababcurve=false; MVCcurve=false; %Need to set other paras e.g. Burst number correctly in config template file

if (~(ISIcurve + ampcurve + gababcurve  +MVCcurve) ==1)
    disp('One and only one of the options must be true')
    stop
end

if (ISIcurve)
my_string = 'Burst Frequency:';
length_to_cover = 5;  %how many character elements need covering with the writing

my_values = [2000 1000 500 333 250 200 143 100 67 50 33 20 10 5 3 2];  %time intervals for paper
%my_values=[1000 500 333 250 200 167 143 125 111 100 91 83 77 71 67 63 59 56 53 50];  %integer delta t steps from 1 to 20 ms
my_values=[500];

end

if(ampcurve)
    my_string = 'Amplitude:';
    %my_values=[450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 850 900 950 1000 1100 1200 1300 1400 1500 1600];
       my_values=[780]
    length_to_cover=5;
    
    my_string_2 = 'Coupling 1:  Map  - nu:';
    my_values_2 = 1.0*1.92e-4*(exp((500-my_values)/100) + 1).^(-1);  %standard 1.0*  500-  /100
    
    length_to_cover_2=10;
    
    my_string_3 = 'Coupling 3:  Map  - nu:';   %external to inhib A
    my_values_3 = 1.0*(-0.72e-4)*ones(1,length(my_values));
    length_to_cover_3=11;
    
    my_string_4 = 'Coupling 4:  Map  - nu:';   %external to inhib B
    my_values_4 = my_values_3;  
    length_to_cover_4=11;
    
    my_string_5 = 'Coupling 16:  Map  - nu:';  %external to copy of e.
    my_values_5 = my_values_2;
    length_to_cover_5=10;
    
    my_string_6 = 'Coupling 21:  Map  - nu:';   %external to layer 5.
    my_values_6 = 1.0*0.2/1.92*my_values_2;
    length_to_cover_6=10;
    
    
    
end

if (gababcurve)
    
    
    
    %select alpha and beta values
    
    
    %Find the list of strings to locate
    my_common_string = 'Dendrite: '
    dendrites_to_change = [3 8 13 18 23];
   
    for k=1:length(dendrites_to_change)
        stringy=[my_common_string num2str(dendrites_to_change(k))];
        if (dendrites_to_change(k) < 10)
            %add a space at the end
            stringy=[stringy ' ']; 
        end
        dend_strings(k,:) = stringy(:);   %Now load up the dend_strings character matrix
    end
  
end

if (MVCcurve)
    my_string = 'Mean:'
    length_to_cover=2;   %how many spots to cover
    my_values=[0:2:12];
end

%Read in a template
%fid=fopen('neurofield_MEPmeasure.conf','r');
fid=fopen(template_to_read_from,'r');    %Reads in the neurofield template.
configfile=fread(fid,'char=>char')' ;  %read in all the data



 if ( ampcurve||ISIcurve||MVCcurve )   %if we do one of these we'll work on the first my_string bit 
    
    %Now we need to find what we want to write
 i=1; 
    comparing=[0];  %null to start
    while ( ~(sum(comparing) == length(comparing)) ) ;   %if they are all ones then we have a matc     
        i=i+1; 
       if (i > length(configfile) ) 
           disp('End of file reached')
           stop
       end
      comparing=(configfile(i:i+length(my_string)-1) == (my_string(:))');  %This does a letter-by-letter comparison with the config file and my selected string.

    end
    disp(['match at i=' num2str(i)]);
    
%Now need to write out the new value in the config file
%Write out next five points

    where_we_write=i+length(my_string);   %j becomes the index of the value we want to write
    
  
 end
    
 
    if (ampcurve);  %Now repeat for the other bits.
      

 %Now we need to find what we want.
 i=1; 
    comparing=[0];  %null to start
    while ( ~(sum(comparing) == length(comparing)) ) ;   %if they are all ones then we have a matc     
        i=i+1; 
       if (i > length(configfile) ) 
           disp('End of file reached')
           stop
       end
      comparing=(configfile(i:i+length(my_string_2)-1) == (my_string_2(:))');  %This does a letter-by-letter comparison with the config file and my selected string.

    end
    disp(['match at i=' num2str(i)]);
    
%Now need to write out the new value in the config file
%Write out next five points

    where_we_write_2=i+length(my_string_2);   %j becomes the index of the value we want to write
    
    
%Now we need to find what we want.
 i=1; 
    comparing=[0];  %null to start
    while ( ~(sum(comparing) == length(comparing)) ) ;   %if they are all ones then we have a matc     
        i=i+1; 
       if (i > length(configfile) ) 
           disp('End of file reached')
           stop
       end
      comparing=(configfile(i:i+length(my_string_3)-1) == (my_string_3(:))');  %This does a letter-by-letter comparison with the config file and my selected string.

    end
    disp(['match at i=' num2str(i)]);
    
%Now need to write out the new value in the config file
%Write out next five points

    where_we_write_3=i+length(my_string_3);   %j becomes the index of the value we want to write

 
    
    
    %Now we need to find what we want.
 i=1; 
    comparing=[0];  %null to start
    while ( ~(sum(comparing) == length(comparing)) ) ;   %if they are all ones then we have a matc     
        i=i+1; 
       if (i > length(configfile) ) 
           disp('End of file reached')
           stop
       end
      comparing=(configfile(i:i+length(my_string_4)-1) == (my_string_4(:))');  %This does a letter-by-letter comparison with the config file and my selected string.

    end
    disp(['match at i=' num2str(i)]);
    
%Now need to write out the new value in the config file
%Write out next five points

    where_we_write_4=i+length(my_string_4);   %j becomes the index of the value we want to write   
    
 
     %Now we need to find what we want.
 i=1; 
    comparing=[0];  %null to start
    while ( ~(sum(comparing) == length(comparing)) ) ;   %if they are all ones then we have a matc     
        i=i+1; 
       if (i > length(configfile) ) 
           disp('End of file reached')
           stop
       end
      comparing=(configfile(i:i+length(my_string_5)-1) == (my_string_5(:))');  %This does a letter-by-letter comparison with the config file and my selected string.

    end
    disp(['match at i=' num2str(i)]);
    
%Now need to write out the new value in the config file
%Write out next five points

    where_we_write_5=i+length(my_string_5);   %j becomes the index of the value we want to write     
    
    
    
    %Now we need to find what we want.
 i=1; 
    comparing=[0];  %null to start
    while ( ~(sum(comparing) == length(comparing)) ) ;   %if they are all ones then we have a matc     
        i=i+1; 
       if (i > length(configfile) ) 
           disp('End of file reached')
           stop
       end
      comparing=(configfile(i:i+length(my_string_6)-1) == (my_string_6(:))');  %This does a letter-by-letter comparison with the config file and my selected string.

    end
    disp(['match at i=' num2str(i)]);
    
%Now need to write out the new value in the config file
%Write out next five points

    where_we_write_6=i+length(my_string_6);   %j becomes the index of the value we want to write      
    
    
    
    end   %if its an amp curve.


  
  %Loop around the variable so we can change it

if (ampcurve||ISIcurve||MVCcurve)

  %loop over these values
for count=1:length(my_values)
  
    
    
    
 disp(['Working on run for ' my_string ' ' num2str(my_values(count)) ])


chars_to_write = num2str(my_values(count));   %These need writing, but we want to do the right number of spaces before.  
spaces_needed = length_to_cover - length(chars_to_write);  %find out how many spaces are needed.

for k=1:spaces_needed;
    chars_to_write = [' ' chars_to_write];  %add a space in front.
end

%Now write this on the file
configfile(where_we_write:where_we_write+length_to_cover-1) = chars_to_write';


if (ampcurve)

    
    %repeat all this several times
    
    %loop over these values


    disp(['Working on run for ' my_string ' ' num2str(my_values_2(count)) ])


    chars_to_write = num2str(my_values_2(count));   %These need writing, but we want to do the right number of spaces before.  
    spaces_needed = length_to_cover_2 - length(chars_to_write);  %find out how many spaces are needed.

    for k=1:spaces_needed;
       chars_to_write = [' ' chars_to_write];  %add a space in front.
    end

    %Now write this on the file
    configfile(where_we_write_2:where_we_write_2+length_to_cover_2-1) = chars_to_write';
    
    


    disp(['Working on run for ' my_string ' ' num2str(my_values_3(count)) ])


    chars_to_write = num2str(my_values_3(count));   %These need writing, but we want to do the right number of spaces before.  
    spaces_needed = length_to_cover_3 - length(chars_to_write);  %find out how many spaces are needed.

    for k=1:spaces_needed;
        chars_to_write = [' ' chars_to_write];  %add a space in front.
    end

    %Now write this on the file
    configfile(where_we_write_3:where_we_write_3+length_to_cover_3-1) = chars_to_write';
    




    disp(['Working on run for ' my_string ' ' num2str(my_values_4(count)) ])


    chars_to_write = num2str(my_values_4(count));   %These need writing, but we want to do the right number of spaces before.  
    spaces_needed = length_to_cover_4 - length(chars_to_write);  %find out how many spaces are needed.

    for k=1:spaces_needed;
        chars_to_write = [' ' chars_to_write];  %add a space in front.
    end

    %Now write this on the file
    configfile(where_we_write_4:where_we_write_4+length_to_cover_4-1) = chars_to_write';
    

    
    disp(['Working on run for ' my_string ' ' num2str(my_values_5(count)) ])


    chars_to_write = num2str(my_values_5(count));   %These need writing, but we want to do the right number of spaces before.  
    spaces_needed = length_to_cover_5 - length(chars_to_write);  %find out how many spaces are needed.

    for k=1:spaces_needed;
        chars_to_write = [' ' chars_to_write];  %add a space in front.
    end

    %Now write this on the file
    configfile(where_we_write_5:where_we_write_5+length_to_cover_5-1) = chars_to_write';
    
    
    disp(['Working on run for ' my_string ' ' num2str(my_values_6(count)) ])


    chars_to_write = num2str(my_values_6(count));   %These need writing, but we want to do the right number of spaces before.  
    spaces_needed = length_to_cover_6 - length(chars_to_write);  %find out how many spaces are needed.

    for k=1:spaces_needed;
        chars_to_write = [' ' chars_to_write];  %add a space in front.
    end

    %Now write this on the file
    configfile(where_we_write_6:where_we_write_6+length_to_cover_6-1) = chars_to_write';
    

end


%Now we should write the file and launch NFTsim
fclose(fid);
status=system(sprintf('rm written_file.conf'))
fid = fopen('written_file.conf','w');
                        
%Now save a new configfile
fwrite(fid,configfile','char');
fclose(fid);

status=system(sprintf('rm example.output'))
status=system(sprintf('../bin/neurofield -i written_file.conf -o example.output'))
output_args =status;
pause(30);   %make sure it runs

marcus_reader_sept17;   %read in, evaluate meps
motomodel;

%read in ave_mep

load 'mep_file.mat'
maxmep(count)=max(ave_mep)

if (ISIcurve)
%make sure we delay to the test pulse
test_index = (1000/my_values(count))*10;  %this is the index time of test pulse
%take 90% of this time. 
maxmep(count) = max(ave_mep(floor(0.9*test_index):end))
end


%Need to find the silent period length. Wait till the MEP is finished, then
%find time when activity comes back.

isactivity = (abs(ave_mep)>0.015);   %are we above 10 mV?
  %find where, after a min of 50ms, isactivity comes in:  Needs to be
  %silent for 25ms at least to call it a CSP. 
  
  %first find where the mep is
  maxfound=0; jmax=0;
  for j=1:10000;
      valfound=ave_mep(j);
          if (valfound > maxfound) 
              maxfound=valfound;
              jmax=j;  %index of the maximum bit
          end
  end
              
  
  j=300;   %start counting at 30msFind where we go quiet
  flaggedit=0;   %haven't found it yet. 
  
  while((flaggedit==0) & (j<10000))
      
      %Go over all j points
  
      while ( ~(isactivity(j) == 0) & (j<10000)); %while there's some activity, go to the next bit. 
          j=j+1;    
      end
  
     %Now we have no activity. Count the quiet bits.
     p=0;
     while  ((isactivity(j+p) == 0) & (j+p < 10000)); %while there's no activity, count up p.
         p=p+1;
     end
  
     %Now we have activity back. How big is p?
     if (p>250);
         %We have a sufficent gap
         CSPtime(count) = (j+p-jmax)/10;   %CSP in ms
         flaggedit=1;  %we have identified a silent period.
     else
         %We don't have a long enough gap. Keep going with the counting. 
         j=j+p;
     end
  end
  
    if (flaggedit==0)
      CSPtime(count)=NaN
  end
  
  
  
  end
  
  


end






% if (gababcurve)
%     
%    % We work on the gabab curve here
%    
%    %Write the config file. First find the appropriate dendrite term
%    
%    for k=1:length(dendrites_to_change)
%        
%        %search for this text
%         i=1; 
%     comparing=[0];  %null to start
%     while ( ~(sum(comparing) == length(comparing)) ) ;   %if they are all ones then we have a matc     
%         i=i+1; 
%        if (i > length(configfile) ) 
%            disp('End of file reached')
%            stop
%        end
%       comparing=(configfile(i:i+length(dend_strigs)-1) == (dend_strings(k,:))');  %This does a letter-by-letter comparison with the config file and my selected string.
% 
%     end
%     disp(['match at i=' num2str(i)]);
%     
% %Now need to write out the new values in the config file
% 
% %work out the spaces needed
% length_needed = 5;
%    text_to_write = ['  alpha:'   'beta:'   ]
%    
%    %Still working on this. Not sure we need it for paper....
%    
%   
%    end
% 
% 
% 
% end

%Now we have a variable maxmep for each of the options done
%close all

figure(21)
if (ISIcurve)
ISI=1000*my_values.^(-1);   
    plot(ISI,maxmep,'kx-'); grid on; xlabel('ISI (ms)'); ylabel('MEP (mV)')
end

if (ampcurve)  
    plot(my_values,maxmep,'kx'); grid on; xlabel('amplitude (/s)'); ylabel('MEP (mV)')   
end

if (MVCcurve)
    plot(my_values*2,maxmep,'kx'); grid on; xlabel('%MVC'); ylabel('MEP (mV)')
end


figure(22)
if (ISIcurve)
    plot(ISI,CSPtime,'kx-'); grid on; xlabel('ISI (ms)'); ylabel('CSP (ms)')
end

if (ampcurve)  
    plot(my_values,CSPtime,'kx-'); grid on; xlabel('amplitude (/s)'); ylabel('CSP (ms)')   
end

if (MVCcurve)
    plot(my_values*2,CSPtime,'kx-'); grid on; xlabel('%MVC'); ylabel('CSP (ms)')
end

if (do_we_save)
save(output_filename, 'my_values', 'maxmep', 'CSPtime');
end



