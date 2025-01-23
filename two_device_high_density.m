% Copyright, M.Bencsik, M.Bisele L.D.Hughes, 2025
clear
close all

Disk_Letter_Name = 'D';


path(path,[Disk_Letter_Name,':\Liam Passport\Study 3\data\Useful_Matlab_Functions']);

% This code builds a training data base of 500 2DFT's per sensor, per
% individual, so as to train the algorithm sensor by sensor, for one
% specific individual:
Individual_No_array = ['01'];

time_increment = 0.1;

temporal_resolution = 0.07;
feature_length = 2.8;
multiplication_factor = 2;
    
  for S_R = 40
    
    for index = 1:size(Individual_No_array,1)
        
        Individual_No = Individual_No_array(index,:);
                
        % Upload the corresponding data set:
        load([Disk_Letter_Name,':\Liam Passport\Study 3\data\control\P0',Individual_No,'\corrected_data_at_',num2str(S_R),'_Hz.mat'])
        % Upload the starting times for the TDB:
        starting_times = csvread([Disk_Letter_Name,':\Liam Passport\Study 3\data\control\P0',Individual_No,'\TDB_starting_times.csv']);
        
        time_axis = (1/S_R)*(0:(length(A_Sacrum_interpolated)-1));
        
        %** process the Sacrum **
        
        counter = 1;
        for uu = 1:size(starting_times,1)
            for vv = 1:size(starting_times,2)
                %find the index of the data according to the starting time:
                [a b] = min(abs(time_axis - starting_times(uu,vv)));
                [a bmax] = min(abs(time_axis - (starting_times(uu,vv) + 4 - feature_length)));
                while (b < bmax)
                    % compute the 2DFT of the relevant section of accelerometer:
                    temp = two_D_FT_Gaussian(A_Sacrum_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
                    TDB_1(:,counter) = temp(:);
                    counter = counter + 1;
                    %increase time by 'time_increment'
                    b = round(b + time_increment*S_R);
                end
            end
            boundary(uu) = counter -1;
        end
        
        
        
        %** process the Thigh L **
        
        counter = 1;
        for uu = 1:size(starting_times,1)
            for vv = 1:size(starting_times,2)
                %find the index of the data according to the starting time:
                [a b] = min(abs(time_axis - starting_times(uu,vv)));
                [a bmax] = min(abs(time_axis - (starting_times(uu,vv) + 4 - feature_length)));
                while (b < bmax)
                    % compute the 2DFT of the relevant section of accelerometer:
                    temp = two_D_FT_Gaussian(A_Thigh_L_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
                    TDB_2(:,counter) = temp(:);
                    counter = counter + 1;
                    %increase time by 'time_increment'
                    b = round(b + time_increment*S_R);
                end
            end
        end
        
        
        %** process the Thigh R **
        
        counter = 1;
        for uu = 1:size(starting_times,1)
            for vv = 1:size(starting_times,2)
                %find the index of the data according to the starting time:
                [a b] = min(abs(time_axis - starting_times(uu,vv)));
                [a bmax] = min(abs(time_axis - (starting_times(uu,vv) + 4 - feature_length)));
                while (b < bmax)
                    % compute the 2DFT of the relevant section of accelerometer:
                    temp = two_D_FT_Gaussian(A_Thigh_R_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
                    TDB_3(:,counter) = temp(:);
                    counter = counter + 1;
                    %increase time by 'time_increment'
                    b = round(b + time_increment*S_R);
                end
            end
        end
        
        
        %** process the Shank L **
        
        
        counter = 1;
        for uu = 1:size(starting_times,1)
            for vv = 1:size(starting_times,2)
                %find the index of the data according to the starting time:
                [a b] = min(abs(time_axis - starting_times(uu,vv)));
                [a bmax] = min(abs(time_axis - (starting_times(uu,vv) + 4 - feature_length)));
                while (b < bmax)
                    % compute the 2DFT of the relevant section of accelerometer:
                    temp = two_D_FT_Gaussian(A_Shank_L_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
                    TDB_4(:,counter) = temp(:);
                    counter = counter + 1;
                    %increase time by 'time_increment'
                    b = round(b + time_increment*S_R);
                end
            end
        end
        
        %** process the Shank R **
        
        counter = 1;
        for uu = 1:size(starting_times,1)
            for vv = 1:size(starting_times,2)
                %find the index of the data according to the starting time:
                [a b] = min(abs(time_axis - starting_times(uu,vv)));
                [a bmax] = min(abs(time_axis - (starting_times(uu,vv) + 4 - feature_length)));
                while (b < bmax)
                    % compute the 2DFT of the relevant section of accelerometer:
                    temp = two_D_FT_Gaussian(A_Shank_R_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
                    TDB_5(:,counter) = temp(:);
                    counter = counter + 1;
                    %increase time by 'time_increment'
                    b = round(b + time_increment*S_R);
                end
            end
        end
        
     for fig_No = 1:5
                subplot(2,3,fig_No)
                eval(['imagesc(log10(TDB_',num2str(fig_No),'))'])
                title(['Sampling Rate = ',num2str(S_R)])
                colormap(jet(256))
     end
            pause(1)
        
        save([Disk_Letter_Name,':\Liam Passport\Study 3\data\control\processing_code\training_data_bases\high_density\high_density_TDB_',Individual_No,'optimised',num2str(S_R),'.mat'],'TDB_1','TDB_2','TDB_3','TDB_4','TDB_5','boundary')
        clear TDB_1 TDB_2 TDB_3 TDB_4 TDB_5 boundary
    end
    
end
