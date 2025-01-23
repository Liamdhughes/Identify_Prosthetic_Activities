clear
close all

Disk_Letter_Name = 'D';

path(path,[Disk_Letter_Name,':\Liam Passport v2\Study 3\data\Useful_Matlab_Functions']);

% This code builds a training data base of 30 2DFT's per sensor, per
% individual, so as to train the algorithm sensor by sensor, for one
% specific individual:

 Individual_No_array = ['01';'02';'03';'04';'05';'06';'07';'08';'09';'10'];
%Individual_No_array = ['10'];

temporal_resolution = 0.02; 
feature_length = 3;
multiplication_factor = 2;
S_R = 100
    
    for index = 1:size(Individual_No_array,1)
        Individual_No = Individual_No_array(index,:);       
        
        % Upload the corresponding data set:
        load([Disk_Letter_Name,':\Liam Passport v2\Study 3\data\amputee\P0',Individual_No,'\EXT\corrected_data_two.mat'])
        % Upload the starting times for the TDB:
        starting_times = csvread([Disk_Letter_Name,':\Liam Passport v2\Study 3\data\amputee\P0',Individual_No,'\EXT\TDB_starting_times2.csv']);
        
        time_axis = (1/S_R)*(0:(length(A_Sacrum_interpolated)-1));
        
        %** process the Sacrum **
        
        counter = 1;
        for uu = 1:size(starting_times,1)
            for vv = 1:size(starting_times,2)
                %find the index of the data according to the starting time:
                [a b] = min(abs(time_axis - starting_times(uu,vv)));
                % compute the 2DFT of the relevant section of accelerometer:
                temp = two_D_FT_Gaussian(A_Sacrum_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
                dim1=size(temp,1);
                dim2=size(temp,2);
                TDB_1(:,counter) = temp(:);
                counter = counter + 1;
            end
        end
        
        
        %** process the Thigh L **
        
        
        counter = 1;
        for uu = 1:size(starting_times,1)
            for vv = 1:size(starting_times,2)
                %find the index of the data according to the starting time:
                [a b] = min(abs(time_axis - starting_times(uu,vv)));
                % compute the 2DFT of the relevant section of accelerometer:
                temp = two_D_FT_Gaussian(A_Pros_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
                TDB_2(:,counter) = temp(:);
                counter = counter + 1;
            end
        end
        
        %** process the Thigh R **
        
        % counter = 1;
        % for uu = 1:size(starting_times,1)
        %     for vv = 1:size(starting_times,2)
        %         %find the index of the data according to the starting time:
        %         [a b] = min(abs(time_axis - starting_times(uu,vv)));
        %         % compute the 2DFT of the relevant section of accelerometer:
        %         temp = two_D_FT_Gaussian(A_Thigh_R_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
        %         TDB_3(:,counter) = temp(:);
        %         counter = counter + 1;
        %     end
        % end
        % 
        % %** process the Shank L **
        % counter = 1;
        % for uu = 1:size(starting_times,1)
        %     for vv = 1:size(starting_times,2)
        %         %find the index of the data according to the starting time:
        %         [a b] = min(abs(time_axis - starting_times(uu,vv)));
        %         % compute the 2DFT of the relevant section of accelerometer:
        %         temp = two_D_FT_Gaussian(A_Shank_L_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
        %         TDB_4(:,counter) = temp(:);
        %         counter = counter + 1;
        %     end
        % end
        % 
        % %** process the Shank R **
        % counter = 1;
        % for uu = 1:size(starting_times,1)
        %     for vv = 1:size(starting_times,2)
        %         %find the index of the data according to the starting time:
        %         [a b] = min(abs(time_axis - starting_times(uu,vv)));
        %         % compute the 2DFT of the relevant section of accelerometer:
        %         temp = two_D_FT_Gaussian(A_Shank_R_interpolated(b:round(b+feature_length*S_R)),multiplication_factor,temporal_resolution,S_R,0.5*feature_length);
        %         TDB_5(:,counter) = temp(:);
        %         counter = counter + 1;
        %     end
        % end
        % 
        
            % monitor the outcome of the training data bases:
            for fig_No = 1:2
                subplot(2,1,fig_No)
                eval(['imagesc(log10(TDB_',num2str(fig_No),'))'])
                title(['Sampling Rate = ',num2str(S_R)])
                colormap(jet(256))
                print(['D:\Liam Passport v2\Study 3\data\amputee\processing_code_ext\low density\Ind_',num2str(Individual_No),'.tif'],'-dtiff','-r300')
            end
        
        save([Disk_Letter_Name,':\Liam Passport v2\Study 3\data\amputee\processing_code_ext\training_data_bases\TDB_',Individual_No,'optimised',num2str(S_R),'.mat'],'TDB_1','TDB_2','dim1','dim2')
        clear TDB_1 TDB_2 
        
        saveas(gcf,'Barchart.png')
    end
