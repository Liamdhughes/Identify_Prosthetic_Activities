% Copyright, M.Bencsik, M.Bisele L.D.Hughes, 2025

participant_No = '1';

path_name = ['D:\Liam Passport v2\Study 3\data\amputee\P001\EXT\'];


files_name = ['P00',participant_No,'_EXT_one.xlsx'];

S_R = 100;

data = xlsread([path_name,files_name]);




ruler_drop = 4*60 + 0;

% A_Sacrum = sqrt((data(:,2)).^2 + (data(:,3)).^2 + (data(:,4)).^2);
A_Pros = sqrt((data(:,2)).^2 + (data(:,3)).^2 + (data(:,4)).^2);

% A_Thigh_R = sqrt((data(:,10)).^2 + (data(:,11)).^2 + (data(:,12)).^2);
% 
% A_Shank_L = sqrt((data(:,14)).^2 + (data(:,15)).^2 + (data(:,16)).^2);
% %
% A_Shank_R = sqrt((data(:,18)).^2 + (data(:,19)).^2 + (data(:,20)).^2);
% %
% 
% subplot(2,1,1)
% time_A_Sacrum = (60/60)*(1/S_R)*(0:(length(A_Sacrum)-1));
% plot(time_A_Sacrum,A_Sacrum)
% xlim([ruler_drop-8 ruler_drop-1])


subplot(2,1,2)
time_A_Pros = (60/60)*(1/S_R)*(0:(length(A_Pros)-1));
plot(time_A_Pros,A_Pros)
% xlim([ruler_drop-8 ruler_drop-1])

% subplot(5,1,3)
% time_A_Thigh_R = (60/60)*(1/S_R)*(0:(length(A_Thigh_R)-1));
% plot(time_A_Thigh_R,A_Thigh_R)
% xlim([ruler_drop-8 ruler_drop-1])
% 
% subplot(5,1,4)
% 
% time_A_Shank_L = (60/60)*(1/S_R)*(0:(length(A_Shank_L)-1));
% plot(time_A_Shank_L,A_Shank_L)
% xlim([ruler_drop-8 ruler_drop-1])
% 
% subplot(5,1,5)
% time_A_Shank_R = (60/60)*(1/S_R)*(0:(length(A_Shank_R)-1));
% plot(time_A_Shank_R,A_Shank_R)
% xlim([ruler_drop-8 ruler_drop-1])

save('raw_one_sensor.mat','A_Pros','time_A_Pros')
