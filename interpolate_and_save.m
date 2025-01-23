% Copyright, M.Bencsik, M.Bisele L.D.Hughes, 2025

clear
close all

load cropped_data_two
load start_timings_two

S_R = 100;

time_duration = (105*60+45)*(101.299/101.75);

time_axis = 0:(1/100):time_duration;

% Undertake the interpolation, assuming a perfect sampling rate of 100 Hz:
A_Sacrum_interpolated = interp1(0:((time_duration)/(length(A_Sacrum_cropped) - 1)):time_duration,A_Sacrum_cropped,time_axis,'pchip');
A_Pros_interpolated = interp1(0:((time_duration)/(length(A_Pros_cropped) - 1)):time_duration,A_Pros_cropped,time_axis,'pchip');

figure
subplot(2,1,1)
plot(time_axis,A_Sacrum_interpolated)
title(['new sampling rate estimate: ',num2str(100*size(A_Sacrum_interpolated,2)./size(A_Sacrum_cropped,1))])
axis tight
subplot(2,1,2)
plot(time_axis,A_Pros_interpolated)
title(['new sampling rate estimate: ',num2str(100*size(A_Pros_interpolated,2)./size(A_Pros_cropped,1))])
axis tight

save('corrected_data_two.mat','A_Sacrum_interpolated','A_Pros_interpolated')
