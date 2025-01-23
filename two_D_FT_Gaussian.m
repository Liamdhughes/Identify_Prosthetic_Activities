function y = two_D_FT_Gaussian(track,MF,TR,SR,FWHM);


window = exp(-(1/((FWHM*SR)/(2*sqrt(log(2))))^2)*((1:length(track)) - 0.5*length(track)).^2);

temp = Spectrogram_revised_by_Martin_Gaussian(track.*window,MF,SR,TR,TR/2);
tdft = abs(fft((temp')))';

y = tdft(:,1:round(size(tdft,2)/2));
