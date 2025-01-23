function y = Spectrogram_revised_by_Martin_Gaussian(sound,Multiplicateur,S_R,s_l,FWHM)

% this function expects the raw sound, the additional density of spectra,
% the sampling rate and the sample length (in seconds) for the FFT.

FFT_length = round(s_l*S_R);

total_seconds = length(sound)/S_R;

% Calculate the number of individual spectra that can be obtained from
% within the file:
total_No_spectra = floor(total_seconds/(FFT_length/S_R));

% Initialise the memory to accomodate the spectrogram with overlap:
Sp_with_overlap = zeros(floor(FFT_length/2),(total_No_spectra+(total_No_spectra-1)*(Multiplicateur-1)));

% Compute the first and last spectra:
C = sound(1:(total_No_spectra*FFT_length));
window = exp(-(1/((FWHM*S_R)/(2*sqrt(log(2))))^2)*((1:FFT_length) - 0.5*FFT_length).^2)'*ones(1,total_No_spectra);
B = abs(fft(window.*reshape(C,FFT_length,total_No_spectra)));
% Store one half of the points (because the power spectrum is
% symmetrical)
temp = B(1:floor(size(B,1)/2),:);

Sp_with_overlap(:,1:Multiplicateur:(total_No_spectra+(total_No_spectra-1)*(Multiplicateur-1))) = temp;

LL = round(FFT_length/Multiplicateur);
window = exp(-(1/((FWHM*S_R)/(2*sqrt(log(2))))^2)*((1:FFT_length) - 0.5*FFT_length).^2)'*ones(1,total_No_spectra-1);
for shift_loop = 1:(Multiplicateur-1)
    
    C = sound((1+LL):((total_No_spectra-1)*FFT_length+LL));
    LL = LL + round(FFT_length/Multiplicateur);
    B = abs(fft(window.*reshape(C,FFT_length,(total_No_spectra-1))));
    % Store one half of the points (because the power spectrum is
    % symmetrical)
    temp = B(1:floor(size(B,1)/2),:);      
    % Perform the spectrogram:
    %Sp_with_overlap(:,(shift_loop+1):Multiplicateur:(shift_loop + (total_No_spectra-1)*(Multiplicateur-1))) = temp;
    Sp_with_overlap(:,(1+shift_loop):Multiplicateur:((total_No_spectra-1+shift_loop)+(total_No_spectra-2)*(Multiplicateur-1))) = temp;
    
end


y = Sp_with_overlap;



