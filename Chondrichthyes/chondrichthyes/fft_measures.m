function[f,abs_FFTX] = fft_measures(x,fs,z)   
% Input: Signal, frequency sample and an integer 0 or 1 to determine if a
% filter is applied to the output spectrum.
% Output: Measures to be used to plot the spectrum (f and originated
% spectrum).
% -----------------------------------------------------------------------
% [f,spectrum] = fft_measures(signal,frequency sample, 0 or 1)

fn=fs/2;				% Nyquist frequency
NFFT=2^(nextpow2(length(x)));
NumUniquePts= ceil((NFFT+1)/2);	% Unique points
FFTX = calculate_fft(x);
FFTX=FFTX(1:NumUniquePts);		% DFT unique points
abs_FFTX= abs(FFTX);		% DFT absolute value
abs_FFTX= abs_FFTX * 2;		% auxiliary calculus
abs_FFTX(1)= abs_FFTX(1);		% DC frequency is unique
if ~rem(NFFT,2)			% if odd, break
abs_FFTX(length (abs_FFTX))= abs_FFTX(length(abs_FFTX))/2;		% fn is unique
end
abs_FFTX= abs_FFTX/length(x);	% DFT absolute value
if z == 1
    abs_FFTX = my_filter(abs_FFTX);
else
    abs_FFTX = abs(abs_FFTX).^2; % power spectrum
end
f=(0:NumUniquePts-1)*2*fn/NFFT;	% DFT frequency vector with NFFT/2 samples

