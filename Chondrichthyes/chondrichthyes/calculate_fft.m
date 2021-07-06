function[seq] = calculate_fft(seq)
% Input: Sequence in time domain
% Output: Sequence in frequence domain
% ------------------------------------
% seq = calculate_fft(seq)
seq = detrend(seq);
NFFT=2^(nextpow2(length(seq))); % Next highest power of 2 greater or equal to length(x)
seq = fft(seq,NFFT);
end