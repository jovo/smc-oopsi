% the format of the waveform files to be read by Adam and Mor's program
% is as follows
% 
% Extension: .dat
% Precision: double
% Machine Format: Little Endian
% Vector: Single horizontal row, 1 data point = 1 time point = 1/5000 second
% Gain: "1" in output = 1 picoamp to be injected
% magnitude should range between 0 and 1

% Example code:

clear all, close all, clc
data = [-11000:-1100 1100:11000];
fid = fopen('testbin5.dat','wb','ieee-le')
fwrite(fid,data,'double');
fclose(fid)
