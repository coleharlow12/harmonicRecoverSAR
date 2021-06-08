% Created By: Cole Harlow (harlow.74@osu.edu)
% Purpose: This script take input data and applies a window to it.
% Inputs:
%           data: This is a vector which contains the data that is to be
%                 windowed. The input data should be of the format so that
%                 the windowed data is indexed along the first dimension
%                 (rows) of the matrix data. That is to say that the windowed 
%                 data is organized in columns. The length of that data is
%                 the number of rows in the matrix. Data can be a 3d matrix
%                 or 2d but 4d and above is not supported currently
%       
%        percHam: The current function takes the input data and applies a
%                 combined rectangular/hamming window to the data. The ends
%                 of the data are where the hamming window is applied while
%                 the middle of the data is where the rectangular window is
%                 applied. The percHam determines what percent of the data
%                 is windowed using the ham function. The range of this
%                 variable should be 0:1
%                 
%          disp:  This variable should be set to 1 if the window applied to
%                 should be displayed and 0 if the window shouldn't be
%                 applied
%                  
%        figNum:  If the data is to be displayed this is the figure number
%                 that will be used. 

function [data] = applyWindowTukey(data,fracHam,disp1,figNum)

    lengthData = size(data,1);  %Calculates the length of the data (how many rows)
    dim2 = size(data,2);        %Calculates the length of the 2nd dimension (how many columns)
    dim3 = size(data,3);        %Calculates the length of the 3rd dimension

    wind0 = tukeywin(lengthData,fracHam);

    if (disp1==1)
        n = 0:1:(lengthData-1);
        figure(figNum)
        plot(n,wind0)

        M = 10*lengthData;
        fftdata= fft(wind0,M);      %Calculates the FFT
        ifftdata= ifft(wind0,M);    %Calculates the IFFT
        fftdata = fftdata/max(fftdata); %normalization
        ifftdata = ifftdata/max(ifftdata); %normalization

        n = 0:1:(M-1);
        figure(figNum+1)
        plot(n,fftshift(20*log10(abs(fftdata))))
        title('FFT of window')
        axis([0 M-1 -30 0])

        figure(figNum+2)
        plot(n,fftshift(20*log10(abs(ifftdata))))
        title('IFFT of window')
        axis([0 M-1 -30 0])
    end

    for i = 1:dim2
        for j = 1:dim3
         data(:,i,j) = wind0.*data(:,i,j); 
        end
    end

end



