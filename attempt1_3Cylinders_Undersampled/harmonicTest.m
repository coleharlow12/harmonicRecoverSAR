%% Assumptions
% If chirping multiple antennas I chirp lowest number to highest i.e (TX0
% then TX1 then TX2) or (TX1 then TX2) or (TX0 then TX1). This is how the
% algorithm assumes it is given data

clc
clear
close all

%% Addpath for Imaging Functions
curDir = pwd; %Gets directory of current file
curDir1 = strsplit(curDir,'/');
helpPath = string(strjoin(curDir1(1:(end-1)),'/'))+'/ThresholdImageHelpers3D';
addpath(helpPath);

%%
measFold = string(curDir)+'/Measurement47'; %Appends location of data files
addpath(measFold); %Adds measurement data to the path

useBackSub = false;
backFold = string(curDir)+'/Measurement74'; %Appends location of background data
if useBackSub, addpath(backFold); end %Adds background data to current path


locs = readtable('measLocs.xlsx'); %Loads the location of each file
locs = locs{:,:}; %Converts from table to array
measLocs = locs(2:end,:); %Gets rid of column number from excel
measLocs = (measLocs-measLocs(1,:))/1000; %Make it so first stop is at the origin and mm-> m

measLocs(:,1) = -measLocs(:,1); %Positive x for printer is negative x in image space
measLocs(:,2) = measLocs(:,3); %Z for printer is y for image space
measLocs(:,3) = 0;

%% Skip (Creates Smaller Aperture)
skipVertBottom = 3; %Makes it so only one row is used
skipVertTop = 4;
skipHorRight = 0;
skipHorLeft = 0;

numHor = length(unique(measLocs(:,1)));
numVert = length(unique(measLocs(:,2)));

if ((skipVertBottom+skipVertTop) >= numVert || (skipHorRight+skipHorLeft) >= numHor )
    fprintf('You cannot skip that many Rows/Columns: Quitting Execution')
   return
end

if (numHor*numVert ~= length(measLocs))
    fprintf('You have Missing Measurements: Quitting Execution')
    return
end

fileInd = 1:length(measLocs); %All of the files available if I used them all
fileIndS = reshape(fileInd,numHor,numVert)';
fileIndS = fileIndS((1+skipVertBottom):(numVert-skipVertTop),(1+skipHorRight):(numHor-skipHorLeft)); %Removes the files I am not interested in using
fileInd = reshape(fileIndS',1,[]);

figure(1)
plot(measLocs(fileInd,1)*100,measLocs(fileInd,2)*100,'r^')
xlabel('x-axis(cm)')
ylabel('y-axis(cm)')

%% Selet TX to Enable
txen = [1 0 0]; %Specifies which TX to use in imaging

tx1RxEn = [3]; %Specifies which Rx data to use with each tx when imaging
tx2RxEn = [];
tx3RxEn = [];

numTXchirped = 3; %The number of TX antennas chirped in each measurement
numTXen = sum(txen); %Calculates the number of tx we are using for imaging

%                 RX4 RX1
% TX2    TX1      RX3 RX2
% 
% TX3

%% Radar Parameters Used For Measurement IWR6843

fstrt = 60.25e9; %starting frequency Hz
slope = 144.984e12; %Hz/s
TxStrtTm = 0; %Transmit start time in seconds
AdcStrtTm = 5.28e-6; %ADC start time in seconds
numSamps = 200; %The number of samples taken per chirp
fs = 10000*10^3; %The sample rate in Hz
rampEndTm = 25.80e-6; %The ramp end time (s)

%Rx locations in mil
rxLoc = [927.934,1677.027,0;...
         927.934,1581.5,0;...
         1023.457,1581.5,0;...
         1023.457,1677.027,0];

%Tx locations in mil
txLoc = [1509.926,1576.600,0;...
         1700.976,1576.600,0;...
         1700.976,1385.559,0];

txLoc = txLoc*0.0000254; %Converts mils to meters
rxLoc = rxLoc*0.0000254; %Converts mils to meters

rxLoc = rxLoc-txLoc(1,:); %Places txLoc1 at the origin
txLoc = txLoc-txLoc(1,:);

rxLoc(:,1) = -rxLoc(:,1);
txLoc(:,1) = -txLoc(:,1);


%_________________Constants Used for Imaging_________________________
c = 3e8;
ts = AdcStrtTm:(1/fs):(AdcStrtTm+1/fs*(numSamps-1));
fTxS = fstrt+ts*slope;
lambda = c./fTxS'; %Calculate the wavelength

%_________________Used to Load the Data______________________________
base = 'ADCdata'; %The base name of all measured data
fnum = 0; %Selects file number to load
append = '_Raw_0.bin'; %Used for loading data
chirpSelect=0; % Of all the chirps measured for a tx which to select

arrayVal = zeros(1,length(fileInd));
arrayInd = 1;

sNum = 30; %The sample number to use

tic
for iM = fileInd
    %load data
    %Load the imaging Data
    fnum = iM-1; %Save files are indexed from 0 
    
    fprintf('Load Data %.0f \n',fnum)
    loadFile = measFold+'/'+base+string(fnum)+append;
    measDataAll = readDCA1000(loadFile);
    if useBackSub, loadBack = backFold+'/'+base+string(fnum)+append; end %loaded background data
    if useBackSub, backDataAll = readDCA1000(loadBack); end %subtract background data
    
    %Account for phase shift in radar architecture
    measDataAll([2,3],:)=exp(1i*pi)*measDataAll([2,3],:);
    if useBackSub, backDataAll([2,3],:)=exp(1i*pi)*backDataAll([2,3],:); end   
    if useBackSub,measDataAll=measDataAll-backDataAll;end 
    
    for txI = 1:length(txen)
        if txen(txI)
            if(txI==1)
                rxEn = tx1RxEn;
            end
            measData = measDataAll(rxEn,((1:numSamps)+(numSamps*(chirpSelect*numTXchirped+(txI-1)))));
            %Smooth the data using moving average filter
            mIm = 1i*smoothdata(imag(measData),2,'movmean',90);
            mRe = smoothdata(real(measData),2,'movmean',90);
            measData= measData-mIm-mRe; 
            
            arrayVal(arrayInd) = measData(sNum); %Takes the first sample
            arrayInd = arrayInd+1;
        end
    end
end

%% Determine Frequencies 
k = 2*pi*fTxS(sNum)/c;
d = 2.5e-3
fs_x = c/d %Spatial Sampling Frequency
b = -[arrayVal(3);arrayVal(4)];
A = [arrayVal(2) arrayVal(1); ...
     arrayVal(3) arrayVal(2);];
 
 % Solve for b1 and b2
 h1 = linsolve(A,b);    %h1 and h2 Give the same answer (not surprising)
 h2 = pinv(A'*A)*A'*b;
 
 h1 = [1;h1];     %1z^(0) + h[1]z^-1 + h[2]z^-2 => 1z^2 + h[1]z + h[2]
 rb = roots(h1);  %Finds the roots of rb
 wr = angle(rb);  %Angles should correspond to normalized frequencies [0 2pi] of signals
 
 fspat = wr/(2*pi)*fs_x;
 fprintf('A signal with frequency %.3f GHz is present \n',fspat/1e9)
 
%% Calculate Angle Given the Frequency
theta = asin(fspat*pi/(k*c))*180/pi
toc

%% Create an FFT of the Data
M = 200;
fftArray = fftshift(fft(arrayVal,M));
fspatfft = (-M/2:(M-1)/2)/M*fs_x;

plot(fspatfft/1e9,abs(fftArray))
xlabel('GHz')