%% Purpose
% The purpose of this script is to use Near Field Focusing (NFF) to
% generate an image of a given 3D volume. 

%% Inputs
% E: This vector contains the complex electric field values at each of the
%    RX antennas and for each of the frequencies used for imaging. The
%    dimensions of this vector is NF x NR. Where NF is the number of
%    frequencies in the sweep and NR is the number of RX antennas in the
%    array

% txLoc: This vector contains the locations of the TX antennas for imaging.
%        The dimensions of this array are 1 x 3. For this version it is
%        easiest to assume that each set of RX data is associated with only
%        one TX antenna

% rxLoc: This vector contains the locations of the RX Antennas for imaging.
%        The dimension of this array are NR x 3. NR is the number of RX
%        antennas in the array associated with the single TX Antenna. Note
%        that the field Data associated with rxLoc(i,:) for all frequencies
%        is E(:,i). For a specific frequency index j and rx i, the field
%        would be E(j,i)

% lamba: This vector contains the wavelength for each of the frequencies in
%        the imaging sweep. The dimensions of this variable are NF x 1.
%        Where NF is the number of frequencies in the sweep

% nCells: This 3 element vector contains the number of imaging cells in
%         each of the 3 cartesian directions x,y,z. nCells(1) gives the
%         number of cells in the x direction (NX), nCells(2) gives the
%         number of cells in the y direction (NY), ncells(3) gives the
%         number of cells in the z direction (NZ).

% ix: This array contains the x location of each image cell. It has
%     dimensions of NX x NY x NZ. 

% iy: This array contains the y location of each image cell. It has
%     dimensions of NX x NY x NZ. 

% iz: This array contains the z location of each image cell. It has
%     dimensions of NX x NY x NZ. 
%% Outputs

%% Script
function [I,dRnd] = AbsImageVector(E,txLoc,rxLoc,lambda,nCells,ix,iy,iz)
    tic
    NR = size(rxLoc,1); %Number of receive antennas per tx antenna location
    NT = size(txLoc,1); %Number of transmist antennas
    NF = size(lambda,1); %Number of testing frequencies
    NC = numel(ix); %number of cell points (NX*NY*NZ)
    
    %% Calculate Distance From TX to Image Cell To RX
    % Converts the Array of locations into a matrix of XYZ locations has
    % dimensions of (NC) x 3
    cellLoc = [reshape(ix,1,[],1).', reshape(iy,1,[],1).', reshape(iz,1,[],1).'];
    % Distance from TX to each of the image cells dimensions of 1 x NC
    dTX_cell = sqrt(sum((cellLoc-txLoc).^2,2)).';
    cosdAng = (sqrt(cellLoc(:,3).^2)./sqrt(cellLoc(:,1).^2+cellLoc(:,2).^2+cellLoc(:,3).^2)).';
    %cosdAng = cosdAng+0.4;
    cosdAng(:)=1;
    
    %Distance from each RX to each of the cells, dimensions of NR x NC
    dRX_cell = squeeze(sqrt(sum((repmat(permute(cellLoc,[3,2,1]),NR,1)-repmat(rxLoc,1,1,NC)).^2,2)));
    %Total Distance traveled to go from the TX to a specific image cell and
    %back to a specific RX dimensions of NR x NC
    dRnd = dRX_cell+repmat(dTX_cell,NR,1);  
    cosdAng = repmat(cosdAng,NR,1);
    
    %% Calculate the wavenumber K (radians/meter)
    k = 2*pi./(lambda); %dimensions of NF x 1
    k = permute(k,[3,2,1]); %dimensions of 1 x 1 x NF
    
    %% Calculates the phase change for each rnd trip at each frequency
    phsChng = dRnd.*k; %phase change in radians dimensions NR x NC x NF
    cField = exp(1i*phsChng); %Calculates field using just phase
    cField = 4*pi*repmat(1./(cosdAng.^2),1,1,NF).*repmat(dRnd,1,1,NF).*cField; %Scales field by the distance travelled
    %cField = 4*pi*cField;
    I = squeeze(sum(sum((permute(cField,[3,1,2]).*repmat(E,1,1,NC)),1),2));
    I = reshape((I),nCells);
    %dRnd = reshape(squeeze(dRnd),nCells);
end