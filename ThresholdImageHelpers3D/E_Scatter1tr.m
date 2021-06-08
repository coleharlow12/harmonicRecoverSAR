function [Er,drnd] = E_Scatter1tr(pt,pr,ps,lambda,sCoeff)
% Inputs
%   pt - Ntx3 matrix containing x-, y-, and z- coordinates of Nt
%        transmitter antennas
%   pr - Nrx3 matrix containing x-, y-, and z- coordinates of Nr receiver antennas
%   ps - Nsx3 matrix containing x-, y-, and z- coordinates of Ns target
%        scatter points
%   lambda - Nfx1 array containing free space wavelengths for Nf frequency
%            points
%   taper_rx - array containing array taper coefficients for RX array
%   taper_tx - array containing array taper coefficients for TX array
%
% Outputs
%   Er - NfxNrxNt matrix containing complex values for scattered
%        electric fields for each RX and TX antenna and frequency summed over all target scatter points 
%        COLE: Seems like it is the electric field at antenna RX which
%        originated at antenna TX and scattered off of a specific scatter
%   Ar - 1xNr matrix containing amplitudes for the scattered
%        electric fields for each receiver antenna summed over all target scatter points 
%
% E_Scatter2 calculates the scattered electric field at each RX and TX for
% all frequency points. The output Er is a matrix containing complex values 
% for each RX and TX and frequency. The total scattered field for each RX and TX 
% antenna is the sum of all scattered fields due to each target scatterer. 
% The output Ar is the amplitudes for the electric fields in Er.
%
% Taper applied to RX and TX arrays by multiplying amplitude with
% coefficients in taper_rx and taper_tx
numArgs = nargin;
% Initialize outputs
Ert = ones(length(lambda),size(pr,1),size(pt,1),size(ps,1));
Art = ones(1,size(pr,1),size(pt,1),size(ps,1));
Er = ones(length(lambda),size(pr,1),size(pt,1));
Ar = ones(1,size(pr,1),size(pt,1));

% Calculate the scattered field  at each RX and TX antenna for all 
% Frequency points
if(numArgs == 4)
    j = sqrt(-1);
    for xx=1:size(pt,1) % for loop over TX antennas
        for kk=1:size(pr,1) % for loop over RX antennas
            for tt=1:size(ps,1) % for loop over target scatter points
                % Vector dsr from RX to scatterer
                dsr = pr(kk,:) - ps(tt,:); %x,y,z coordinates to travel from 
                                           %the scatter to RX
                % Vector dts from TX to target scatterer
                dts = ps(tt,:) - pt(xx,:); %x,y,z coordinates to travel from 
                                           %the scatter to TX
                dtsr(1,kk,xx,tt)=norm(dts)+norm(dsr); % calculates the distance to go from 
                                            % a given TX to the scatterer and 
                                            % then back to the RX from the
                                            % scatterer
                %Amplitude note the 1/r dependence of radiation(Griffiths 467)                        
                Art(1,kk,xx,tt) = (1./(4*pi*dtsr(1,kk,xx,tt)));
                % Calculates the scattered electric field including the phase
                % Note this section also takes into account the wavelength due 
                % to its effect on the phase
                Ert(:,kk,xx,tt) = Art(1,kk,xx,tt).*exp((-j*2*pi./lambda).*dtsr(1,kk,xx,tt));
            end
        end
    end

    % Sum scattered fields and amplitudes over all target scatter
    % points
    Er = squeeze(sum(Ert,4)); %sums along dimension four and then removes any one dimensions
    drnd = squeeze(dtsr);
    
elseif(numArgs==5)
    j = sqrt(-1);
    for xx=1:size(pt,1) % for loop over TX antennas
        for kk=1:size(pr,1) % for loop over RX antennas
            for tt=1:size(ps,1) % for loop over target scatter points
                % Vector dsr from RX to scatterer
                dsr = pr(kk,:) - ps(tt,:); %x,y,z coordinates to travel from 
                                           %the scatter to RX
                % Vector dts from TX to target scatterer
                dts = ps(tt,:) - pt(xx,:); %x,y,z coordinates to travel from 
                                           %the scatter to TX
                dtsr(1,kk,xx,tt)=norm(dts)+norm(dsr); % calculates the distance to go from 
                                            % a given TX to the scatterer and 
                                            % then back to the RX from the
                                            % scatterer
                %Amplitude note the 1/r dependence of radiation(Griffiths 467)                        
                Art(1,kk,xx,tt) = sCoeff(tt)*(1./(4*pi*dtsr(1,kk,xx,tt)));
                % Calculates the scattered electric field including the phase
                % Note this section also takes into account the wavelength due 
                % to its effect on the phase
                Ert(:,kk,xx,tt) = Art(1,kk,xx,tt).*exp((-j*2*pi./lambda).*dtsr(1,kk,xx,tt));
            end
        end
    end

    % Sum scattered fields and amplitudes over all target scatter
    % points
    Er = squeeze(sum(Ert,4)); %sums along dimension four and then removes any one dimensions
    drnd = squeeze(dtsr);
    
end
end
    

