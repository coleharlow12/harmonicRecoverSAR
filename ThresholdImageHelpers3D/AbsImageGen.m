function [I] = AbsImageGen(Ern,txLoc,rxLoc,lambda,ix,iy,iz,nCells,plnCfg)
    NR = size(rxLoc,1); %Number of receive antennas per tx antenna location
    NT = size(txLoc,1); %Number of transmist antennas
    
    if(plnCfg(1) == 1) %xy plane
        for xIm = 1:nCells(1)
            for yIm = 1:nCells(2)
                for indTX = 1:NT
                    for indRX = 1:NR
                        %RX antenna k coordinates - location of the image cell
                        dir = rxLoc(indRX,:)-[ix(xIm,yIm),iy(xIm,yIm),iz(xIm,yIm)];
                        %image cell location - TX antenna xx coordinates
                        dti = [ix(xIm,yIm),iy(xIm,yIm),iz(xIm,yIm)]-txLoc(indTX,:);
                        %total distance from tx to images cell plus image cell to
                        %rx
                        dtir=(norm(dti)+norm(dir));
                        %Calculates phases of expected signal. If there is a
                        %distance term why is it not squared?
                        EI(:,indRX,indTX) = 4*pi*dtir*exp((1i*2*pi./lambda)*dtir);
                    end
                end
                %Calculate the multiplication of EI and the measured fields
                I(xIm,yIm)=(sum(sum(sum(Ern(:,:,:).*EI)))); %sum(A) sums across the first array dimension of A
            end
        end
      
    %IGNORE BELOW   
    elseif(plnCfg(1) == 2) %xz plane
         for xIm = 1:nCells(1)
            for zIm = 1:nCells(3)
                for indTX = 1:NT
                    for indRX = 1:NR
                        %RX antenna k coordinates - location of the image cell
                        dir = rxLoc(indRX,:)-[ix(xIm,zIm),iy(xIm,zIm),iz(xIm,zIm)];
                        %image cell location - TX antenna xx coordinates
                        dti = [ix(xIm,zIm),iy(xIm,zIm),iz(xIm,zIm)]-txLoc(indTX,:);
                        %total distance from tx to images cell plus image cell to
                        %rx
                        dtir=(norm(dti)+norm(dir));
                        %Calculates phases of expected signal. If there is a
                        %distance term why is it not squared?
                        EI(:,indRX,indTX) = dtir*exp((1i*2*pi./lambda).*dtir);
                    end
                end
                %Calculate the multiplication of EI and the measured fields
                I(xIm,zIm)=abs(sum(sum(sum(Ern(:,:,:).*EI)))); %sum(A) sums across the first array dimension of A
            end
         end
         
     elseif(plnCfg(1) == 3) %yz plane
         for yIm = 1:nCells(2)
            for zIm = 1:nCells(3)
                for indTX = 1:NT
                    for indRX = 1:NR
                        %RX antenna k coordinates - location of the image cell
                        dir = rxLoc(indRX,:)-[ix(yIm,zIm),iy(yIm,zIm),iz(yIm,zIm)];
                        %image cell location - TX antenna xx coordinates
                        dti = [ix(yIm,zIm),iy(yIm,zIm),iz(yIm,zIm)]-txLoc(indTX,:);
                        %total distance from tx to images cell plus image cell to
                        %rx
                        dtir=(norm(dti)+norm(dir));
                        %Calculates phases of expected signal. If there is a
                        %distance term why is it not squared?
                        EI(:,indRX,indTX) = dtir*exp((1i*2*pi./lambda).*dtir);
                    end
                end
                %Calculate the multiplication of EI and the measured fields
                I(yIm,zIm)=abs(sum(sum(sum(Ern(:,:,:).*EI)))); %sum(A) sums across the first array dimension of A
            end
         end  
    end
end