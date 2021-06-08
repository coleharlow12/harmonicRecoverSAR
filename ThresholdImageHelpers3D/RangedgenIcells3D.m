function [ix,iy,iz,nCells] = RangedgenIcells3D(CellDim,SpaceStart,SpaceEnd)

%Generates the x,y,z locations for each of the image cells. The variables
%are as follows
    %(1)  SpaceDim: 
    %Holds the x,y,z lengths of the imaging space. The center of the
    %rectangular prism-shaped imaging space is x=0, y=SpaceDim(2)/2, z = 0.
    %This was done so that the antennas are always located on one end of
    %the prism's imaging space
    
    %(2) CellDim:
    %Holds the x,y,z lengths of each image cell. SpaceDim./CellDim must
    %always be an integer
    
    %(3) plnCfg:
    %Tells what plan we will be imaging in plnCfg(1) as well as
    %The offset in the imaging plane plnCfg(2)
    SpaceDim = SpaceEnd-SpaceStart;

    nCells = round(SpaceDim./CellDim)+1; %The total cells in the x,y,z directions
    hCells = round(SpaceDim./CellDim)/2; %Cells needed to cover half the image in each direction 
                                    %(not inluding origin cell)
                                    
    centX = SpaceStart(1)+SpaceDim(1)/2;
    centY = SpaceStart(2)+SpaceDim(2)/2;   %Y-coord in the center of the image plane
    centZ = SpaceStart(3)+SpaceDim(3)/2;  %Z-coord in the center of the image plane

    ix =[(centX-CellDim(1)*hCells(1)):CellDim(1):(centX+CellDim(1)*hCells(1))]; % Creates x vector with all x coordinates
    ix = repmat(ix',1,nCells(2),nCells(3)); %Creates a matrix so that the x value is constant over columns but varies with row

    iy = [(centY-CellDim(2)*hCells(2)):CellDim(2):(centY+CellDim(2)*hCells(2))]; %Creates y vector holding all y coordinates
    iy = repmat(iy,nCells(1),1,nCells(3)); %Creates a matrix so that the y value varies over columns and is constant across rows

    iz = [(centZ-CellDim(3)*hCells(3)):CellDim(3):(centZ+CellDim(3)*hCells(3))]; %Creates z vector holding all z coordinates
    iz = permute(iz,[3,1,2]); %Rotates iz so that z varies along the 3rd dimension
    iz =repmat(iz,nCells(1),nCells(2),1); %Repeats z so that it is a 3D object which varies along
