%% FAUSTINE'S FILE
%% PRODUCE AN EXCEL FILE WITH DATA NECESSARY FOR COMPARING JEFFERY VS. BROWNIAN MOTION
%~~~~ VARIABLES THAT ARE CALCULATED ARE:
%
% * In the xy plane:
%
% ** Projected length Lp (in microns) based on coordinates (can be improved using the least squares method)
% ** Projected arc length (in microns), to compare with Lp
% ** Angle phi : atan(delta y / delta x)
% ** nx and ny: deduced from Lp
%
% * In the xz plane (out-of-plane):
%
% ** Angle theta: depends on Lp, average filament, length, and aspect ratio of the filament (length / helix pitch), according to Gao et al., Phys. Rev. Lett. 2015.
% ** nz: can be approximated as nz ~ sin(theta)
%
% * INDICATORS TO COMPARE JEFFERY VS. BROWNIAN MOTION
% ** Jeffery constant C (depends on nx, ny, nz, and the filament aspect ratio)
% ** Modified Jeffery constant Cm = sign(C)/(1 + abs(C))
% ** Autocorrelation function of Cm --> gives a Jeffery decay time tau
% ** Jeffery oscillation period tJ (depends on the filament aspect ratio and the shear rate)

%% CODE
% Name of the file when created (created in the code folder by default)
filename = 'additionaldata.xlsx';

% Info on the filament in microns
% * equivalent of Jeffery width ("small" radius of the ellipsoid) = helix pitch
pitch = 0.879;
% * equivalent of Jeffery length ("long" radius of the ellipsoid) = length of a straight line going from one end of the filament to the other
longueur = 8.301;
% * Aspect ratio
lambda = longueur/pitch;

% Variables which initial value need to be set up before launching the loop
% * Average length on all frames Lmean
Lmean = 0;

% Incremental variables calculations
for j = 1 : xy.nframe
% Projected length Lp in the xy plane in microns
    Lp(j) = sqrt((xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1))^2+(xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))^2);
    Lpinmicrons(j) = Lp(j)/10.24; %10.24 is the conversion factor pixel/degrees for imgs of 1024x1024 of 100x100microns
% Projected arc length Arclen in the xy plane in microns
    arclenMic(j) = xy.arclen(j)*100/1024;
% Angle Phi (projection in the xy plane) in degrees
    phi(j) = atan( (xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))/(xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1)) );
    phiindeg(j) = phi(j) * 180 / pi;
% Angle Theta (approximation of the projection in the xz plane), using Gao et al., Phys. Rev. Lett. 2015
% * Average length Lmean on all frames of a given experimental set
    Lmean = (xy.arclen(j) + Lmean)/xy.nframe;
% * Theta angle based on Gao et al.
    theta(j) = asin( ((lambda * Lp(j)) / (Lmean - 1)) / (lambda - 1));
% Z coordinate approximated by the theta angle according to Gao et al.
    nz(j) = sin(theta(j));
% X and Y coordinates approximated using Lp
    nx(j) = xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1);
    ny(j) = xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2);
% Jeffery constant C and modified constant Cm
    C(j) = sqrt(nx(j)^2 + (nz(j)^2/lambda^2))/ny(j);
    Cm(j) = sign(C(j))/(1+abs(C(j)));
end

% Writing out the variables of interest in the Excel file
writematrix(Lpinmicrons,filename,'Sheet',1);
writematrix(arclenMic,filename,'Sheet',2);
writematrix(phiindeg,filename,'Sheet',3);
writematrix(nz,filename,'Sheet',4);
writematrix(C,filename,'Sheet',5);
writematrix(Cm,filename,'Sheet',6); 

% Giving a name to the sheets that allows recognizing easily the variables they stand for
e = actxserver('Excel.Application'); % # open Activex server
ewb = e.Workbooks.Open('C:\Users\Faustine\Documents\POSTDOC\Image treatment\Francesco - Matlab\Modified code\additionaldata.xlsx'); % # open file (enter full path!)
ewb.Worksheets.Item(1).Name = 'Length (microns)'; % # rename 1st sheet
ewb.Worksheets.Item(2).Name = 'Arclen (microns)'; % # rename 2nd sheet
ewb.Worksheets.Item(3).Name = 'Phi (xy) in deg'; % # rename 3rd sheet
ewb.Worksheets.Item(4).Name = 'Nz (coordinate)'; % # rename 4th sheet
ewb.Worksheets.Item(5).Name = 'C (Jeff. ct)'; % # rename 5th sheet
ewb.Worksheets.Item(6).Name = 'Cm (modif. Jeff. ct)'; % # rename 6th sheet
ewb.Save % # save to the same file
ewb.Close(false)
e.Quit

% 3D plot with the centroid coordinates (x,y) in microns as a function of the frame number
%for j = 1 : length(xy.frame)
%    CENTROID_X(j) = xy.centroid{j}(1,1)*100/1024;
%    CENTROID_Y(j) = xy.centroid{j}(1,2)*100/1024;
%end
%plot3(xy.frame,CENTROID_X,CENTROID_Y);
% * Allows to space the tick units equally along each axis
%axis equal
% * Axes labels
%xlabel('Frame nb')
%ylabel('Centroid X')
%zlabel('Centroid Y')