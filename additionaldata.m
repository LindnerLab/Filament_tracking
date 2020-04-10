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
% ** Rotational diffusion time tau_r
% ** Jeffery oscillation period tJ (depends on the filament aspect ratio and the shear rate)

%% CODE
% Name of the file when created (created in the code folder by default)
filename = 'additionaldata.xlsx';

% INFO ON THE FILAMENT IN MICRONS
% * equivalent of Jeffery width ("small" radius of the ellipsoid) = helix diameter (between the top of one side of the helix to the other side of the helix)
diameter = 0.879;
% * equivalent of Jeffery length ("long" radius of the ellipsoid) = length of a straight line going from one end of the filament to the other
fil_length = 8.301;
% * Aspect ratio
lambda = fil_length/diameter;
% * Average length based on arclen in microns (arclenMic)
for j = 1 : xy.nframe
    % Projected arc length Arclen in the xy plane in microns
    arclenMic(j) = xy.arclen(j)*100/1024;
end
Lmean = sum(arclenMic)/xy.nframe;

% INCREMENTAL VARIABLES CALCULATION
for j = 1 : xy.nframe
% Projected length Lp in the xy plane in microns
    Lp(j) = sqrt((xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1))^2+(xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))^2);
    Lpinmicrons(j) = Lp(j)/10.24; %10.24 is the conversion factor pixel/degrees for imgs of 1024x1024 of 100x100microns
% Angle Phi (projection in the xy plane) in degrees
    phi(j) = atan( (xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))/(xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1)) );
    phiindeg(j) = phi(j) * 180 / pi;
% Angle Theta (approximation of the projection in the xz plane), using Gao et al., Phys. Rev. Lett. 2015
    theta(j) = asin( ((lambda * Lp(j)) / (Lmean - 1)) / (lambda - 1));
% Z coordinate approximated by the theta angle according to Gao et al.
    nz(j) = sin(theta(j));
% X and Y coordinates approximated using Lp
    nx(j) = xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1);
    ny(j) = xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2);
% JEFFERY CONSTANT C AND MODIFIED CONSTANT Cm
    C(j) = sqrt(nx(j)^2 + (nz(j)^2/lambda^2))/ny(j); %!!!! need to fix the time instead of the frame nb
    Cm(j) = sign(C(j))/(1+abs(C(j))); %!!!! need to fix the time instead of the frame nb
end

% AUTOCORRELATION FUNCTION autocorr(function,'Numlags',number of lags)
% * Y coordinates stored within a variable
% * X coordinates are 1,2,3... as numerous as the number of lags since the function is discretized
Ycorr = real(autocorr(Cm,'Numlags',40)); %real() is the real part of Cm, which are complex numbers
Xcorr = (1:41);
% * For the fit, Xcorr and Ycorr must be column vectors
T_Ycorr = transpose(Ycorr);
T_Xcorr = transpose(Xcorr);
% Exponential fit for a simple exponential a * exp(b * t) = a * exp(-t/tau)
% For fits, coordinates must be column vectors
expofit = fit(T_Xcorr,T_Ycorr,'exp1');
plot(expofit,T_Xcorr,T_Ycorr);
% Harvesting expofit coefficients
a = expofit.a;
tau = -1/expofit.b;

% JEFFERY OSCILLATION PERIOD
%!!!! Need to determine gammadot for my experiments
gammadot = 18; %shear rate (in s-1); 18 s-1 is the value given by Z�ttl et al.
tJ = (2*pi*(lambda + 1/lambda))/gammadot; %Jeffery oscillation period (according to Z�ttl et al., 2019)

% ROTATIONAL DIFFUSION TIME 
kb = 1.38064852 * 10^(-23); % Boltzmann constant in m2 kg s-2 K-1
T = 25 + 273.15; % Temperature in Kelvin
eta = 10^-3; % dynamic viscosity in Pa.s (dynamic viscosity of water as a first pass)
R = 8 * 10^(-6); % radius of an equivalent sphere in m (taken: length of the filament; !!!! can be optimized if the equivalent radius is calculated!)
%Drsphere = (kb * T) / (8*pi*eta*R^3); % rotational diffusion coefficient (Stokes-Einstein relationship)
% * Adaptation of the Stokes-Einstein relationship by Perrin (1936) to prolate bodies 
a =(fil_length/2)*10^(-6); % half-length of the filament in m
b = (diameter/2)*10^-6; % half-width of the filament in m
p = a/b; %new aspect ratio needed for Nuris' diffusion coefficient formula (cf. her thesis p.47).
%Dr = (kb * T) / (6*pi*eta*(((a^2-b^2)^(1/2)) / log((a+(a^2-b^2)^(1/2))/b^2))); % Saverio's adapted diffusion coefficient to prolates (cf. p. 40 of my notebook) 
Dr = (kb * T) / (6*eta*V*g); % Nuris' formula (eq. 10)
V = (4*pi*a*b^2)/3; % Volume of a prolate ellipsoid (Nuris' thesis)
S = (1 / sqrt(p^2-1)) * log(p+sqrt(p^2-1)); % Nuris' formula (eq. 12)
g = (2*(p^4-1)) / (3 * p * ((2*p^2-1)*S - p)); % Nuris' formula (eq. 11)
tau_r = 1/(2*Dr); % rotational diffusion time

% WRITING OUT THE VARIABLE OF INTEREST IN AN EXCEL FILE
% * DATA (SHEET 1)
% 
% ** Column titles for the Data sheet
xlswrite(filename,{'Lp (�m)'},'Feuil1','A1');
xlswrite(filename,{'Arclen (�m)'},'Feuil1','B1');
xlswrite(filename,{'Phi (deg)'},'Feuil1','C1');
xlswrite(filename,{'nz (�m)'},'Feuil1','D1');
xlswrite(filename,{'Jeff. C'},'Feuil1','E1');
xlswrite(filename,{'Modif. Jeff. Cm'},'Feuil1','F1');
xlswrite(filename,{'Expofit coeff a'},'Feuil1','G1');
xlswrite(filename,{'Expofit coeff tau'},'Feuil1','H1');
xlswrite(filename,{'Jeff. period tJ'},'Feuil1','I1');
xlswrite(filename,{'Rot. diff. time tau_r'},'Feuil1','J1');
% ** Writing out data by columns on the Excel sheet 1
% Using transpose() because data are by default organized in rows instead of columns
writematrix(transpose(Lpinmicrons),filename,'Sheet',1, 'Range', 'A2');
writematrix(transpose(arclenMic),filename,'Sheet',1, 'Range', 'B2');
writematrix(transpose(phiindeg),filename,'Sheet',1, 'Range', 'C2');
writematrix(transpose(nz),filename,'Sheet',1, 'Range', 'D2');
writematrix(transpose(C),filename,'Sheet',1, 'Range', 'E2');
writematrix(transpose(Cm),filename,'Sheet',1, 'Range', 'F2');
writematrix(a,filename,'Sheet',1,'Range','G2');
writematrix(tau,filename,'Sheet',1,'Range','H2');
writematrix(tJ,filename,'Sheet',1,'Range','I2');
writematrix(tau_r,filename,'Sheet',1,'Range','J2');

%
% * PARAMETERS (SHEET 2)
%
% ** Row titles for the Parameters sheet
xlswrite(filename,{'FRAMES INFO'},'Feuil2','A1');
xlswrite(filename,{'First frame treated'},'Feuil2','A2');
xlswrite(filename,{'Frame step'},'Feuil2','A3');
xlswrite(filename,{'Final frame treated'},'Feuil2','A4');
xlswrite(filename,{'Nb of treated frames'},'Feuil2','A5');
xlswrite(filename,{'Nb of empty frames'},'Feuil2','A6');
%
xlswrite(filename,{'FLAGELLUM INFO'},'Feuil2','A7');
xlswrite(filename,{'Diameter (�m)'},'Feuil2','A8');
xlswrite(filename,{'Length (�m)'},'Feuil2','A9');
xlswrite(filename,{'Aspect ratio lambda'},'Feuil2','A10');
xlswrite(filename,{'Average length Lmean (�m)'},'Feuil2','A11');
%
xlswrite(filename,{'CODE PARAMETERS'},'Feuil2','A12');
xlswrite(filename,{'basepath'},'Feuil2','A13');
xlswrite(filename,{'tifname'},'Feuil2','A14');
xlswrite(filename,{'Nb of filaments'},'Feuil2','A15');
xlswrite(filename,{'Fibermetric thickness (px)'},'Feuil2','A16');
xlswrite(filename,{'Fibermetric structsensitivity'},'Feuil2','A17');
xlswrite(filename,{'Gaussian blur noise lengthscale lnois (px)'},'Feuil2','A18');
xlswrite(filename,{'Gaussian blur object size lobject (px)'},'Feuil2','A19');
xlswrite(filename,{'Gaussian blur threshold'},'Feuil2','A20');
xlswrite(filename,{'Binarization sensitivity'},'Feuil2','A21');
xlswrite(filename,{'Skeletonization MinBranchLength (px)'},'Feuil2','A22');
xlswrite(filename,{'Bspline ds'},'Feuil2','A23');
xlswrite(filename,{'Bspline npnts'},'Feuil2','A24');
%
% * Writing out parameters in rows
% ** First frame
writematrix(xy.frame(1),filename,'Sheet',2,'Range','B2');
% ** Frame step chosen
for j = 1 : xy.nframe
    if xy.frame(j) == j
        if ismember(j+1,xy.emptyframe) == 0 %ismember tells if the j+1 element belongs to the xy.frame array
            writematrix(xy.frame(j+1)-xy.frame(j),filename,'Sheet',2,'Range','B3');
        end
    end %!!!! sometimes the loop ends and no value was found for the step if the first if is not verified. Need to improve this.
end
% Last frame treated (this number + nb of empty frames = total number of frames asked to treat)
writematrix(xy.frame(xy.nframe),filename,'Sheet',2,'Range','B4');
% Nb of frames treated
writematrix(xy.nframe,filename,'Sheet',2,'Range','B5');
% Nb of empty frames
writematrix(length(xy.emptyframe),filename,'Sheet',2,'Range','B6');
%
writematrix(diameter,filename,'Sheet',2,'Range','B8');
writematrix(fil_length,filename,'Sheet',2,'Range','B9');
writematrix(lambda,filename,'Sheet',2,'Range','B10');
writematrix(Lmean,filename,'Sheet',2,'Range','B11');
%
writematrix(basepath,filename,'Sheet',2,'Range','B13');
writematrix(tifname,filename,'Sheet',2,'Range','B14');
writematrix(FilNum,filename,'Sheet',2,'Range','B15');
writematrix(thickness,filename,'Sheet',2,'Range','B16');
writematrix(structsensitivity,filename,'Sheet',2,'Range','B17');
writematrix(lnoise,filename,'Sheet',2,'Range','B18');
writematrix(lobject,filename,'Sheet',2,'Range','B19');
writematrix(threshold,filename,'Sheet',2,'Range','B20');
writematrix(sensitivity,filename,'Sheet',2,'Range','B21');
writematrix(MinBranchLength,filename,'Sheet',2,'Range','B22');
writematrix(ds,filename,'Sheet',2,'Range','B23');
writematrix(npnts,filename,'Sheet',2,'Range','B24');
%
% GIVING NAMES TO THE EXCEL SHEETS
e = actxserver('Excel.Application'); % # open Activex server
ewb = e.Workbooks.Open('C:\Users\Faustine\Documents\POSTDOC\Image treatment\Francesco - Matlab\Modified code\additionaldata.xlsx'); % # open file (enter full path!)
ewb.Worksheets.Item(1).Name = 'Data'; % # rename 1st sheet
ewb.Worksheets.Item(2).Name = 'Parameters'; % # rename 2nd sheet
ewb.Worksheets.Item(2).Range('A1:B1').Interior.Color=hex2dec('F0F4C3'); % # color row A1
ewb.Worksheets.Item(2).Range('A7:B7').Interior.Color=hex2dec('F0F4C3'); % # color row A6
ewb.Worksheets.Item(2).Range('A12:B12').Interior.Color=hex2dec('F0F4C3'); % # color row A11
ewb.Save;
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