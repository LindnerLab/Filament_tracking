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
% ** Ratio tau/tJ (number of Jeffery oscillations a filament performs before losing memory about its Jeffery orbit state)

%% CODE
% Name of the file when created (created in the code folder by default)
filename = 'additionaldata.xlsx';

% FLAGELLUM INFO
% * diameter (equivalent of Jeffery width i.e. between the top of one side of the helix to the other side of the helix)
% * fil_length (length of a straight line going from one end of the filament to the other)
% * lambda = aspect ratio
diameter = 0.879;
fil_length = 8.301;
lambda = fil_length/diameter;
% * Average length based on the projected arclen in the xy plane in microns (arclenMic)
for j = 1 : xy.nframe
    arclenMic(j) = xy.arclen(j)*100/1024;
end
Lmean = sum(arclenMic)/xy.nframe;

% FRAME INFO
% * in_frame = First frame (real frame number)
% * fin_frame = Last frame (real frame number)
% * step = frame step (a frame is analyzed every ... frames)
% * total = chosen number of frames (not taking into account missing frames)
% * F Frequence ("frame per second") in Hz
% * FS frame step in seconds (time elapsed between two frames)
% * FSav average frame step in seconds (taking into account missing frames)
in_frame = xy.frame(1);
fin_frame = xy.frame(xy.nframe);
total = (xy.nframe + length(xy.emptyframe));
step = frame_step;
F = 30;
FS = 1/F*step;
FSav = (total / xy.nframe)*FS;

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
% Exponential fit for a simple exponential a * exp(b * k) where k is a nb of frames such as t = k*FSav (t=time)
% Therefore, converting frames in times, the time constant tau such as exp(-t/tau) -> tau = -FSav/b
expofit = fit(T_Xcorr,T_Ycorr,'exp1');
plot(expofit,T_Xcorr,T_Ycorr);
% Harvesting expofit coefficients
a = expofit.a;
tau = (-FSav/expofit.b);

% JEFFERY OSCILLATION PERIOD
%!!!! Need to determine gammadot for my experiments
% * gammadot is the shear rate (in s-1); 18 s-1 is the value given by Z�ttl et al (default value here)
% * tJ is the Jeffery oscillation period (in s) according to Z�ttl et al., 2019
% * Losingmemory is the number of Jeffery oscillations a filament performs before losing memory about its Jeffery orbit state
gammadot = 18; 
tJ = (2*pi*(lambda + 1/lambda))/gammadot;
Losingmemory = tau/tJ;

% ROTATIONAL DIFFUSION TIME FOR PROLATE BODIES
% * Boltzmann constant kb in m2 kg s-2 K-1
% * T Temperature in Kelvin
% * eta dynamic viscosity in Pa.s (dynamic viscosity of water used as a first pass)
kb = 1.38064852 * 10^(-23);
T = 25 + 273.15; 
eta = 10^-3; 
% * a half-length of the filament in m
% * b half-width of the filament in m
% * p new aspect ratio needed for Nuris' diffusion coefficient formula (cf. her thesis p.47).
% * V Volume of a prolate ellipsoid in m3 (Nuris' thesis)
% * g and S coefficients adapted to prolate ellipsoids (Nuris'formula eq. 11 and 12)
% * Dr diffusion coefficient for a prolate ellipsoid; Nuris' formula (eq. 10)
a =(fil_length/2)*10^(-6); 
b = (diameter/2)*10^-6; 
p = a/b;
V = (4*pi*a*b^2)/3; 
S = (1 / sqrt(p^2-1)) * log(p+sqrt(p^2-1));
g = (2*(p^4-1)) / (3 * p * ((2*p^2-1)*S - p));
Dr = (kb * T) / (6*eta*V*g);
% Rotation diffusion time
tau_r = 1/(2*Dr);

% WRITING OUT THE VARIABLES OF INTEREST IN AN EXCEL FILE
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
xlswrite(filename,{'Losing memory ratio tau/tJ'},'Feuil1','K1');
%
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
writematrix(Losingmemory,filename,'Sheet',1,'Range','K2');
%
% * FRAMES AND FLAGELLUM (SHEET 2)
%
% ** Row titles for the Frames and Flagellum sheet
xlswrite(filename,{'FRAMES INFO'},'Feuil2','A1');
xlswrite(filename,{'Frequence (Hz)'},'Feuil2','A2');
xlswrite(filename,{'Frame step (s)'},'Feuil2','A3');
xlswrite(filename,{'Average frame step (s)'},'Feuil2','A4');
xlswrite(filename,{'First frame treated'},'Feuil2','A5');
xlswrite(filename,{'Final frame treated'},'Feuil2','A6');
xlswrite(filename,{'Chosen number of frames'},'Feuil2','A7');
xlswrite(filename,{'Nb of treated frames'},'Feuil2','A8');
xlswrite(filename,{'Nb of empty frames'},'Feuil2','A9');
xlswrite(filename,{'FLAGELLUM INFO'},'Feuil2','A10');
xlswrite(filename,{'Diameter (�m)'},'Feuil2','A11');
xlswrite(filename,{'Length (�m)'},'Feuil2','A12');
xlswrite(filename,{'Aspect ratio lambda'},'Feuil2','A13');
xlswrite(filename,{'Average length Lmean (�m)'},'Feuil2','A14');
%
% Writing out data by rows on the Excel sheet 2
writematrix(F,filename,'Sheet',2,'Range','B2');
writematrix(FS,filename,'Sheet',2,'Range','B3');
writematrix(FSav,filename,'Sheet',2,'Range','B4');
writematrix(in_frame,filename,'Sheet',2,'Range','B5');
writematrix(fin_frame,filename,'Sheet',2,'Range','B6');
writematrix(total,filename,'Sheet',2,'Range','B7');
writematrix(xy.nframe,filename,'Sheet',2,'Range','B8');
writematrix(length(xy.emptyframe),filename,'Sheet',2,'Range','B9');
%
writematrix(diameter,filename,'Sheet',2,'Range','B11');
writematrix(fil_length,filename,'Sheet',2,'Range','B12');
writematrix(lambda,filename,'Sheet',2,'Range','B13');
writematrix(Lmean,filename,'Sheet',2,'Range','B14');
%
% * CODE PARAMETERS (SHEET 3)
%
% ** Row titles for the Code Parameters sheet
xlswrite(filename,{'CODE PARAMETERS'},'Feuil3','A1');
xlswrite(filename,{'basepath'},'Feuil3','A2');
xlswrite(filename,{'tifname'},'Feuil3','A3');
xlswrite(filename,{'Nb of filaments'},'Feuil3','A4');
xlswrite(filename,{'Fibermetric thickness (px)'},'Feuil3','A5');
xlswrite(filename,{'Fibermetric structsensitivity'},'Feuil3','A6');
xlswrite(filename,{'Gaussian blur noise lengthscale lnois (px)'},'Feuil3','A7');
xlswrite(filename,{'Gaussian blur object size lobject (px)'},'Feuil3','A8');
xlswrite(filename,{'Gaussian blur threshold'},'Feuil3','A9');
xlswrite(filename,{'Binarization sensitivity'},'Feuil3','A10');
xlswrite(filename,{'Skeletonization MinBranchLength (px)'},'Feuil3','A11');
xlswrite(filename,{'Bspline ds'},'Feuil3','A12');
xlswrite(filename,{'Bspline npnts'},'Feuil3','A13');
%
% Writing out data by rows on the Excel sheet 2
writematrix(basepath,filename,'Sheet',3,'Range','B2');
writematrix(tifname,filename,'Sheet',3,'Range','B3');
writematrix(FilNum,filename,'Sheet',3,'Range','B4');
writematrix(thickness,filename,'Sheet',3,'Range','B5');
writematrix(structsensitivity,filename,'Sheet',3,'Range','B6');
writematrix(lnoise,filename,'Sheet',3,'Range','B7');
writematrix(lobject,filename,'Sheet',3,'Range','B8');
writematrix(threshold,filename,'Sheet',3,'Range','B9');
writematrix(sensitivity,filename,'Sheet',3,'Range','B10');
writematrix(MinBranchLength,filename,'Sheet',3,'Range','B11');
writematrix(ds,filename,'Sheet',3,'Range','B12');
writematrix(npnts,filename,'Sheet',3,'Range','B13');
%
% GIVING NAMES TO THE EXCEL SHEETS
e = actxserver('Excel.Application'); % # open Activex server
ewb = e.Workbooks.Open('C:\Users\Faustine\Documents\POSTDOC\Image treatment\Francesco - Matlab\Modified code\additionaldata.xlsx'); % # open file (enter full path!)
ewb.Worksheets.Item(1).Name = 'Data'; % # rename 1st sheet
ewb.Worksheets.Item(2).Name = 'Frames and Flagellum'; % # rename 2nd sheet
ewb.Worksheets.Item(3).Name = 'Code Parameters'; % # rename 3rd sheet
ewb.Worksheets.Item(2).Range('A1:B1').Interior.Color=hex2dec('F0F4C3'); % # color row A1 in sheet 2
ewb.Worksheets.Item(2).Range('A7:B7').Interior.Color=hex2dec('F0F4C3'); % # color row A6 in sheet 2
ewb.Worksheets.Item(3).Range('A1:B1').Interior.Color=hex2dec('F0F4C3'); % # color row A1 in sheet 3
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