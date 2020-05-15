%% FAUSTINE'S FILE
%% PRODUCE PLOTS & AN EXCEL FILE WITH DATA NECESSARY FOR COMPARING JEFFERY VS. BROWNIAN MOTION

%~~~~ PARAMETERS THAT ARE RECORDED ARE:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% * FRAME PARAMETERS:
%
% ** in_frame = First frame (real frame number)
% ** fin_frame = Last frame (real frame number)
% ** step = frame step (a frame is analyzed every ... frames)
% ** total = chosen number of frames (missing frames included)
% ** F Frequence ("frame per second") in Hz
% ** FS frame step in seconds (time elapsed between two frames)
% ** FSav average frame step in seconds (missing frames excluded)
% ** Frameindex (number of treated frame, missing frame excluded)
% ** Framenumber (number of treated frame, missing frames included)
% ** Time = frame index * average frame step (in seconds)
%
% * FLAGELLUM PARAMETERS
%
% ** diameter (equivalent of Jeffery width i.e. between the top of one side of the helix to the other side of the helix)
% ** fil_length (length of a straight line going from one end of the filament to the other)
% ** lambda = aspect ratio

%~~~~ VARIABLES THAT ARE CALCULATED ARE:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% * IN THE XY PLANE (observation plane):
%
% ** nx and ny (nxMIC and nyMIC) : delta x and delta y on the endpoints of the detected filament;
% ** Lp (LpMIC): projected lengths (in pixels and in microns) based on coordinates
% ** arclen (arclenMIC): projected arc length (in pixels and in microns)
% ** Xsi (xsiindeg): projected angle = atan(delta y / delta x) in radians and degrees
% ** REAL LENGTH is estimated via LmeanMIC (mean of LpMIC) and Lav5MIC;
%    Lav5MIC corresponds to the average of the 5% maximal values taken by Lp over the
%    course of all frames treated. 5%maximal values are calculated with
%    the maxk(A,k) function taken the k maximal values of the array A.
%    round() approximates to the nearest integer.
%    --> The length Lav5MIC replaces fil_length when calculated, if
%    superior to the initial fil_length calculated by hand.
%
% * OUT-OF-PLANE (IN 3D)
% 
% ** theta (thetaindeg): angle between the filament and the y axis
%
% * IN THE XZ PLANE (out-of-plane):
%
% ** Lpz (LpzMIC) = projected lengths (in pixels and in microns) such as Lpz = cos(phi)*Lav5MIC
% ** nz (nzMIC) such as nzMIC = sin(theta) * LpzMIC
% ** phi (phiindeg): angle between Lpz and the x axis
%
% * INDICATORS TO COMPARE JEFFERY VS. BROWNIAN MOTION
% ** Jeffery constant C (depends on nx, ny, nz, and the filament aspect ratio)
% ** Modified Jeffery constant Cm = sign(C/(1 + abs(C)) ou plus simplement C/(1+C)
% ** Autocorrelation function of Cm:
%       autocorr(function, 'Numlags', number of lags) --> gives a Jeffery decay time tau
%       * Y coordinates stored within a variable
%       * X coordinates are 1,2,3... as numerous as the number of lags since the function is discretized
%       * For the fit, Xcorr and Ycorr must be column vectors --> T_Xcorr and T_Ycorr
%       * Exponential fit exp1 = a * exp(b * k) where k is a nb of frames such as t = k*FSav (t=time)
%       * Tau (time constant) such as exp(-t/tau) -> tau = -FSav/b
% ** Rotational diffusion time tau_r
% ** Jeffery oscillation period tJ (depends on the filament aspect ratio and the shear rate)
%       * gammadot is the shear rate (in s-1); 18 s-1 is the value given by Zöttl et al
%       * tJ is the Jeffery oscillation period (in s) according to Zöttl et al., 2019
%       * Losingmemory is the number of Jeffery oscillations a filament performs before losing memory about its Jeffery orbit state

%% CODE
% Name of the file when created (created in the code folder by default)
filename = 'additionaldata.xlsx';

% FRAME INFO
in_frame = xy.frame(1);
fin_frame = xy.frame(xy.nframe);
total = (xy.nframe + length(xy.emptyframe));
step = frame_step;
F = 30;
FS = 1/F*step;
FSav = (total / xy.nframe)*FS;
Frameindex = transpose([1:length(xy.frame)]);
Framenumber = transpose(xy.frame);
Time = transpose([1:length(xy.frame)]*FSav);

% FLAGELLUM INFO
% Batch 1: diameter = 0.879; fil_length=8.301;
% Batch 2: diameter = 1.17; fil_length = 8.35; mesured on 22971: d=0.912; l=8.076;
% Batch 3: diameter = 0.69; fil_length = 8.77;
% Batch 4: diameter = 0.869; fil_length = 10.215;
% Batch 5: diameter = 0.977; fil_length = 7.887;
% Batch 6: diameter = 0.834; fil_length = 4.844;
% Batch 7: diameter = 0.886; fil_length = 7.852;
in_diameter = 0.912;
in_fil_length = 8.076;
fil_length = in_fil_length;
lambda = fil_length/diameter;

% ON THE XY PLANE
% COORDINATES & LENGTHS (px and microns) AND ANGLE
for j = 1 : xy.nframe
    nx(j) = xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1);
    ny(j) = xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2);
    Lp(j) = sqrt(nx(j)^2 + ny(j)^2);
    nxMIC(j) = nx(j) * 100 / 1024;
    nyMIC(j) = ny(j) * 100 / 1024;
    LpMIC(j) = Lp(j) * 100 / 1024;
    arclenMic(j) = xy.arclen(j)*100/1024;
    xsi(j) = ny(j)/nx(j);
    xsiindeg(j) = xsi(j) * 180 / pi;
end
% REAL LENGTH
fivepercent = round((5/100) * xy.nframe);
L5 = maxk(transpose(Lp),fivepercent);
Lav5 = sum(L5)/fivepercent;
Lav5MIC = Lav5 * 100 / 1024;
%
Lmean = sum(Lp)/xy.nframe;
LmeanMIC = sum(LpMIC)/xy.nframe;

% NEW CALCULATION OF FILAMENT LENGTH AND ASPECT RATIO BASED ON THE 5% MAXIMAL VALUES OF Lp
if Lav5MIC > fil_length
    fil_length = Lav5MIC;
    lambda = fil_length/diameter;
    disp('Initial filament length was replaced by Lav5MIC.')
else
    fil_length = fil_length;
    lambda = fil_length/diameter;
    disp('Initial filament length was conserved.')
end

% OUT-OF-PLANE (IN 3D)
for j = 1 : xy.nframe
    theta(j) = acos(nyMIC(j)/fil_length);
    thetaindeg(j) = theta(j) * 180 / pi;
end

% ON THE XZ PLANE
% COORDINATE AND PROJECTED LENGTH (px and microns)
for j = 1 : xy.nframe
        % PROJECTED LENGTH
    Lpz(j) = (fil_length/100*1024) * sin(theta(j));
    LpzMIC(j) = fil_length * sin(theta(j));
        % ANGLE
    phi(j) = acos(nx(j)/Lpz(j));
    phi(j) = real(phi(j)); % removing imaginary phi in case Lp > fil_length (under- or over-estimate of length)
    phiindeg(j) = phi(j) * 180 / pi;
        % COORDINATE
    nz(j) = sin(phi(j)) * Lpz(j);
    nzMIC(j) = sin(phi(j)) * LpzMIC(j);
end

% NEW UNIT COORDINATES AND ANDREAS COORDINATES (UNY, UNZ)
for j = 1 : xy.nframe
    Unx(j) = nxMIC(j)/fil_length;
    Uny(j) = nyMIC(j)/fil_length;
    Unz(j) = nzMIC(j)/fil_length;
        % Andreas' coordinates
    UNY(j) = Unz(j);
    UNZ(j) = Uny(j);
end
    
% JEFFERY CONSTANT C AND MODIFIED CONSTANT Cm
for j = 1 : xy.nframe
    CAndreas(j) = sqrt(Unx(j)^2 + (UNZ(j)^2/lambda^2))/UNY(j);
    Cm(j) = sign(CAndreas(j))/(1+abs(CAndreas(j)));
    Cprime(j) = CAndreas(j)/(1+CAndreas(j));
        % IF HORIZONTAL HELE-SHAW CELL
    %Chorizontal(j) = sqrt(Unx(j)^2 + (UNY(j)^2/lambda^2))/Uny(j);
    %Cm(j) = sign(Chorizontal(j))/(1+abs(Chorizontal(j)));
end

% AUTOCORRELATION OF Cm AND Tau
NbofLags = xy.nframe-1;
Ycorr = autocorr(Cm,'Numlags',NbofLags);
Xcorr = (0:xy.nframe-1)*FSav;
T_Ycorr = transpose(Ycorr);
T_Xcorr = transpose(Xcorr);
% CORRELATION TIME
for j = 1 : xy.nframe
    if -0.005 < Ycorr(j) && Ycorr(j) < 0.005
        Corrframe(j) = j;
    else
        Corrframe(j) = 0;
    end
    CORR(j) = Corrframe(j) > 0; %Logical matrix: 1 = Ycorr ~ 0
    Indexones = Corrframe(CORR); %Records the frame indexes when Ycorr ~ 0
    Sum = 0;
    for i = 1 : length(Indexones)-1
        if Indexones(i+1) == Indexones(i) + 1 %If there is no frame left between two frame indexes where Ycorr ~ 0, then start incrementing 
            Sum = 1 + Sum;
        else
            Sum = Sum;
        end
        if Sum == fivepercent-1 %When Ycorr has been ~ 0 for a while (5% of total number of frames), decorrelation is supposed to have happened.
            Limitcorrframe = Indexeones(i-fivepercent); %Stores the frame index at which Ycorr starts being 0 for a long time
        elseif Sum > fivepercent+1
            Difference = Sum - fivepercent+1;
            Limitcorrframe = Indexones(i-fivepercent-Difference);
        else
            Limitcorrframe = 0;
        end
    end
end
Maximaldurationwhenautocorrisnul = (Sum+1)*FSav;
if Limitcorrframe == 0
    Limitcorrtime = 0;
    input_plotautocorr = input('Decorrelation did not occur. Do you still want to plot the autocorr function? Press: \n 1 = yes \n 2 = no \n');
    if input_plotautocorr == 1
        input_NbofLags = input('Choose your decorrelation frame : \n');
        Ycorr = autocorr(Cm,'Numlags',input_NbofLags);
        Xcorr = (1:input_NbofLags+1)*FSav;
        T_Ycorr = transpose(Ycorr);
        T_Xcorr = transpose(Xcorr);
        expofit = fit(T_Xcorr,T_Ycorr,'exp1','StartPoint',[T_Xcorr(100),T_Ycorr(100)]); %Limitcorrtime with smoothing of the curve
        fita = expofit.a;
        tau = (-FSav/expofit.b);
    else
        input_NbofLags = 0;
        disp('No autocorr function plotted. Tau is not valid.')
        tau = 0;
    end
else
    Limitcorrtime = Limitcorrframe * FSav; %equivalent to tau (without smoothing)
    input_NbofLags = input(strcat('The decorrelation frame is ', num2str(Limitcorrframe),' over ', num2str(xy.nframe),' frames. Choose your decorrelation frame : \n'));
    Ycorr = autocorr(Cm,'Numlags',input_NbofLags);
    Xcorr = (1:input_NbofLags+1)*FSav;
    T_Ycorr = transpose(Ycorr);
    T_Xcorr = transpose(Xcorr);
    expofit = fit(T_Xcorr,T_Ycorr,'exp1','StartPoint',[T_Xcorr(100),T_Ycorr(100)]); %Limitcorrtime with smoothing of the curve
    fita = expofit.a;
    tau = (-FSav/expofit.b);
end

% ROTATIONAL DIFFUSION TIME FOR PROLATE BODIES
%       * Boltzmann constant kb in m2 kg s-2 K-1
%       * T Temperature in Kelvin (ambient temperature)
%       * eta dynamic viscosity in Pa.s (dynamic viscosity of water used as a first pass)
kb = 1.38064852 * 10^(-23);
T = 20 + 273.15;
eta = 10^-3;
%       * a half-length of the filament in m
%       * b half-width of the filament in m
%       * p new aspect ratio needed for Nuris' diffusion coefficient formula (cf. her thesis p.47).
%       * V Volume of a prolate ellipsoid in m3 (Nuris' thesis)
%       * g and S coefficients adapted to prolate ellipsoids (Nuris'formula eq. 11 and 12)
%       * Dr diffusion coefficient for a prolate ellipsoid; Nuris' formula (eq. 10)
a =(fil_length/2)*10^(-6); 
b = (diameter/2)*10^-6; 
p = a/b;
V = (4*pi*a*(b^2))/3; 
S = (1 / sqrt(p^2-1)) * log(p+sqrt(p^2-1));
g = (2*(p^4-1)) / (3 * p * ((2*p^2-1)*S - p));
Dr = (kb * T) / (6*eta*V*g);
%       * ROTATIONAL DIFFUSION TIME
tau_r = 1/(2*Dr);

% JEFFERY OSCILLATION PERIOD tJ
% Batch 1: gammadot = 6.3 using shear_y[250,37] because the filament was determined to be 37microns away from the closest channel wall using ImageJ, and apprx. 250microns away from top and bottom.
% Batch 2: gammadot = 16.8 using shear_y[250,20]
% Batch 3: gammadot = 8.797 using shear_y[250,33]
% Batch 4: gammadot = 7.58 using shear_y[250,35]
% Batch 5: gammadot = 13.04 using shear_y[250,26]
% Batch 6: gammadot = 12.43 using shear_Y[250,27]
% Batch 7: inconnu 
% !!!! in some cases, gammadot changes with time as the filament deviates from a straight line trajectory
gammadot = 16.8;
tJ = (2*pi*(lambda + 1/lambda))/gammadot;
Losingmemory = tau/tJ;

%~~~~PLOTTING FIGURES
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

close all

% PLOTTING FUNCTIONS OF TIME
% ** (1) Phi and CAndreas (which is C)
% ** (2) Phi and Cm
% ** (3) Phi and LpMIC
% ** (4) CAndreas and Cm
% ** (5) Lp and Unx, Uny, Unz
%
% ** Figure (1): Phi and CAndreas
%
figure()
title('Phi(t) and C(t)')
%Plotting on the same figure
hold on
xlabel('Time (s)')
%Left vertical axis
yyaxis left
ylabel('Phi (deg)')
plot(Time, phiindeg);
% Right vertical axis
yyaxis right
ylabel('CAndreas');
plot(Time, CAndreas);
%End of plotting on the same figure
hold off
saveas(gcf,'Fig1_Phi-and-C','pdf');
saveas(gcf,'Fig1_Phi-and-C','fig');
%
% ** Figure (2): Phi(t) and Cm(t)
%
figure()
title('Phi(t) and Cm(t)')
hold on
xlabel('Time (s)')
yyaxis left
ylabel('Phi (deg)')
plot(Time, phiindeg);
yyaxis right
ylabel('Cm');
plot(Time, Cm);
hold off
saveas(gcf,'Fig2_Phi-and-Cm','pdf');
saveas(gcf,'Fig2_Phi-and-Cm','fig');
%
% ** Figure (3): Phi and LpMIC
%
figure()
title('Phi(t) and Lp(t)')
hold on
xlabel('Time (s)')
yyaxis left
ylabel('Phi')
plot(Time, phiindeg);
yyaxis right
ylabel('Lp (µm)');
plot(Time, LpMIC);
hold off
saveas(gcf,'Fig3_Phi-and-LpMIC','pdf');
saveas(gcf,'Fig3_Phi-and-LpMIC','fig');
%
% ** Figure (4): C and Cm
%
figure()
title('C(t) and Cm(t)')
hold on
xlabel('Time (s)')
yyaxis left
ylabel('C')
plot(Time, CAndreas);
yyaxis right
ylabel('Cm');
plot(Time, Cm);
hold off
saveas(gcf,'Fig4_C-and-Cm','pdf');
saveas(gcf,'Fig4_C-and-Cm','fig');
%
% ** Figure (5): LpMIC, Unx, Uny, UNZ
figure()
title('Lp(t), Unx(t), Uny(t), and Unz(t)')
hold on
xlabel('Time (s)')
yyaxis left
ylabel('Lp (µm)')
plot(Time, LpMIC);
yyaxis right
ylabel('Unx, Uny, Unz');
plot(Time, Unx);
plot(Time, Uny);
plot(Time, Unz);
hold off
legend('Lp','Unx','Uny','Unz')
saveas(gcf,'Fig5_LpMIC-Unx-Uny-Unz','pdf');
saveas(gcf,'Fig5_LpMIC-Unx-Uny-Unz','fig');

% PLOTTING THE AUTOCORRELATION FUNCTION
if input_NbofLags == 0
    disp('Decorrelation did not occur.')
else
    figure()
    hold on
    if ne(Limitcorrtime,0) == 1 %If Limitcorrtime is not equal to 0
        xline(Limitcorrtime,'r:','Decorrelation time');
    end
    plot(expofit,T_Xcorr,T_Ycorr);
    hold off
    saveas(gcf,'Fig6_autocorr-Cm','pdf');
    saveas(gcf,'Fig6_autocorr-Cm','fig');
end

% PLOT PROBABILITY DENSITY FUNCTIONS (PDF) FOR Lp, PHI, THETA, Cm
figure()
PDF_Lp(LpMIC)
saveas(gcf,'Fig7_PDF-Lp','pdf');
saveas(gcf,'Fig7_PDF-Lp','fig');
figure()
PDF_phi(phi)
saveas(gcf,'Fig8_PDF-Phi','pdf');
saveas(gcf,'Fig8_PDF-Phi','fig');
figure()
PDF_theta(theta)
saveas(gcf,'Fig9_PDF-Theta','pdf');
saveas(gcf,'Fig9_PDF-Theta','fig');
figure()
PDF_xsi(xsi)
saveas(gcf,'Fig10_PDF-Xsi','pdf');
saveas(gcf,'Fig10_PDF-Xsi','fig');
figure()
PDF_Cm(Cm)
saveas(gcf,'Fig11_PDF-Cm','pdf');
saveas(gcf,'Fig11_PDF-Cm','fig');
%D'après Martyna:
%figure;
%histogram(phi,nbdebarres,'Normalization','pdf')

%~~~~WRITING OUT DATA IN AN EXCEL FILE
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% WRITING OUT THE VARIABLES OF INTEREST IN AN EXCEL FILE
% * DATA (SHEET 1)
% 
% ** Column titles for the Data sheet
xlswrite(filename,{'Frame index'},'Feuil1','A1');
xlswrite(filename,{'Frame n°'},'Feuil1','B1');
xlswrite(filename,{'Time (s)'},'Feuil1','C1');
xlswrite(filename,{'Unx'},'Feuil1','D1');
xlswrite(filename,{'Uny'},'Feuil1','E1');
xlswrite(filename,{'Unz'},'Feuil1','F1');
xlswrite(filename,{'Lp (µm)'},'Feuil1','G1');
xlswrite(filename,{'Lpz (µm)'},'Feuil1','H1');
xlswrite(filename,{'Phi (deg)'},'Feuil1','I1');
xlswrite(filename,{'Theta (deg)'},'Feuil1','J1');
xlswrite(filename,{'Andreas Jeff. C'},'Feuil1','K1');
xlswrite(filename,{'Modif. Jeff. Cm'},'Feuil1','L1');
xlswrite(filename,{'Maximal duration when autocorr=0 (s)'},'Feuil1','M1');
xlswrite(filename,{'Decorrelation time (s)'},'Feuil1','N1');
xlswrite(filename,{'Expofit coeff tau (s)'},'Feuil1','O1');
xlswrite(filename,{'Rot. diff. time tau_r (s)'},'Feuil1','P1');
xlswrite(filename,{'Jeff. period tJ (s)'},'Feuil1','Q1');
xlswrite(filename,{'Losing memory ratio tau/tJ'},'Feuil1','R1');
% FOR HORIZONTAL HELE-SHAW CELLS
%xlswrite(filename,{'Horiz. Jeff. C'},'Feuil1','K1');
%
% ** Writing out data by columns on the Excel sheet 1
% Using transpose() because data are by default organized in rows instead of columns
writematrix(Frameindex,filename,'Sheet',1, 'Range', 'A2');
writematrix(Framenumber,filename,'Sheet',1, 'Range', 'B2');
writematrix(Time,filename,'Sheet',1, 'Range', 'C2');
writematrix(transpose(Unx),filename,'Sheet',1, 'Range', 'D2');
writematrix(transpose(Uny),filename,'Sheet',1, 'Range', 'E2');
writematrix(transpose(Unz),filename,'Sheet',1, 'Range', 'F2');
writematrix(transpose(LpMIC),filename,'Sheet',1, 'Range', 'G2');
writematrix(transpose(LpzMIC),filename,'Sheet',1, 'Range', 'H2');
writematrix(transpose(phiindeg),filename,'Sheet',1, 'Range', 'I2');
writematrix(transpose(thetaindeg),filename,'Sheet',1, 'Range', 'J2');
writematrix(transpose(CAndreas),filename,'Sheet',1, 'Range', 'K2');
writematrix(transpose(Cm),filename,'Sheet',1, 'Range', 'L2');
writematrix(Maximaldurationwhenautocorrisnul,filename,'Sheet',1, 'Range', 'M2');
writematrix(Limitcorrtime,filename,'Sheet',1,'Range','N2');
writematrix(tau,filename,'Sheet',1,'Range','O2');
writematrix(tau_r,filename,'Sheet',1,'Range','P2');
writematrix(tJ,filename,'Sheet',1,'Range','Q2');
writematrix(Losingmemory,filename,'Sheet',1,'Range','R2');
% FOR HORIZONTAL HELE-SHAW CELLS
%writematrix(transpose(Chorizontal),filename,'Sheet',1, 'Range', 'J2'); %for horizontal Hele-Shaw cells
%
% * FRAMES AND FLAGELLUM (SHEET 2)
%
% ** Row titles for the Frames and Flagellum sheet
xlswrite(filename,{'FRAMES INFO'},'Feuil2','A1');
xlswrite(filename,{'Frequence (Hz)'},'Feuil2','A2');
xlswrite(filename,{'Frame step (dimensionless)'},'Feuil2','A3');
xlswrite(filename,{'Frame step (s)'},'Feuil2','A4');
xlswrite(filename,{'Average frame step (s)'},'Feuil2','A5');
xlswrite(filename,{'First frame treated'},'Feuil2','A6');
xlswrite(filename,{'Final frame treated'},'Feuil2','A7');
xlswrite(filename,{'Chosen number of frames'},'Feuil2','A8');
xlswrite(filename,{'Nb of treated frames'},'Feuil2','A9');
xlswrite(filename,{'Nb of empty frames'},'Feuil2','A10');
%
xlswrite(filename,{'FLAGELLUM INFO'},'Feuil2','A11');
xlswrite(filename,{'Initial diameter (µm)'},'Feuil2','A12');
xlswrite(filename,{'Initial length (µm)'},'Feuil2','A13');
xlswrite(filename,{'Average length Lmean (µm)'},'Feuil2','A14');
xlswrite(filename,{'Average maximal length Lav5 over 5% of Lp (µm)'},'Feuil2','A15');
xlswrite(filename,{'Choosen length (µm)'},'Feuil2','A16');
xlswrite(filename,{'Aspect ratio lambda'},'Feuil2','A17');
%
% Writing out data by rows on the Excel sheet 2
writematrix(F,filename,'Sheet',2,'Range','B2');
writematrix(step,filename,'Sheet',2,'Range','B3');
writematrix(FS,filename,'Sheet',2,'Range','B4');
writematrix(FSav,filename,'Sheet',2,'Range','B5');
writematrix(in_frame,filename,'Sheet',2,'Range','B6');
writematrix(fin_frame,filename,'Sheet',2,'Range','B7');
writematrix(total,filename,'Sheet',2,'Range','B8');
writematrix(xy.nframe,filename,'Sheet',2,'Range','B9');
writematrix(length(xy.emptyframe),filename,'Sheet',2,'Range','B10');
%
writematrix(diameter,filename,'Sheet',2,'Range','B12');
writematrix(in_fil_length,filename,'Sheet',2,'Range','B13');
writematrix(LmeanMIC,filename,'Sheet',2,'Range','B14');
writematrix(Lav5MIC,filename,'Sheet',2,'Range','B15');
writematrix(fil_length,filename,'Sheet',2,'Range','B16');
writematrix(lambda,filename,'Sheet',2,'Range','B17');
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

% GIVING NAMES TO THE EXCEL SHEETS
e = actxserver('Excel.Application'); % # open Activex server
ewb = e.Workbooks.Open('C:\Users\Faustine\Documents\POSTDOC\Image treatment\Francesco - Matlab\Modified code\additionaldata.xlsx'); % # open file (enter full path!)
ewb.Worksheets.Item(1).Name = 'Data'; % # rename 1st sheet
ewb.Worksheets.Item(2).Name = 'Frames and Flagellum'; % # rename 2nd sheet
ewb.Worksheets.Item(3).Name = 'Code Parameters'; % # rename 3rd sheet
ewb.Worksheets.Item(2).Range('A1:B1').Interior.Color=hex2dec('F0F4C3'); % # color row A1 - FRAMES in sheet 2
ewb.Worksheets.Item(2).Range('A11:B11').Interior.Color=hex2dec('F0F4C3'); % # color row A11 - FLAGELLUM in sheet 2
ewb.Worksheets.Item(3).Range('A1:B1').Interior.Color=hex2dec('F0F4C3'); % # color row A1 - CODE in sheet 3
ewb.Save;
ewb.Close(false)
e.Quit