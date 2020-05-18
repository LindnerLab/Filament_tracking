%% FAUSTINE'S FILE
%% ESSENTIAL DATA ARE SAVED TO AN EXCEL FILE
% This file is a short (and essential) version of additionaldata.m
% It can only be run after having run additionaldata.m

%~~~~ VARIABLES AND PARAMETERS THAT ARE RECORDED ARE:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% * FRAMES AND FLAGELLUM
% ** Missingpercentage (%) = length(xy.emptyframe)/total (percentage of missing frames on total chosen number of frames)
% ** Diameter (µm) = initial filament diameter (measured by hand)
% ** Chosenlength (µm) = fil_length (at the end of the additionaldata.m process)
% ** Lav5MIC (µm) = Average maximal length on 5% maximal Lp values
%
% * DATA
% ** gammadot (s-1) = shear rate
% ** Decorrthreshold (s) = 5% IN A ROW of the autocorr values of total number of frames included in [-0.005;0.005]); sufficient to say that decorrelation has occurred.
% ** Maxdecorrduration (s) = maximal time when autocorr values remain IN A ROW in [-0.005;0.005].
% ** Decorrelation (Y/N) = decorrelation has occurred or not (based on Decorrthreshold and Maxdecorrduration)
% ** Limitcorrtime (s) = starting decorrelation time (if one was found in agreement with Decorrthreshold and Maxdecorrduration)
% ** Limitcorrtimesmooth (s) = starting decorrelation time; corresponds to tau (calculated from expofit.b, the exponential fit of autocorr(Cm) which is very sensitive to noise (bias)
% ** tau_r (s) = rotational diffusion time
% ** tJ (s) = Jeffery oscillation period (in ideal case with no noise)

%% CODE
% Name of the file when created (created in the code folder by default)
filename = 'essentialdata.xlsx';

% VARIABLES & PARAMETERS
Missingpercentage = length(xy.emptyframe)/total*100;
Diameter = in_diameter;
Chosen_length = fil_length;
Fivepercent_averagelength = Lav5MIC;
Shear_rate = gammadot;
Decorrthreshold = fivepercent * FSav;
Maxdecorrduration = Sum * FSav;
if Maxdecorrduration >= Decorrthreshold
    Decorrelation = 'YES';
else
    Decorrelation = 'NO';
end
Limitcorrtime = Limitcorrtime;
Limitcorrtimesmooth = tau;
tau_r = tau_r;
tJ = tJ;

% WRITING OUT NAMES
xlswrite(filename,{'FRAMES AND FLAGELLUM'},'Feuil1','A1');
xlswrite(filename,{'Missing percentage (%)'},'Feuil1','A2');
xlswrite(filename,{'Diameter (µm)'},'Feuil1','A3');
xlswrite(filename,{'Chosen length (µm)'},'Feuil1','A4');
xlswrite(filename,{'Max 5% average length (µm)'},'Feuil1','A5');
%
xlswrite(filename,{'DATA'},'Feuil1','A6');
xlswrite(filename,{'Decorrelation'},'Feuil1','A7');
xlswrite(filename,{'Shear rate (s-1)'},'Feuil1','A8');
xlswrite(filename,{'Decorrelation threshold (s)'},'Feuil1','A9');
xlswrite(filename,{'Maximal decorrelation duration vs. threshold (%)'},'Feuil1','A10');
xlswrite(filename,{'Decorrelation time (s)'},'Feuil1','A11');
xlswrite(filename,{'Smooth decorrelation time tau (s)'},'Feuil1','A12');
xlswrite(filename,{'tau_r (s)'},'Feuil1','A13');
xlswrite(filename,{'tJ (s)'},'Feuil1','A14');

% WRITING OUT VALUES
writematrix(Missingpercentage,filename,'Sheet',1,'Range','B2');
writematrix(Diameter,filename,'Sheet',1,'Range','B3');
writematrix(Chosen_length,filename,'Sheet',1,'Range','B4');
writematrix(Fivepercent_averagelength,filename,'Sheet',1,'Range','B5');
%
writematrix(Decorrelation,filename,'Sheet',1,'Range','B7');
writematrix(Shear_rate,filename,'Sheet',1,'Range','B8');
writematrix(Decorrthreshold,filename,'Sheet',1,'Range','B9');
writematrix(Maxdecorrduration,filename,'Sheet',1,'Range','B10');
writematrix(Limitcorrtime,filename,'Sheet',1,'Range','B11');
writematrix(Limitcorrtimesmooth,filename,'Sheet',1,'Range','B12');
writematrix(tau_r,filename,'Sheet',1,'Range','B13');
writematrix(tJ,filename,'Sheet',1,'Range','B14');

% GIVING NAMES TO THE EXCEL SHEETS
e = actxserver('Excel.Application'); % # open Activex server
ewb = e.Workbooks.Open('C:\Users\Faustine\Documents\POSTDOC\Image treatment\Francesco - Matlab\Modified code\essentialdata.xlsx'); % # open file (enter full path!)
hWorksheet = ewb.Worksheets.Item(1);
hWorksheet.Range('A1:B1').Interior.Color=hex2dec('F0F4C3'); % # color row A1 - FRAMES & FLAGELLUM
hWorksheet.Range('A6:B6').Interior.Color=hex2dec('F0F4C3'); % # color row A6 - DATA
hWorksheet.Columns.Item(1).columnWidth = 50; % width of the first column
ewb.Save;
ewb.Close(false)
e.Quit