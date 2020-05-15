%% Movie of a filament under flow &/or brownian motion
%% last modification: 02/04/2020
% allows multiple filaments to be shown in the same frame

%% previous modification: 30/03/2020
% uses directly xy(i).frame(j) & pathintif to find the frame where the filament exists

%% read mat file where coordinates are stored and define path where tiff file is stored
basepath='E:\Helicies in flow-Faustine&Martyna\Martyna\';
path = strcat(basepath,'results\');

disp('select the *mat file')
[file,~] = uigetfile(path);
load(strcat(path,file));
[inext,~] = regexp(file,'.mat');
moviename = file(1:inext-1);

lzero = max(lobject,ceil(5*lnoise)); % size of each edges where gaussian_blur set values to 0


Uframe = 1 : imtot;

%% initialize the video object writer
writerObj = VideoWriter(strcat(path,moviename,'.avi')); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj);  

figure;

for j = 1 : length(Uframe)
    
  pause(0.0001);
    
  img=imread(pathintif,Uframe(j)); % read original frame
  imshow(img(lzero+1:end-lzero,lzero+1:end-lzero),[]); % plot frame
  hold on
   
    for i = 1 : N_fil
        h = ismember(xy(i).frame,Uframe(j));
        kk = find(h==1);        
        if not(isempty(kk))
        % plot reconstructed centerline 
        plot(xy(i).spl{kk}(:,1),size(L,2)-xy(i).spl{kk}(:,2),'color','w','linewidth',0.5);
        hold on
        end
     end
    
xlim([0 size(L,2)])
ylim([0 size(L,1)])

    axis off
    axis equal
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame); 
    hold off % hold off to create a new frame at each call
     
end


hold off
close(writerObj); % Saves the movie.