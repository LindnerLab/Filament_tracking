%% This is Faustine's file
filename = 'additionaldata.xlsx';
%% calculate the angle

for j = 1 : length(xy.frame)
angle(j) = (xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))/(xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1));
angleindeg(j) = angle(j) * 180 / pi;
end
writematrix(angleindeg,filename,'Sheet',1);

%% calculate the length in microns based on coordinates

for j = 1 : length(xy.frame)
    L(j) = sqrt((xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1))^2+(xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))^2);
    Linmicrons(j) = L(j)/10.24; %10.24 is the conversion factor pixel/degrees for imgs of 1024x1024 of 100x100microns.
end
writematrix(Linmicrons,filename,'Sheet',2);

%% print a 3D plot with the centroid coordinates (x,y) as a function of the frame number
% trajectory in microns (1024x1024 px -> 100x100 microns)

for j = 1 : length(xy.frame)
    CENTROID_X(j) = xy.centroid{j}(1,1)*100/1024;
    CENTROID_Y(j) = xy.centroid{j}(1,2)*100/1024;
end
plot3(xy.frame,CENTROID_X,CENTROID_Y);
%axis equal
xlabel('Frame nb')
ylabel('Centroid X')
zlabel('Centroid Y')