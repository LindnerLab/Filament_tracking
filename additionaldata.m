filename = 'additionaldata.xlsx';
%% calculate the angle (Faustine)

for j = 1 : length(xy.frame)
angle(j) = (xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))/(xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1));
angleindeg(j) = angle(j) * 180 / pi;
end
writematrix(angleindeg,filename,'Sheet',1);

%% calculate the length in microns based on coordinates (Faustine)

for j = 1 : length(xy.frame)
    L(j) = sqrt((xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1))^2+(xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))^2);
    Linmicrons(j) = L(j)/10.24; %10.24 is the conversion factor pixel/degrees for imgs of 1024x1024 of 100x100microns.
end
writematrix(Linmicrons,filename,'Sheet',2);