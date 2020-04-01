filename = 'additionaldata.xlsx';
%% calculate the angle (Faustine)

for j = 1 : length(xy.frame)
angle(j) = (xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))/(xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1));
end
writematrix(angle,filename,'Sheet',1);

%% calculate the length in microns based on coordinates (Faustine)

for j = 1 : length(xy.frame)
    L(j) = sqrt((xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1))^2+(xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))^2);
    Linmicrons(j) = L(j)/10.24;
end
writematrix(Linmicrons,filename,'Sheet',2);