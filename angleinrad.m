%% calculate the angle (Faustine)

for j = 1 : length(xy.frame)
angle(j) = (xy.spl{j}(length(xy.spl{j}),2)-xy.spl{j}(1,2))/(xy.spl{j}(length(xy.spl{j}),1)-xy.spl{j}(1,1));
end
filename = 'angle.xlsx';
writematrix(angle,filename,'Sheet',1);