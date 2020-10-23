function next_loc = getNextLocation(pt, xl, yl, step)

% Update movement
% determine range of possible next coordinates (binary dynamics)
xbound = [max(pt(:,1)-step, xl(1)), min(pt(:,1)+step, xl(2))]; 
ybound = [max(pt(:,2)-step, yl(1)), min(pt(:,2)+step, yl(2))];
% Randomly choose next coordinate
num_par = size(pt,1);
row = (1:num_par)';
colx = randi([1,2],num_par,1);
coly = randi([1,2],num_par,1);
sz = size(xbound);
indx = sub2ind(sz,row,colx);
indy = sub2ind(sz,row,coly);

next_loc = [xbound(indx), ybound(indy)];
end