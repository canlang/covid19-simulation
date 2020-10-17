clc;clearvars;close all;
% Set parameters
num_par = 50;            % number of particles

xl = [0 100];           % x axis limit
yl = [0 100];           % y axis limit
start = [randi(xl,num_par,1),randi(yl,num_par,1)];  % starting coordinate (x,y)
gain = 0.01;            % amount of motion at each iteration (% of axis range)
time = 30;              % run time (seconds)
day = 0;
limit_day = 100;
showTrace = true;      % when true, the path is recorded
% epi parameter
infectP = 0.3;
%%
% initialize figure

% figure
% subplot(2,1,1)
ah = axes; 
% h = plot(start(:,1), start(:,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); 
h = scatter(start(:,1), start(:,2),50,repmat([0 1 0],num_par,1),'filled');

% xlim(xl);
% xlim(yl);
xlim(ah, xl)
ylim(ah, yl)
% axis equal
hold(ah, 'on')
% grid on
% hold on
%%
% randomly jitter for set amount of time
t1 = tic; 
pt = [start,zeros(num_par,1)];
pt(1:2,3) = 1;
h.CData(1:2,:) = repmat([1 0 0],2,1);

stepX = range(xl) * gain; 
stepY = range(yl) * gain; 

%%
while toc(t1) < time && ~all(pt(:,3))
% while day < limit_day
   % determine range of possible next coordinates (binary dynamics)
   xbound = [max(pt(:,1)-stepX, xl(1)), min(pt(:,1)+stepX, xl(2))]; 
   ybound = [max(pt(:,2)-stepY, yl(1)), min(pt(:,2)+stepY, yl(2))];
   % Randomly choose next coordinate
   row = 1:num_par;
   colx = randi([1,2],1,num_par);
   coly = randi([1,2],1,num_par);
   sz = size(xbound);
   indx = sub2ind(sz,row,colx);
   indy = sub2ind(sz,row,coly);
   pt(:,1) = xbound(indx); 
   pt(:,2) = ybound(indy);
   % Apply infection transition
%    pt(:,3) = updateInfection(pt);
   infectedIdx = find(pt(:,3)==1);
   noninfectedIdx = find(pt(:,3)==0);
   sus_loc = pt(noninfectedIdx,1:2);
   inf_loc = pt(infectedIdx,1:2);
   [i_idx,i_dist] = knnsearch(inf_loc,sus_loc,'K',1);
   if sum((i_dist<5) > 0)           % get index neighbor in distance 20 m
       if sum(rand(sum(i_dist<20),1)<infectP)>0
           infectionIdx = noninfectedIdx(rand(sum(i_dist<5),1)<infectP);
           pt(infectionIdx,3) = 1;
           h.CData(infectionIdx,:) = repmat([1 0 0],length(infectionIdx),1);
       end
   end
   
   % plot segment
   if showTrace
       plot(ah,[h.XData(1), pt(1,1)], [h.YData(1), pt(1,2)], '-k')
   end
   % Update position of dot on axis
   h.XData = pt(:,1); 
   h.YData = pt(:,2);
   
   drawnow
   day=day+1;
end
