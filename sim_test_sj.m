clearvars;clc;
% Set parameters
num_par = 100;            % number of particles

xl = [0 100];           % x axis limit
yl = [0 100];           % y axis limit
start = [randi(xl,num_par,1),randi(yl,num_par,1)];  % starting coordinate (x,y)
gain = 0.01;            % amount of motion at each iteration (% of axis range)
time = 10;              % run time (seconds)
day = 0;
limit_day = 100;

% epi parameter
r_infect = 10;
infectP = 0.8;
i_period = 13;

%%
close all;
% initialize figure
f1 = figure;
set(gcf,'units','points','position',[500,500,600,400])
% h = plot(start(:,1), start(:,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); 
% subplot(2,1,1)
ax = gca;
h = scatter(start(:,1), start(:,2),25,repmat([0 1 0],num_par,1),'filled');
% h = scatter(start(:,1), start(:,2),25,zeros(num_par,1),'filled');

% axis image
xlim(ax, xl)
ylim(ax, yl)
hold(ax, 'on')
% grid on
%%
% randomly jitter for set amount of time
t1 = tic; 
pt = [start,zeros(num_par,1)];
% first infection
pt(1:2,3) = 1;
h.CData(1:2,:) = ones(2,1)*[1 0 0];

stepX = range(xl) * gain; 
stepY = range(yl) * gain; 
sdf(gcf,'paper_f150')

video_flag = 1;
switch video_flag
    case 1
        video_filename = sprintf('n100-idot8-period%d',i_period);
        v = VideoWriter(strcat('vids/',video_filename),'MPEG-4');
        v.FrameRate = 10;
        v.Quality = 100;
        open(v);
        frame = getframe(f1);
        writeVideo(v,frame);
end

figure
set(gcf,'units','points','position',[1100,500,600,400])
h2 = animatedline;
ax2 = gca;
ylim(ax2, yl)

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
   infectiousIdx = find(pt(:,3)>=1&pt(:,3)<=i_period);
   susceptibleIdx = find(pt(:,3)==0);
   sus_loc = pt(susceptibleIdx,1:2);
   inf_loc = pt(infectiousIdx,1:2);
   [i_idx,i_dist] = knnsearch(inf_loc,sus_loc,'K',1);
   if any(i_dist<r_infect)           % get index neighbor in distance r m
       rndSel = rand(sum(i_dist<r_infect),1)<infectP;
       contactIdx = find(i_dist<r_infect);
       infectionIdx = susceptibleIdx(contactIdx(rndSel)); 
       if ~isempty(infectionIdx)           
           pt(infectionIdx,3) = 1;
           h.CData(infectionIdx,:) = repmat([1 0 0],length(infectionIdx),1);           
       end
       day_infect=sum(rndSel);
       addpoints(h2,day,day_infect);
   else
       addpoints(h2,day,0);
   end
   
   if sum(pt(:,3)>0)
       pt(pt(:,3)>0,3) = pt(pt(:,3)>0,3)+1;
       h.CData(pt(:,3)>i_period,:) = ones(sum(pt(:,3)>i_period),1)*[0 0 1];
   end
   
      
   % Update position of dot on axis
   h.XData = pt(:,1); 
   h.YData = pt(:,2);
   
   drawnow
   day=day+1;
   
   switch video_flag
        case 1
        frame = getframe(f1);
        writeVideo(v,frame);
   end
end
switch video_flag, case 1, close(v);end