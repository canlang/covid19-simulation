clearvars;clc;
% Setup parameters
num_par = 500;            % number of particles

xl = [0 100];           
yl = [0 100];           
start = [randi(xl,num_par,1),randi(yl,num_par,1)];  % starting coordinate (x,y)
gain = 0.01;              % movement 
time = 30;        
day = 1;
limit_day = 100;

% Epi parameter
r_infect = 5;
infectP = 0.8;
i_period = 5;

%% % Initialize figure
close all;

figure
set(gcf,'units','points','position',[500,500,1200,400])
ax1 = subplot(121);
h = scatter(start(:,1), start(:,2),35,...
    repmat([0 1 0],num_par,1),'filled','MarkerFaceAlpha',.5);
xlim(ax1, xl)
ylim(ax1, yl)
% hold(ax1, 'on')


ax2 = subplot(122);
h2 = animatedline(0,0,'Color','m','LineWidth',1,'Marker','o');
hold on
h3 = animatedline(0,0,'Color','k','LineWidth',1);
% ylim(ax2, [0 num_par])
% xlim(ax2, [0 limit_day])
xlabel('Time (day)');
ylabel('# of infect');
legend('Today','Current');

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

video_flag = 0;
switch video_flag
    case 1
        fprintf('Saving vids...\n');
        video_filename = sprintf('pt%d-idot%.1f-period%d',num_par,infectP,i_period);
        v = VideoWriter(strcat('vids/',video_filename),'MPEG-4');
        v.FrameRate = 10;
        v.Quality = 100;
        open(v);
        frame = getframe(gcf);
        writeVideo(v,frame);
end
%%
while toc(t1) < time && any(pt(:,3)>=1&pt(:,3)<=i_period)            
   pt(pt(:,3)>0,3) = pt(pt(:,3)>0,3)+1;
   h.CData(pt(:,3)>i_period,:) = ones(sum(pt(:,3)>i_period),1)*[0 0 1];
   
   % Update movement
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
   % Smooth animation
   dt = 5;
   fun1 = @(a,b) linspace(a,b,dt);
   dx = arrayfun(fun1,pt(:,1),xbound(indx)','UniformOutput',false);
   dy = arrayfun(fun1,pt(:,2),ybound(indy)','UniformOutput',false);
   dx = vertcat(dx{:});
   dy = vertcat(dy{:});
   for i=1:dt
       h.XData = dx(:,i); 
       h.YData = dy(:,i);       
       drawnow 
   end
   
   pt(:,1) = xbound(indx); 
   pt(:,2) = ybound(indy);
   h.XData = pt(:,1); 
   h.YData = pt(:,2);    
    
   % Apply infection transition`
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
   addpoints(h3,day,sum(pt(:,3)>=1&pt(:,3)<=i_period));       
   axis tight
   drawnow

   day=day+1;
      
   switch video_flag
        case 1
        frame = getframe(gcf);
        writeVideo(v,frame);
   end
end
switch video_flag, case 1, close(v);end