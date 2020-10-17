clc;clearvars;close all;
% Set parameters
xl = [-45 45];          % x axis limit
yl = [-45 45];          % y axis limit
start = [0,0];          % starting coordinate (x,y)
gain = 0.01;            % amount of motion at each iteration (% of axis range)
time = 30;              % run time (seconds)
showTrace = false;      % when true, the path is recorded
% initialize figure
figure
ah = axes; 
h = plot(start(1), start(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); 
xlim(ah, xl)
ylim(ah, yl)
hold(ah, 'on')
% randomly jitter for set amount of time
t1 = tic; 
pt = start; 
stepX = range(xl) * gain; 
stepY = range(yl) * gain; 
while toc(t1) < time
   % determine range of possible next coordinates
   xbound = [max(pt(1)-stepX, xl(1)), min(pt(1)+stepX, xl(2))]; 
   ybound = [max(pt(2)-stepY, yl(1)), min(pt(2)+stepY, yl(2))];
   % Randomly choose next coordinate
   pt(1) = xbound((rand < 0.5) + 1); 
   pt(2) = ybound((rand < 0.5) + 1); 
   % plot segment
   if showTrace
       plot([h.XData, pt(1)], [h.YData, pt(2)], '-k')
   end
   % Update position of dot on axis
   h.XData = pt(1); 
   h.YData = pt(2); 
   drawnow
end
