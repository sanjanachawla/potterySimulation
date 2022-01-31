[X,Y] = meshgrid(linspace(-5, 5, 50));
fcn = @(x,y,k) k*x.^2 + y.^2;
v = [1:-0.05:-1;  -1:0.05:1];
for k1 = 1:2
    for k2 = v(k1,:)
        surfc(X, Y, fcn(X,Y,k2))
        axis([-5  5    -5  5    -30  50])
        drawnow
        
    end
end

%%  trail for rotating elipse?
 clc;
xc=50; %xCenter
yc=50; %yCenter
a=25; %xRad
c=100; %yRad
m = 1000;
x = zeros(m,1);
y = zeros(m,1);
theta = linspace(0,2*pi,m);
for k = 1:m
        x(k) = a * cos(theta(k));
        y(k) = c * sin(theta(k));
end
alpha = input('Enter the rotation angle');
R  = [cos(alpha) -sin(alpha); ...
      sin(alpha)  cos(alpha)];
rCoords = R*[x' ; y'];   
xr = rCoords(1,:)';      
yr = rCoords(2,:)';      
plot(x+xc,y+yc,'r');
grid on;
hold on;
axis equal;
plot(xr+xc,yr+yc,'b');
%% an animated rotating ellipse
step1 = linspace(50,200,100);

figure
hAx = axes('XLim',[-250 250], 'YLim',[-250 250], ...
    'Drawmode','fast', 'NextPlot','add');
axis(hAx, 'equal')

p = calculateEllipse(0, 0, step1(1), step1(end), step1(1));
hLine = line('XData',p(:,1), 'YData',p(:,2), 'EraseMode','xor',  ...
    'Color','r', 'LineWidth',3);

for i=1:numel(step1)
    p = calculateEllipse(0, 0, step1(i), step1(numel(step1)-i+1), step1(i));
    set(hLine,'XData',p(:,1), 'YData',p(:,2))  %# update X/Y data
%     pause(.10)                                 %# slow down animation
    drawnow                                    %# force refresh
    if ~ishandle(hLine), return; end           %# in case you close the figure
end

%% rim (THIS ONE WORKS*****) 

wn = 0.06;
z = 0.1;
tfr = tf([5*wn], [1 2*z*wn wn^2]);
step1 = linspace(40,1200,600);
[r, t] = step(tfr, step1);
%r = r*4;
a = r;
c = 2*dcgain(tfr) - r;
figure
hAx = axes('XLim',[-250 250], 'YLim',[-250 250], ...
    'Drawmode','fast', 'NextPlot','add');
axis(hAx, 'equal')

p = calculateEllipse(0, 0, a(1), c(1), step1(1));
hLine = line('XData',p(:,1), 'YData',p(:,2), 'EraseMode','xor',  ...
    'Color','r', 'LineWidth',3);

for i=1:numel(step1)
    p = calculateEllipse(0, 0, a(i), c(i), step1(i));
    set(hLine,'XData',p(:,1), 'YData',p(:,2))  %# update X/Y data
%     pause(.10)                                 %# slow down animation
    drawnow                                    %# force refresh
    if ~ishandle(hLine), return; end           %# in case you close the figure
end

%% Now to make it 3d!!!!

wn = 0.06;
z = 0.1;
tfr = tf([5*wn], [1 2*z*wn wn^2]);

step1 = linspace(20,1200,600);
[r, t] = step(tfr, step1);
%r = r*4;
a = r;
c = 2*dcgain(tfr) - r;
figure
hAx = axes('XLim',[-250 250], 'YLim',[-250 250], ...
    'Drawmode','fast', 'NextPlot','add');
axis(hAx, 'equal')

p = calculateEllipse(0, 0, a(1), c(1), step1(1));
hLine = line('XData',p(:,1), 'YData',p(:,2), 'EraseMode','xor',  ...
    'Color','r', 'LineWidth',3);

for i=1:numel(step1)
    p = calculateEllipse(0, 0, a(i), c(i), step1(i));
    set(hLine,'XData',p(:,1), 'YData',p(:,2))  %# update X/Y data
%     pause(.10)                                 %# slow down animation
    drawnow                                    %# force refresh
    if ~ishandle(hLine), return; end           %# in case you close the figure
end
%% Now to make it 3d x2!!!!

wn = 0.06;
z = 0.1;
tfr = tf([5*wn], [1 2*z*wn wn^2]);
step1 = linspace(20,1200,600);
[r, t] = step(tfr, step1);
%r = r*4;
a = r;
c = 2*dcgain(tfr) - r;
% [l, l2, H] = cylinder(r, 35);
figure
hAx = axes('XLim',[-250 250], 'YLim',[-250 250], ...
    'Drawmode','fast', 'NextPlot','add');
axis(hAx, 'equal')
g= (1:10)';
h = repmat(g,1,36);
% for j = 1:36
    p = calculateEllipse(0, 0, a(1), c(1), step1(1));
    hLine = line('XData',p(:,1), 'YData',p(:,2),'Zdata', h(1,:),'EraseMode','xor',  ...
        'Color','r', 'LineWidth',3);
    view(3)
    % surf(p(:,1),p(:,2), H)

    for i=1:numel(step1)
        for j = 1:9
            hLine = line('XData',p(:,1), 'YData',p(:,2),'Zdata', h(j,:),'EraseMode','xor',  ...
            'Color','r', 'LineWidth',3);
            p = calculateEllipse(0, 0, a(i), c(i), step1(i));
            set(hLine,'XData',p(:,1), 'YData',p(:,2))  %# update X/Y data
            view(3)
        %     pause(.10)                                 %# slow down animation
        %     surf(p(:,1),p(:,2), H)
            drawnow                                    %# force refresh
            if ~ishandle(hLine), return; end   %# in case you close the figure
        end
    end
% end

%% Now to make it 3d x3!!!!

wn = 0.06;
z = 0.1;
tfr = tf([5*wn], [1 2*z*wn wn^2]);
step1 = linspace(20,1200,600);
[r, t] = step(tfr, step1);
%r = r*4;
a = r;
c = 2*dcgain(tfr) - r;
% [l, l2, H] = cylinder(r, 35);
figure
hAx = axes('XLim',[-250 250], 'YLim',[-250 250],'ZLim', [-250 250], ...
    'Drawmode','fast', 'NextPlot','add');
axis(hAx, 'equal')
view(3)
g= (1:10)';
h = repmat(g,1,36);
% for j = 1:36
    p = calculateEllipse(0, 0, a(1), c(1), step1(1));
    plot(p(:,1), p(:,2),'Color','r', 'LineWidth',3);
    for i=1:numel(step1)
%         for j = 1:9
            hLine = line('XData',p(:,1), 'YData',p(:,2),'EraseMode','xor',  ...
             'Color','r', 'LineWidth',3);
            p = calculateEllipse(0, 0, a(i), c(i), step1(i));
            plot(p(:,1),p(:,2));
            drawnow                                    %# force refresh
            if ~ishandle(hLine), return; end   %# in case you close the figure
%         end
    end
% end


%% cylinder w 3d plot
% clc; clear all ;
k = p(:,1) ;
l = p(:,2);
z = (1:36)';
z = repmat(z,1,36);


for i = 1:36
    plot3(k,l,z(i,:));
    hold on
end
drawnow
 pause(0.2)

 %% 3d x4??
wn = 0.06;
z = 0.1;
tfr = tf([5*wn], [1 2*z*wn wn^2]);
step1 = linspace(0,1200,600);
[r, t] = step(tfr, step1);
%r = r*4;
a = r;
c = 2*dcgain(tfr) - r;

figure
hAx = axes('XLim',[-250 250], 'YLim',[-250 250], ...
    'Drawmode','fast', 'NextPlot','add');
axis(hAx, 'equal')

p = calculateEllipse(0, 0, a(1), c(1), step1(1));
x = p(:,1) ;
y = p(:,2);
z = (1:36)';
z = repmat(z,1,36);

for i=1:numel(step1)
    clear var p;
    p = calculateEllipse(0, 0, a(i), c(i), step1(i));
    x = p(:,1) ;
    y = p(:,2);%# update X/Y data
%     pause(.10)                                 %# slow down animation
    for j = 1:36
        plot3(x,y,z(j,:));
        hold on
    end
%     pause(0.2)
%     drawnow %# force refresh
    
    if ~ishandle(hLine), return; end           %# in case you close the figure
end

%% 3dx5
z = (1:36)';
z = repmat(z,1,36);
figure
for i=1:numel(step1)
    p = calculateEllipse(0, 0, a(i), c(i), step1(i));
    x = p(:,1) ;
    y = p(:,2);%# update X/Y data
%     pause(.10)                                 %# slow down animation
    for j = 1:36
        plot3(x,y,z(j,:)); 
    end
    plot(x,y)
    pause(0.2)
    drawnow %# force refresh
    if ~ishandle(hLine), return; end           %# in case you close the figure
end

%% try a symmetric cylinder
% z = (1:36)';
% z = repmat(z,1,36);
figure
for i=1:numel(step1)
    p = calculateEllipse(0, 0, a(i), c(i), step1(i));
    x = p(:,1) ;
    y = p(:,2);%# update X/Y data
%     pause(.10)                                 %# slow down animation
    [c1, c2, c3] = cylinder(x);
    surf(c1, c2, c3)
%     plot(x,y)
    pause(0.2)
    drawnow %# force refresh
    if ~ishandle(hLine), return; end           %# in case you close the figure
end

