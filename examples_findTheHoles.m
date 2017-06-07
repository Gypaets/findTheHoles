%%-------------------------------------------------------------------------
%% Example 1
% Alice ran a cfd-simulation around two cylinders and gave Bob a *.csv file 
% with the points coordinates and the airflow speed at each point. Bob 
% wants to calculate the particle acceleration at each point and therefore 
% needs to mesh the points set. Simply calculating a delaunay triangulation 
% won't work because points on the upper and lower part of the cylinder 
% would be connected, which physically does't make any sense. Furthermore, 
% Bob needs to determine which points are on the surface of each cylinder.
%
% Run the following sections to see the raw data, run the program and see
% the results.
%% Create Alice's data
[X, Y] = meshgrid(1:100,1:90);
X = X(:)+2*(0.5-rand(9000,1))/20;
Y = Y(:)+2*(0.5-rand(9000,1))/20;
MAT = 9000;
x0 = 50;
y0 = 50;
r0 = 20;
c = zeros(MAT,1);
for i = 1:MAT
    if ((X(i)-x0)^2+(Y(i)-y0)^2)<r0^2
        X(i) = 0;
        Y(i) = 0;
        c(i) = 1;
    end
end
x1 = 80;
y1 = 20;
r1 = 10;
for i = 1:MAT
    if ((X(i)-x1)^2+(Y(i)-y1)^2)<r1^2
        X(i) = 0;
        Y(i) = 0;
        c(i) = 1;
    end
end
X(find(c)) = [];
Y(find(c)) = [];
T = delaunay(X,Y);
plot(X,Y,'.')

%% Bob runs findTheHoles
% The resulting mesh is saved in 'mesh'.
% The indices of each knot at the cylinder i are saved in cylinders{i}
[mesh, cylinders] = findTheHoles([X Y]);

%% Plot results
% Triangulation
figure();
triplot(mesh.ConnectivityList,mesh.Points(:,1),mesh.Points(:,2));

% Automatically recognized cylinders
figure();
hold on;
for i = 1:length(cylinders)
    lX = X(cylinders{i});
    lY = Y(cylinders{i});
    plot([lX; lX(1)],[lY; lY(1)]);
end
clear lX lY i

%% Manual update
% Bob doesn't like some of the calculated values and runs findTheHoles
% again in manual mode.
[mesh2, cylinders2] = findTheHoles([X Y],1.5,1);