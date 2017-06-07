function [triangulation, holes] = findTheHoles(XY,S,M,T)
%% Definition
% findTheHoles is a 2D mesh reconstruction tool which automatically
%    identifies holes in a points cloud. 
%
%% Usage
% Input:
%  XY= Nx2 matrix with point coordinates.
% Optional arguments:
%  S = Critical area ratio (real and positive number). Knots with a value 
%      higher than the ratio of maximum to minimum neighbouring polygon
%      area are identified as hole-border points. Default value is 3.
%  M = Flag for manual editing. Default value is 0.
%      0: No manual editing.
%      1: Manual editing of hole-border points. After the first calculation 
%         a figure with the point cloud is opened. Automatically recognized 
%         points are marked with an 'x'. The user can add/remove 
%         not/wrongly recognized border points selecting them with the 
%         brush and clicking the appropiate button.
%  T = Triangulation (Nx3 matrix with polygon vertices). If not given the
%      delaunay triangulation of X,Y is used.
% 
% Output:
%  triangulation
%     .Points: nx2 matrix with points' coordinates.
%     .ConnectivityList: nx3 matrix with the triangulation polygons.
%  holes: Cell array. Each cell element is an nx1 matrix with ordered
%     hole border point indices. 
%   
%% Method
% To identify the holes a triangulation of the points cloud is needed. If
% none is given, the delaunay triangulation is used. Points where the
% ratio between maximum and minimum neighbouring polygon area exceeds a
% critical value are identified as hole-border points. Polygons for which 
% all vertices are identified as hole-border points are deleted.
% The adjacency matrix of the resulting mesh is clustered to get the
% holes. 
%
%% Issues
% The clustering algorithm isn't fully developed yet and has issues with 
% interconnected holes among others.
%
%% Examples
% See file examples_findTheHoles.m
%
%% Author
% Gypaets
%
%% Function
% Default values
disp('Running findTheHoles.m...')
switch nargin
    case 4
        true;
    case 3
        T = delaunay(XY);
    case 2
        T = delaunay(XY);
        M = 0;
    case 1
        T = delaunay(XY);
        M = 0;
        S = 3;
    otherwise
        error('Wrong arguments number')
end

X = XY(:,1); % x-coordinate Points
Y = XY(:,2); % y-coordinate Points
nV = length(X); % Number of points/vertices
nP = length(T); % Number of polygons in triangulation

% Global variables
global borderPoints boundaryPlot

%-------------------------------------------------------------------------%
% Mesh creation and boundary points recognition
% Polygon area
areas = polyarea(X(T)',Y(T)')';

% Vertex-Polygon matrix
mVerPol = sparse(T(:),repmat(1:nP,[1 3]),1,nV,nP);
mVPA = bsxfun(@times,mVerPol,areas')';

% Maximum and minimum polygon area at each vertex
maxA = max(mVPA);
minA = max(spfun(@(x) x.^-1,mVPA)).^-1;

% Ratio between maximum and minimum polygon area at each vertex
arRat = maxA./minA;

% Select vertices with areas ratio arRat greater than the threshold
borderPoints = find(arRat > S)';

% Manual editing
if M == 1
    fig = figure();
    set(fig,'Position',[680,558,560,420]);
    boundaryPlot = plot(X(borderPoints),Y(borderPoints),'x');
    hold on
    allPlot = plot(X,Y,'.');
    brush on
    controlpanel = uipanel('Title','Control buttons',...
        'Position',[0.01,0.55,0.27,0.375],...
        'Units','pixels',...
        'parent',fig);
    
    % Button to add brushed points to boundary points list
    addBP_b = uicontrol(1,'Style', 'pushbutton',...
        'String', 'Add boundary points',...
        'Position',[10, 90, 130, 30],...
        'Units','pixel',...
        'Callback',{@addBP, allPlot, X, Y},...
        'parent',controlpanel);
    
    % Button to remove brushed points from boundary points list
    falseBP_b = uicontrol(1,'Style', 'pushbutton',...
        'String', 'False positive',...
        'Position',[10, 50, 130, 30],...
        'Units','pixel',...
        'Callback',{@falseBP, allPlot, X, Y},...
        'parent',controlpanel);
    
    % Button to close the figure
    close_b = uicontrol(1,'Style', 'pushbutton',...
        'String', 'Close figure',...
        'Position',[10, 10, 130, 30],...
        'Units','pixel',...
        'Callback','close(gcf)',...
        'parent',controlpanel);
    
    waitfor(fig);
end

% Sparse array with border points indices
borderPointsS = sparse(1,borderPoints,1,1,nV);

% Number of border vertices in every polygon
bordVert = borderPointsS*mVerPol;

% Get index of polygons with three border vertices
detPol = find(bordVert == 3);

% Remove polygons connected to three border vertices
T(detPol,:) = [];

% Save resulting triangulation
triangulation.Points = XY;
triangulation.ConnectivityList = T;


%-------------------------------------------------------------------------%
% Identify border polygons and cluster them in a cell array
% All points adjacency matrix
[lr, lc] = find(sparse(reshape(T(:,[1 1 2]).',[],1),reshape(T(:,[2 3 3]).',[],1),1,nV,nV));
adjacencyMatrix = sparse(lr,lc,1,nV,nV)+sparse(lr,lc,1,nV,nV).';

% Border-points adjacency matrix
borderAdjacencyMatrix = adjacencyMatrix(borderPoints,borderPoints);
[lr, lc] = find(borderAdjacencyMatrix);
borderAdjacencyMatrix = full(sparse(lr,lc,1,length(borderPoints),length(borderPoints)));

% Remove dead ends from adjacency matrix and actualize border points:
% loop until every border points is connected to exactly two other points.
while not(prod(sum(borderAdjacencyMatrix) == 2))
    % Edges number at each point
    edgesNumber = sum(borderAdjacencyMatrix);
    
    % Dead ends list
    deadEnds = find(edgesNumber <= 1);
    
    % Workaround for interconnected geometries
    % Run after loop repetition without changes
    if isempty(deadEnds)
        deadEnds = find(edgesNumber > 2);
        display('findTheHoles: complex geometries! Hole clustering results probably wrong.')
    end
    
    % New adjacency matrix
    borderAdjacencyMatrix(deadEnds,:) = [];
    borderAdjacencyMatrix(:,deadEnds) = [];
    
    % Actualized border points
    borderPoints(deadEnds) = [];
end

% Initial values to cluster the points
lHoleNumber = 0;
processedKnots = zeros(1,length(borderAdjacencyMatrix));

% Run loop through all knots to determine the hole they're in
while not(prod(processedKnots))
     % Start new hole
    lHoleNumber = lHoleNumber+1;
    
    % Start new hole from first free point found
    notProcessed = find(not(processedKnots));
    lKnot = notProcessed(1);
    
    neighbourKnots = find(borderAdjacencyMatrix(:,lKnot));
    lHoleKnotsList = lKnot;
    
    % Follow picked point until hole is closed
    while not(prod(ismember(neighbourKnots,lHoleKnotsList)))
        n1 = neighbourKnots(1);
        n2 = neighbourKnots(2);
        
        % Check what neighbour point is already in the list and take the
        % other one. If neither of them is in 'lHoleKnotsList' (the case
        % when processing the first point) take the first one.
        if ismember(n2,lHoleKnotsList)
            lHoleKnotsList = [lHoleKnotsList n1];
            processedKnots(lKnot) = 1;
            lKnot = n1;
        elseif ismember(n1,lHoleKnotsList)
            lHoleKnotsList = [lHoleKnotsList n2];
            processedKnots(lKnot) = 1;
            lKnot = n2;
        elseif length(lHoleKnotsList) == 1
            lHoleKnotsList = [lHoleKnotsList n1];
            processedKnots(lKnot) = 1;
            lKnot = n1;
        end

        neighbourKnots = find(borderAdjacencyMatrix(:,lKnot));
    end
    
    % Mark last hole-point as processed too
    processedKnots(lKnot) = 1;
    
    % Write ordered hole knots to result
    holes{lHoleNumber} = borderPoints(lHoleKnotsList); 
end

%------------------------------------------%
% Plot result
if M == 1
    % Plot mesh
    fig = figure();
    triplot(T,X,Y);
    hold on
    plot(X,Y,'.',X(borderPoints),Y(borderPoints),'x');
    
    % Plot hole borders
    figure();
    hold on;
    for i = 1:length(holes)
        lX = X(holes{i});
        lY = Y(holes{i});
        plot([lX; lX(1)],[lY; lY(1)]);
    end
    waitfor(fig)
end

end

%------------------------------------------------------------------%
function addBP(~, ~, lineSeries1, X, Y)
% Adds boundary points to boundary points list
% Handle of invoking object and event data are ignored
%
% Global variables
global borderPoints boundaryPlot

% Get index of brushed points
selectedPoints = find(get(lineSeries1,'BrushData'))';

% Add selected points to boundary points list
borderPoints = unique([borderPoints; selectedPoints]);

% Update boundary point plot
delete(boundaryPlot),
boundaryPlot = plot(X(borderPoints),Y(borderPoints),'x');

end

%-----------------------------------------------------------------%
function falseBP(~, ~, lineSeries1, X, Y)
% Removes selected points from boundary points list
% Handle of invoking object and event data are ignored
%
% Global variables
global borderPoints boundaryPlot

% Get index of brushed points
selectedPoints = find(get(lineSeries1,'BrushData'));

% Remove brushed points from the boundary-points list 
spIndex = find(ismember(borderPoints,selectedPoints));
borderPoints(spIndex) = [];

% Update boundary point plot
delete(boundaryPlot);
boundaryPlot = plot(X(borderPoints),Y(borderPoints),'x');

end