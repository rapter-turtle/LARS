function [ data, g, data0 ] = air3D(accuracy)

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 10.0;                  % End time.
plotSteps = 50;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'global';

doMask = 1; %%% every step avoid docking bay. 
doMin = 0;
%--------------------------------------------------------------------------
targetRadius = 1;
V0 = 2.5;
Vd = 1.0;
Wd = 0.3;
Kp = 1.0;
Kd = 1.0;
%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the 3D reachable set.
% displayType = 'contour';
displayType = 'contour';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% Visualize the angular dimension a little bigger.
aspectRatio = [1 1];

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

%---------------------------------------------------------------------------
% Approximately how many grid cells?
%   (Slightly different grid cell counts will be chosen for each dimension.)
Nx = 31;

% Create the grid.
g.dim = 4;
g.min = [ -10;-5; -5; -5];
g.max = [ +2;+5; +5; +5];
g.bdry = { @addGhostExtrapolate; @addGhostExtrapolate; @addGhostExtrapolate; @addGhostExtrapolate};
% Roughly equal dx in x and y (so different N).
g.N = [ Nx;Nx;Nx;Nx];
% Need to trim max bound in \psi (since the BC are periodic in this dimension).
% g.max(3) = g.max(3) * (1 - 1 / g.N(3));
g = processGrid(g);

%---------------------------------------------------------------------------
% Create initial conditions (cylinder centered on origin).
% data = shapeRectangleByCorners(g, [-1;-1;0;0], [1; 1; 0;0]);

% data = shapeCylinder(g, [3,4], [ 0; 0;0; 0], targetRadius);
% data0 = data;

% data1 = shapeRectangleByCorners(g, [-2.5;-1;0;0], [2.5; 2; 0;0]);
% data2 = shapeRectangleByCorners(g, [-1.5;-1;0;0], [1.5; 1; 0;0]);

dataone =  shapeRectangleByCorners(g, [-1;-1.5;-5;-5], [1; 1.5; 5;5]);
% data =  shapeRectangleByCorners(g, [-1;-1;-5;-5], [1; 1; 5;5]);

% dataone = shapeCylinder(g, [3,4], [ 0; 0;0; 0], targetRadius);
data1 = shapeRectangleByCorners(g, [-1;-5.0;-100;-100], [2; 5.0; 100;100]);
data2 = shapeRectangleByCorners(g, [-1;-1.5;-50;-50], [1; 1.5; 50;50]);
data3 = max(data1,-data2);%docking bay
mask = -data3;
data = max(dataone,mask); % initial masking (every step ) growint set과 avoiding set의 여집합에 대한 교집합.
data0 = data;
%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @air3DHamFunc;
schemeData.partialFunc = @air3DPartialFunc;
schemeData.grid = g;

% The Hamiltonian and partial functions need problem parameters.
schemeData.V0 = V0;
schemeData.Vd = Vd;
schemeData.Wd = Wd;
schemeData.Kp = Kp;
schemeData.Kd = Kd;

%---------------------------------------------------------------------------
% Choose degree of dissipation.

switch(dissType)
 case 'global'
  schemeData.dissFunc = @artificialDissipationGLF;
 case 'local'
  schemeData.dissFunc = @artificialDissipationLLF;
 case 'locallocal'
  schemeData.dissFunc = @artificialDissipationLLLF;
 otherwise
  error('Unknown dissipation function %s', dissFunc);
end

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'veryHigh';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  schemeData.derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  schemeData.derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  schemeData.derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  schemeData.derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end


schemeData.mask = mask(:);
schemeData.doMask = doMask;

schemeData.min = data(:);
schemeData.doMin = doMin;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integratorOptions = odeCFLset(integratorOptions,...
                              'postTimestep',@maskmask); % every time step PI is constrained

%---------------------------------------------------------------------------
% Restrict the Hamiltonian so that reachable set only grows.
%   The Lax-Friedrichs approximation scheme MUST already be completely set up.
innerFunc = schemeFunc;
innerData = schemeData;
clear schemeFunc schemeData;

% Wrap the true Hamiltonian inside the term approximation restriction routine.
schemeFunc = @termRestrictUpdate;
schemeData.innerFunc = innerFunc;
schemeData.innerData = innerData;
schemeData.positive = 0;

schemeData.mask = mask(:);
schemeData.doMask = doMask;

schemeData.min = data(:);
schemeData.doMin = doMin;

integratorOptions = odeCFLset(integratorOptions,...
                              'postTimestep',@maskmask);

%---------------------------------------------------------------------------
% Initialize Display
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

% h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);

% camlight right;  camlight left;
% hold on;
% axis(g.axis);
% % daspect(aspectRatio);
% drawnow;

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);



  
  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end


  % 1. z, r 차원에 대해 최소값 계산
  dataMin = min(min(data, [], 4), [], 3);  % size = [x, y]

  % 2. 조건에 맞는 (x, y) 마스크 생성
  mask2D = dataMin < 0;

  % 3. 해당 위치의 x, y 인덱스 추출
  [xIdx, yIdx] = find(mask2D);  % OK! mask2D는 2D니까 오류 없음

  % 실제 x, y 좌표 값
  xVals = g.xs{1}(sub2ind(size(dataMin), xIdx, yIdx));
  yVals = g.xs{2}(sub2ind(size(dataMin), xIdx, yIdx));

  figure;
  scatter(xVals, yVals, 'filled');
  title('x, y where min(data(x,y,:,:)) < 0');
  xlabel('x'); ylabel('y');
  axis equal;
  xlim([-10, 2]);
  ylim([-5, 5]);

  % Get correct figure, and remember its current view.
  % figure(f);
  % [ view_az, view_el ] = view;
  % 
  % % Delete last visualization if necessary.
  % if(deleteLastPlot)
  %   delete(h);
  % end
  % 
  % % Move to next subplot if necessary.
  % if(useSubplots)
  %   plotNum = plotNum + 1;
  %   subplot(rows, cols, plotNum);
  % end

  % % Create new visualization.
  % h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
  % 
  % % Restore view.
  % view(view_az, view_el);
  
end
% Reshape the data matrix into columns (x, y, z, r, V)
xVals = g.xs{1}; % Get x values from grid
yVals = g.xs{2}; % Get y values from grid
zVals = g.xs{3}; % Get z values from grid
rVals = g.xs{4}; % Get r values from grid

% Flatten the data to match (x, y, z, r, V)
xFlat = xVals(:);
yFlat = yVals(:);
zFlat = zVals(:);
rFlat = rVals(:);
VFlat = data(:); % Assuming 'data' contains the velocity V values

% Create a table with the desired format (x, y, z, r, V)
outputTable = table(xFlat, yFlat, zFlat, rFlat, VFlat);

% Write the table to an Excel file
% Construct filename dynamically including Nx, Kp, Kd
filename = sprintf('output_data_Nx%d_Kp%.2f_Kd%.2f.xlsx', Nx, Kp, Kd);

% Write the table to an Excel file with the dynamic filename
writetable(outputTable, filename);


endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = air3DHamFunc(t, data, deriv, schemeData)


checkStructureFields(schemeData, 'grid', 'V0', 'Vd', ...
                                 'Wd', 'Kp', 'Kd');

grid = schemeData.grid;

% implements equation (3.3) from my thesis term by term
%   with allowances for \script A and \script B \neq [ -1, +1 ]
%   where deriv{i} is p_i
%         x_r is grid.xs{1}, y_r is grid.xs{2}, \psi_r is grid.xs{3}
%         v_a is velocityA, v_b is velocityB, 
%         \script A is inputA and \script B is inputB
hamValue = -((-schemeData.V0+grid.xs{3} ).* deriv{1} ...
	     - schemeData.Kp * grid.xs{1} .* deriv{3} ...
         - schemeData.Kd * (-schemeData.V0+grid.xs{3}+schemeData.Vd.*sign(deriv{1})) .* deriv{3} ...
	     + schemeData.Vd* abs(deriv{1}) ...
	     + schemeData.Wd * abs(deriv{3}) ...
         + grid.xs{4}.* deriv{2} ...
	     - schemeData.Kp * grid.xs{2} .* deriv{4} ...
         - schemeData.Kd * (grid.xs{4}+schemeData.Vd.*sign(deriv{2})) .* deriv{4} ...
	     + schemeData.Vd* abs(deriv{2}) ...
	     + schemeData.Wd * abs(deriv{4}));



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = air3DPartialFunc(t, data, derivMin, derivMax, schemeData, dim)


checkStructureFields(schemeData, 'grid', 'V0', 'Vd', ...
                                 'Wd', 'Kp', 'Kd');

grid = schemeData.grid;

switch dim
  case 1
    alpha = abs(schemeData.V0 + ...
                + schemeData.Vd + abs(grid.xs{3}));

  case 3
    alpha = abs(schemeData.Kp * grid.xs{1}) ...
            + schemeData.Wd;
  case 2
    alpha = abs(schemeData.V0 + ...
                + schemeData.Vd + abs(grid.xs{4}));

  case 4
    alpha = abs(schemeData.Kp * grid.xs{2}) ...
            + schemeData.Wd;
  otherwise
    error([ 'Partials for the game of two identical vehicles' ...
            ' only exist in dimensions 1-3' ]);
end

function [yOut, schemeDataOut] = maskmask(~,yIn,schemeData)
checkStructureFields(schemeData, 'doMask','doMin');
if(schemeData.doMask)
    checkStructureFields(schemeData,'mask');
    yOut = max(yIn, schemeData.mask); %every step compare with masking set, it makes reachset never grow in masking set
else
    yOut = yIn;
end

schemeDataOut = schemeData;
if(schemeData.doMin)
    checkStructureFields(schemeData, 'min');
    schemeDataOut.min = min(yOut, schemeDataOut.min);
end
