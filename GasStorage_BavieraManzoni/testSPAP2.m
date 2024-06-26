% Define the data points
x = [1, 2, 3, 4, 5, 6];
y = [5, 3, 2, 4, 9, 4];

% Define the order of the B-spline (cubic spline)
order = 4;

% Define the number of polynomial pieces (choose this to avoid the error)
numPieces = 3; % This means we will have 3 polynomial pieces in the spline

% Fit the B-spline to the data
sp = spap2(numPieces, order, x, y);

% Evaluate the B-spline at the given x points
val = fnval(sp, x);

% Display the results
disp('Fitted values:');
disp(val);

% Plot the original data and the B-spline fit
figure;
plot(x, y, 'bo', 'DisplayName', 'Original Data'); % Original data points
hold on;
xx = linspace(min(x), max(x), 100); % Generate more points for a smooth plot
yy = fnval(sp, xx); % Evaluate the B-spline at these points
plot(xx, yy, 'r-', 'DisplayName', 'B-spline Fit'); % B-spline fit
legend;
title('B-spline Regression');
xlabel('x');
ylabel('y');
