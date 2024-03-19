% Variables are vectors/matricies:
a = [1 2 3 4 5]
b = [1;2;3;4;5]
c = [1 2; 3 4];

% Shortcuts to creating vectors/matricies:
x = 1:10;
y = 2*x;

%Basic plotting:
figure; plot(x,y);
figure; plot(x,y,'r.');

%Functions:
t = 0:.1:10;
[y,dy] = kinematic(t,0,20,-9.8);
figure; errorbar(t,y,dy);
figure; plot(t,y); hold on;
for i = 0:20:200
    plot(t,kinematic(t,i,20,-9.8));
end

%Importing data:
dataTable = readtable('foam_1click_10min_clamp.xlsx');
times = dataTable{:,1};
heights = dataTable{:,2};

%Advanced plotting:
figure; hist(heights);
xlabel('Height (cm)'); ylabel('Number of trials');
figure; plot(times,heights,'o');
model = fitlm(times,heights) % See https://www.mathworks.com/help/stats/coefficient-standard-errors-and-confidence-intervals.html
figure; plot(model)

function [x,dx] = kinematic(time, x0,v0,a)
    x = x0 + (v0*time) + (.5*a*(time.^2));
    dx = 2;
    dv = 2;
    da = 1;
    dx = sqrt((dx).^2 + ((dv*time).^2) + ((.5*da*(time.^2)).^2));
end
