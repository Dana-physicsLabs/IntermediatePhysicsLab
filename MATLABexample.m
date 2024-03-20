% close all 

%%  Significant Figures Check in 

% What are the significant figures on these functions below? 

a = 1.2*0.00004032; 

b = 100/(58^2);

c = (1/2) * pi;

d = a/c^3;


%clear all

%% Standard deviaton
% https://www.mathworks.com/help/matlab/ref/std.html 
%Importing data:
dataTable = readtable('foam_1click_10min_clamp.xlsx');
times = dataTable{:,1};
heights = dataTable{:,2};

t_stdev = std(times);
h_stdev = std(heights);

%data_stdev = std(dataTable); %Why do you think it breaks? Hint: Python
%also has this issues, especailly if you love Pandas
all_stev = std(times,heights); %Why does this one work? What does it mean physically?

h_mean = mean(heights);

h_stderror = h_stdev/ sqrt(h_mean);
%% Matlab Functions 

% Variables are vectors/matricies:
a = [1 2 3 4 5]
b = [1;2;3;4;5]
c = [1 2; 3 4];

% Shortcuts to creating vectors/matricies:
x = 1:10;
y = 2*x;

%Basic plotting:
figure(1); plot(x,y);
figure(2); plot(x,y,'r.');

%Functions:
t = 0:.1:10;
[y,dy] = kinematic(t,0,20,-9.8);
figure(3); errorbar(t,y,dy);
figure(4); plot(t,y); hold on;
for i = 0:20:200
    plot(t,kinematic(t,i,20,-9.8));
end



%Advanced plotting:
figure(5); hist(heights);
xlabel('Height (cm)'); ylabel('Number of trials');
figure(6); plot(times,heights,'o');
model = fitlm(times,heights) % See https://www.mathworks.com/help/stats/coefficient-standard-errors-and-confidence-intervals.html
figure(7); plot(model)

function [x,dx] = kinematic(time, x0,v0,a)
    x = x0 + (v0*time) + (.5*a*(time.^2));
    dx = 2;
    dv = 2;
    da = 1;
    dx = sqrt((dx).^2 + ((dv*time).^2) + ((.5*da*(time.^2)).^2));
end
