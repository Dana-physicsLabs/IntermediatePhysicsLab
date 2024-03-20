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



%% Symbolic Uncertainty 
% Transposed from Charlie Tribbles excellent python code poorly by Dana.
% Any mistakes are Dana's, any good is Charlie's. 
% Need Symbolic Math toolbox
% To check and see if you have it please type 'ver' into your Command
% Window. WPI does have a license for it, so you can get it. 


syms r
syms l
syms c

r_val = 10.0;         %# resistance
l_val = 1e-3 ;       %# inductance
c_val = 3e-6  ;      %# capacitance

r_error = r_val * 0.05 ;   % # uncertainty of resistance
l_error = l_val * 0.10  ;  % # uncertainty of inductance
c_error = c_val * 0.10   ; % # uncertainty of capacitance

gamma = r/l;
omega = sqrt(1/(l * c) - (r^ 2)/(4 * (l ^ 2))) ;   %# eq for omega of the damping force

partial_gr = diff(gamma, r);     % # partial derivative of gamma with respect to r
partial_gl = diff(gamma, l) ;    % # partial derivative of gamma with respect to l

% So above we just did our math, and now we are going to evaluate it
f(r,l) = gamma;
Gamma_Fun = matlabFunction(f, 'vars', {[r,l]}); %This makes it a function we can evaluate
f(r,l,c) = omega;
Omega_fun = matlabFunction(f, 'vars', {[r,l,c]}); %This makes it a function we can evaluate

f(r,l,c) = partial_gr;
Partial_gr_fun = matlabFunction(f, 'vars', {[r,l,c]});
f(r,l,c) = partial_gl;
Partial_gl_fun = matlabFunction(f, 'vars', {[r,l,c]});


all_error = [r_error,l_error,c_error]; 
gamma_error = ((Partial_gr_fun(all_error) * r_error)^2 + (Partial_gl_fun(all_error)*l_error)^2)^(1/2);

%gamma_error = ((partial_gr.evalf(subs={r: r_val, l: l_val, c: c_val}) * r_error) ^ 2 + 
%(partial_gl.evalf(subs={r: r_val, l: l_val, c: c_val}) * l_error) ^ 2) ^ 0.5


%Back to just the symbolic stuff
partial_or = diff(omega, r) ;    % # partial derivative of omega with respect to r
partial_ol = diff(omega, l) ;    % # partial derivative of omega with respect to l
partial_oc = diff(omega, c);     % # partial derivative of omega with respect to c


%Back to evaluation stuff 
f(r,l,c) = partial_or;
Partial_or_fun = matlabFunction(f, 'vars', {[r,l,c]});

f(r,l,c) = partial_ol;
Partial_ol_fun = matlabFunction(f, 'vars', {[r,l,c]});

f(r,l,c) = partial_oc;
Partial_oc_fun = matlabFunction(f, 'vars', {[r,l,c]});


%Final Evaluation Thing
omega_error = ((Partial_or_fun(all_error) * r_error) ^ 2 + (Partial_ol_fun(all_error) * l_error)^2 + (Partial_oc_fun(all_error) * c_error)^2)^(1/2)

%omega_error = ((partial_or.evalf(subs={r: r_val, l: l_val, c: c_val}) * r_error) ** 2 \
%             + (partial_ol.evalf(subs={r: r_val, l: l_val, c: c_val}) * l_error) ** 2 \
%             + (partial_oc.evalf(subs={r: r_val, l: l_val, c: c_val}) * c_error) ** 2) ** 0.5

%% Functions in Matlab are always at the bottom 

function [x,dx] = kinematic(time, x0,v0,a)
    x = x0 + (v0*time) + (.5*a*(time.^2));
    dx = 2;
    dv = 2;
    da = 1;
    dx = sqrt((dx).^2 + ((dv*time).^2) + ((.5*da*(time.^2)).^2));
end


