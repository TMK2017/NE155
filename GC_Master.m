function [] = GC_Master(userx,usery,userD,userS,usersigma,filename)
%% Green-Corvino Master Function
% Controls Everything
% Times everything
% Runs GUI
tic; % start timer
v = 'v4.0';

% close figure windows
close all



%Write code name, version number, author names, date and time of execution
version_message = GC_VersionData(v);
[D,S,sigma,x,y,flag] = GC_InputData(userx,usery,userD,userS,usersigma);
if flag
   error('Input data is incorrectly formatted and did not pass GC_InputData.m check') 
end
% Print inputs to file
GC_InputEcho(D,S,sigma,x,y, filename, version_message)

% Solve for phi
phi = GC_DiffSolver(x,y,D,sigma,S);

% plot results
figure
hold on
grid on
view(3)
[X,Y]= meshgrid(x,y);
surf(X,Y,(phi)','EdgeColor','None')
xlabel('X')
ylabel('Y')
zlabel('\phi(x,y)')

% print to output file
warning('off','all')
xlswrite(filename,phi','phi')
warning('on','all')

toc
end