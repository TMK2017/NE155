function [version_message] = GC_VersionData(v) 
%%Write code name, version number, author names, date and time of execution
% fprintf('Code Name: 2D Diffusion Equation Solver\n')
% fprintf('Matlab Version Number: %s\n',version) 
% fprintf('GC_Master.m Version Number: %\n',v) %Pull from main module
% fprintf('Author Names: Joseph Corvino and Theresa Green\n')
% fprintf('Date: %s\n',date)
% time = clock;
% fprintf('Time of Execution: %g:%02.0f:%0.4g\n',time(4),time(5),time(6))

version_message = {'Code Name: 2D Diffusion Equation Solver'};
version_message{2,1} = ['Matlab Version Number: ',version];
version_message{3,1} = ['GC_Master.m Version Number: ',v];
version_message{4,1} = 'Author Names: Joseph Corvino and Theresa Green';
version_message{5,1} = ['Date: ',date];
time = fix(clock);
version_message{6,1} = ['Time of Execution: ',num2str(time(4)),':',num2str(time(5)),':',num2str(time(6))];


end
