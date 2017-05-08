function [] = VersionData(version) %Is this a function?
%%Write code name, version number, author names, date and time of execution
fprintf('Code Name: 2D Diffusion Equation Solver\n')
fprintf('Version Number: %s\n',version) %Pull from main module?
fprintf('Author Names: Joseph Corvino and Theresa Green\n')
fprintf('Date: %s\n',date)
time = clock;
fprintf('Time of Execution: %g:%02.0f:%0.4g\n',time(4),time(5),time(6))
end
