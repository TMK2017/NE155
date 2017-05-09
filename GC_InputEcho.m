function [  ] = GC_InputEcho(D,S,sigma,x,y, filename , version_message )
%Print the input data for each cell to the output file

%% Checks if output file alread exists and if it is '.xlsx' format
if exist(filename,'file') == 2
    error('Output file already exists. Please delete or pick a new output file name')
end

if ~all(filename(end-4:end) == '.xlsx')
    error('Output file must be a ''.xlsx'' extension')
end

%% Print version information, D, S, Sigma, X, and Y to output file
warning('off','all')
xlswrite(filename,version_message,'Details')
xlswrite(filename,D,'D')
xlswrite(filename,S,'S')
xlswrite(filename,sigma,'Sigma')
xlswrite(filename,x,'X')
xlswrite(filename,y','Y')
warning('on','all')

end

