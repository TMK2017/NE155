%% Green-Corvino Master script
% Controls Everything
% Times everything
% Runs GUI

% Start afresh :)
clear all
close all
clc

% Verifty versions of each subroutine
[version_flag, error_message ] = GC_Version();

% Stop GC_Master if versions are not up-to-date
if ~version_flag
    display(error_message)
    return
end



display('GC_Master Sucessful')