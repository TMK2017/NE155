%% Green-Corvino Master script
% Controls Everything
% Times everything
% Runs GUI

[version_flag, error_message ] = GC_Version();


if version_flag
    display(error_message)
    return
end