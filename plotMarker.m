%% ------------------------------------------------------------------------
%% the timer callback function definition
function plotMarker(...
    obj, ...            % refers to the object that called this function (necessary parameter for all callback functions)
    eventdata, ...      % this parameter is not used but is necessary for all callback functions
    player, ...         % audioplayer object to the callback function
    h, ...              % marker handle
    fs)                 % sampling frequency

% check if sound is playing, then only plot new marker
if strcmp(player.Running, 'on')

    % update marker
    h.XData = player.CurrentSample*[1 1]/fs;

end