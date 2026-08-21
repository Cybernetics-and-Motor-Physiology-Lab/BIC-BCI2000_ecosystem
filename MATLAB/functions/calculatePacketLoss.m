function packetLoss = calculatePacketLoss(state_info, trace_stim)
    % Function to calculate packet loss in recording
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Inputs: 
    %   Recquired: -> state_info: -struct, obtained by load_bcidat.m functin
    %   Optional: -> trace_stim: -boolean, returns packet loss only for
    %   stimulation period
    % Output:
    % -> packetLoss: -double, packet loss expressed in percent
    % Frederik Lampert 2024
    N = length(state_info.Running); % Length of the recording
    switch nargin
        case 1
            N_lost = sum(state_info.ImplantLostSample);
            packetLoss = (N_lost/N)*100;
        case 2
            N_lost = sum(state_info.ImplantLostSample(state_info.ImplantStimulation==1));
            N_stim = sum(state_info.ImplantStimulation);
            packetLoss = (N_lost/N_stim)*100;
    end
end