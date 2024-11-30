%% Phase and Gain Margin
% Both findings follow the same logic.
% 1) Find the point that gain/phase becomes zero/reaches -180deg.
% 2) Find the corresponding point in the other plot (phase/mag
%   respectively)
% 3) Find the distance between this point and the 180deg/0 gain.
function phase_margin = Phase_margin(Hsys, W)
    W2   = interp(W, 100); 
    [gain, phase] =  bode(Hsys, W2);
    phase = squeeze(phase);
    gain = squeeze(gain);

    below_zero_idx = find(gain <= 1, 1, 'first'); % Obviously the point above

    % Should be in index below_zero_idx - 1.

    % If no phase reaches -180, gain margin is theoretically infinite
    if isempty(below_zero_idx)
        phase_margin = Inf;
        return;
    end
   % Interpolate between the frequencies and get the gain from the
   % system function.
   phase_margin = abs(phase(below_zero_idx) + 180);
end