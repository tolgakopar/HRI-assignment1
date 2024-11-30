%% Phase and Gain Margin
% Both findings follow the same logic.
% 1) Find the point that gain/phase becomes zero/reaches -180deg.
% 2) Find the corresponding point in the other plot (phase/mag
%   respectively)
% 3) Find the distance between this point and the 180deg/0 gain.
function gain_margin = Gain_margin(Hsys, W)
    W2   = interp(W, 100); 
    [gain, phase] =  bode(Hsys, W2);
    phase = squeeze(phase);
    gain = squeeze(gain);

    below_180_idx = find(phase <= -180, 1, 'first'); % Obviously the point above
    % Should be in index below_180_idx - 1.

    % If no phase reaches -180, gain margin is theoretically infinite
    if isempty(below_180_idx)
        gain_margin = Inf;
        return;
    end

   % To calculate the phase margin/gain margin just get the
   % difference between the corresponding boint and the -180deg/0.
   gain_margin = abs(gain(below_180_idx));
   gain_margin = 20*log(gain_margin);
end

