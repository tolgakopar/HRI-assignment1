%% CutOff Frequency.
% Both findings follow the same logic.
% 1) Interpolate both phase and frequency vector.
% 2) Find where phase is below -180 for the first time.
% 3) Get the freq at that point.
function W_cutoff = cutoff_freq(W, Apmn)
    W2   = interp(W, 100); 
    Apmn2 = interp(Apmn, 100);
    below_180_idx = find(Apmn2 <= -180, 1, 'first'); % Obviously the point above
    % Should be in index below_180_idx - 1.

    % If no phase reaches -180, gain margin is theoretically infinite
    if isempty(below_180_idx)
        W_cutoff = Inf;
        return;
    end

   % To calculate the phase margin/gain margin just get the
   % difference between the corresponding boint and the -180deg/0.
   
   W_cutoff = W2(below_180_idx);
end