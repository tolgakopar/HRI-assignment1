%% CutOff Frequency and Phase Margins.
% Both findings follow the same logic.
% 1) Interpolate both gain and frequency vector.
% 2) Find where gain is below 1 for the first time.
% 3) Get the freq at that point.
function [W_cutoff, phase_margin] = cutoff_freq(W, Mpmn, Apmn)
    W2   = interp(W, 100); 
    Mpmn2 = interp(Mpmn, 100);
    Apmn2 = interp(Apmn, 100);
    if (Mpmn2(1) > 1)
        below_zero = find(Mpmn2 <= 1, 1, 'first'); 
    else
        below_zero = find(Mpmn2 >= 1, 1, 'first');
    end

    if isempty(below_zero)
        W_cutoff = Inf;
        phase_margin = 0;
        return;
    end

   % To find the cutoff frequency just get the corresponding frequency for
   % amplitude = 1.
   W_cutoff = W2(below_zero);
   phase_margin = abs(Apmn2(below_zero) + 180);

end
