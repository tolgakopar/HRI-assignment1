%% CutOff Frequency.
% Both findings follow the same logic.
% 1) Interpolate both gain and frequency vector.
% 2) Find where gain is below 1 for the first time.
% 3) Get the freq at that point.
function W_cutoff = cutoff_freq(W, Mpmn)
    W2   = interp(W, 100); 
    Apmn2 = interp(Mpmn, 100);
    below_zero = find(Apmn2 <= 1, 1, 'first'); 

    if isempty(below_zero)
        W_cutoff = Inf;
        return;
    end

   % To find the cutoff frequency just get the corresponding frequency for
   % amplitude = 1.
   W_cutoff = W2(below_zero);
end
