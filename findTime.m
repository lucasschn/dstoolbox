function [t, Ts] = findTime(motion,c)
N = length(motion.alpha);
if isempty(motion.Ts)
    switch class(motion)
        case 'PitchingMotion'
            % k = omega*c/(2*V)
            if ~isempty(motion.freq)
                freq = motion.freq;
            else
                freq = motion.k*2*motion.V/(2*pi*c);
             end
            if ~isempty(motion.f_pts)
                Ts = motion.f_pts/freq; % rev/sampling interval * time/rev = time/sampling intervail
            else
                error('Ts or f_pts must be set on this data beforehand, e.g. using findSinus.')
            end
        case 'RampUpMotion'
            % r = alphadot*c/(2*V)
            if isempty(motion.Ts)
                error('Computing the time vector without the sampling frequency Ts is impossible for a RampUpMotion object.')
            end
    end
else
    Ts = motion.Ts;
end
    t = reshape(0:Ts:Ts*(N-1),size(motion.alpha));
end
