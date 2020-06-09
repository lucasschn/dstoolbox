function alpha_crit = computeAlphaCrit(airfoil,alpha_ds0,r)
% For a given airfoil with y-intercept alpha_ds0 for alpha_ds(r) linearfit, returns the a critical angle
% of attack corresponding to the reduced pitch rates r. alpha_crit can then be used in Sheng's stall criterion, according to Sheng 2008. 

        % initialization of the vectors
        alpha_crit = -ones(size(r));
        
        % definition of alpha_crit depending on r (only for Sheng)
        for k=1:length(r)
            if r >= airfoil.r0
                alpha_crit(k) = alpha_ds0;
            else % if r<r0
                alpha_crit(k) = airfoil.steady.alpha_ss + (alpha_ds0 - airfoil.steady.alpha_ss)*r/airfoil.r0;
            end
        end
end