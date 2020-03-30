function dxdt = ComputeRates_2state(t,x,s)

    % state variables
    x1 = max(x(1),0); 
    x2 = max(x(2),0);
    x3 = max(x(3),0);
    x4 = max(x(4),0);
    x5 = max(x(5),0); 

    % input
    Tnuc = interp1(s.nuclocT, s.nucloc, t); %%% change this!!!!
    
    % parameters
    kon = s.q0(1);
    koff = s.q0(2);
    b1 = s.q0(3);
    g1 = s.q0(4); 
    b2 = s.q0(5); 
    g2 = s.q0(6); 
    b0 = s.q0(7);
    ron = s.q0(8);
    roff = s.q0(9);
    THR2 = s.q0(10);
    
    % implement thresholding
    if Tnuc >= THR2
        roff = 0;
    else
        roff = roff;
    end
    
    % rates CHECKED THESE EQUATIONS - 3-STATE
    dxdt = ones(5,1);    
    dxdt(1) = roff.*x2 - ron.*x1; % dP0
    dxdt(2) = koff.*x3 + ron.*x1 - (roff+kon*Tnuc).*x2; % dPoff
    dxdt(3) = kon*Tnuc*x2 - koff*x3; % dPon
    dxdt(4) = b0 + b1*x3 - g1*x4; % dmRNA/dt
    dxdt(5) = b2*x4 - g2*x5; % dP/dt
    
% if you add Tnuc to Ron, then it's like you put a treshold on it, if you
% just threshold Ron, then it works
% so right now -- what's the minimal model -- add in 1 additional state,
% then make ron dependent on TF (dependence could be either threshold or
% mass action)
end