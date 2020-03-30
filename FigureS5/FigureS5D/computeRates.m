function dxdt = ComputeRates(t,x,s)

    % state variables
    x1 = max(x(1),0); 
    x2 = max(x(2),0);
    x3 = max(x(3),0);
    
    % input
    Tnuc = interp1(s.nuclocT, s.nucloc, t); %%% change this!!!!
    
    % parameters
    kact = s.q0(1);
    kinact = s.q0(2);
    b1 = s.q0(3);
    g1 = s.q0(4); 
    b2 = s.q0(5); 
    g2 = s.q0(6); 
    b0 = s.q0(7);
    
    % rates
    dxdt = zeros(3,1);    

    dxdt(1) = kact*(1-x1)*Tnuc - kinact*x1; % dPon/dt
%    dxdt(1) = kactS*(1-x1) + kact*(1-x1)*Tnuc - kinact*x1; % dPon/dt
    dxdt(2) = b0 + b1*x1 - g1*x2; % dmRNA/dt
    dxdt(3) = b2*x2 - g2*x3; % dP/dt  
    
end