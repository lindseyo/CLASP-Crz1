function [alls] = run_Ode_2state(M,L,tspan,SSvals)
    s = struct();
    
    % initialize
    %opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);
    opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);

    s.t   = tspan;
    
    % number of iterations
    N    = size(L.nucloc,1);
    
    % initial parameters
    s.q0 = M;
    
    alls = cell(1,N);
    %disp(strcat('runningthru:',num2str(j)))
    for k = 1:N
        %display(strcat('k=',num2str(k)))
        % nuc loc
        s.nucloc = L.nucloc(k,:);
        s.nuclocT = L.times(k,:);
        %disp(strcat('gothru:',num2str(k)))  

        % solve pka system
        xinit  = [];
        xinit(1) = SSvals(1).p0SS; 
        xinit(2) = SSvals(1).poffSS;
        xinit(3) = SSvals(1).ponSS;
        xinit(4) = SSvals(1).mrnaSS;
        xinit(5) = SSvals(1).protSS;
        
        %[t_x,x]  = ode45(@(t,x) ComputeRates_2state(t,x,s), tspan, xinit, opts);

        [t_x,x]  = ode113(@(t,x) ComputeRates_2state(t,x,s), tspan, xinit, opts);
        
        s.all = x;

        % store data
        alls{k} = s;
    end
    
end



