function [alls] = run_Ode(M,L,tspan)
    s = struct();
    
    % initialize
    opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1);

    %tspan = linspace(0,300,30e1);

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
        % define TF
        TF = L.nucloc(k,1); 
        % parameters        
        kact = s.q0(1);
        kinact = s.q0(2);
        b1 = s.q0(3);
        g1 = s.q0(4); 
        b2 = s.q0(5); 
        g2 = s.q0(6); 
        b0 = s.q0(7);
        
        % solve pka system
        xinit  = []; % assume here the system starts in the OFF state
        xinit(1) = (kact.*TF)./(kact.*TF+ kinact);%1/(((kinact/kact)*((TF^n+kd)/TF^n))+1); %mrna
        xinit(2) = (b0+b1.*(kact.*TF./(kact.*TF+kinact)))./g1; %(b1/g1)*1/(((kinact/kact)*((TF^n+kd)/TF^n))+1); %mrna
        xinit(3) = (b2./g2).*(b0+b1.*(kact.*TF./(kact.*TF+kinact)))./g1;%(b2/g2)*(b1/g1)*1/(((kinact/kact)*((TF^n+kd)/TF^n))+1); %prot
        xinit(4) = 0;
        xinit(5) = 0;
        [t_x,x]  = ode45(@(t,x) computeRates_ponIntegral(t,x,s), tspan, xinit, opts);
        
        s.pon     = x(:,1);
        s.mrna    = x(:,2); 
        s.prot    = x(:,3);
        s.intPon  = x(:,4);
        s.right   = x(:,5);

        % store data
        alls{k} = s;
        %alls{k} = s.prot;
    end
    
end



