function Bes_comp_n = Jortberg_Bessel(N, x)
    %n>>N, n>>xo
    n=100;
    bes= [];
    for nn =N:1:n
        bes_t = besselj(nn, x);
        bes = [bes, bes_t];
    end
    
    alpha = 0.0001;
    J_next = 0;
    J_now = alpha;
    Bes_comp = [J_next, J_now]; %J100, J99
    for nn=n-2:-1:0
        J_prior = (2*nn/x)*J_now - J_next; %recursion relation
        Bes_comp = [Bes_comp, J_prior];

        J_next = J_now;
        J_now = J_prior; 
    end
    Bes_comp = fliplr(Bes_comp);
    
    %find normalization factor
    sum=0;
    limit = floor(length(Bes_comp)/2);
    for m=2:1:limit
        sum = sum + Bes_comp(2*m); 
    end
    factor = 1/(2*sum+Bes_comp(1)); %1/(2*sum(J_2m(x))+J_0(x))
    Bes_comp = Bes_comp*factor;
    Bes_comp_n = Bes_comp(N+1); %+1 since Bes_comp list begins from 0
    
    SHOW=0;
    if SHOW
        disp(Bes_comp_n);
        disp(bes(1));
        figure; plot(N:100, Bes_comp(N+1:101)); hold on, 
        plot(N:100, bes); title(['J(n,x) = J(', num2str(N),',', num2str(x),')']);
        xlabel('n'); ylabel(['Jn for x=', num2str(x)])
    end
end
