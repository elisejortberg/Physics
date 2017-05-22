function sum = R_f(n)
    %Elise Jortberg PHYS 4606 HW 1
    %input: n (order of Riemann)
    %output: Riemann zeta function sum
    
    format long
    %initialize
    if n ==2
        %for n=2 ML runs out of memory with estimated N. Slow convergence
        N=1000000; 
        sum = 1;
    else
        %part a: estimate N based on order (n) within error of 10^-11
        N= round((10^(-11)*(n-1))^(1/(-n+1))); 
        if n==3
            sum = 0.25; %initialize according to modified series
        else
            sum =0;
        end
    end
    
    for i=[N:-1:1]
        %part b: if special cases of slow convergence for n=2,3 then
        %use modified series to improve rate of convergence
        if n==3
            sum = sum + (3*i+2)/((i^3)*(i+1)*(i+2));
        elseif n==2
            sum = sum + 1/((i^2)*(i+1));
        else
            sum = sum + 1/(i^n); %Riemann zeta function
        end
    end %end for loop
end
