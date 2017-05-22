function [Uout, L, L_avg, C] = Ising_v2(T, sweeps)
    %This function computes the Ising model for ferromagnetism. The lattice
    %consists of magnetic dipole moments (atomic spins). The energy of the
    %system (the Hamiltonian) is dictated by the spin interactions and
    %their resulting phase transitions. The configuration probability draws
    %from the Boltzmann distribution.
    %For more info see: https://en.wikipedia.org/wiki/Ising_model
    %Inputs
    %T = Temperature
    %sweeps = # runs/updates
    %Outputs
    %U = Energy of system (1xsweeps)
    %L = (N+ - N-)/N spin states
    %C = last state of lattice
    
    %T = 2.2;
    m = 100; %mxm lattice
    B = 0;   %magnetic field
    k = 1;   %k = constant
    mu = 1;  %magnetic moment
    J = 1;   %coupling between neighboring moments
    N= m*m;

    beta = 1/(k*T);

    C = 2*(rand(m,m) > 0.5) -1; %create lattice
    figure(1); imagesc(C); colormap(gray); title('initial');
    hold on;

    U = zeros(1,sweeps); %initialize energy vector
    
    figure(2); hold on;
    xlim([0 sweeps]); ylim([-2.5 0]);
    
    Cv = C(:);

    ind = 1:N; %nearest neighbor indices
    Lind = ind-m; 
    Rind = ind+m;
    Bind = ind+1;
    Tind = ind-1;

    %edge cases
    Lind(1:m) = (N-m+1):1:N;
    Rind(N-m+1:N) = [1:m];
    Bind([m:m:N]) = [1:m:(N-m+1)];
    Tind([1:m:(N-m+1)]) = [m:m:N];

    iter = 0;
    Utotal = 0;
    
    
    Sig_i = sum(Cv);
    for i = 1:N
        nn_sum = Cv(i)*sum([Cv(Lind(i)), Cv(Rind(i)), Cv(Tind(i)), Cv(Bind(i))]);
        %current energy 
        U_1 = -J*nn_sum - mu*B*Sig_i;
        Utotal = Utotal + U_1;
    end
    Utotal = Utotal/2; %divide by 2 to account for repeating
    
    while iter<=sweeps
        for i =1:N
            Sig_i = sum(Cv);
            nn_sum = Cv(i)*sum([Cv(Lind(i)), Cv(Rind(i)), Cv(Tind(i)), Cv(Bind(i))]);
            %current energy
            H_qj = -J*nn_sum - mu*B*Sig_i;

            %flip one spin
            Cv(i) = -Cv(i);
            %relcalculate energy
            Sig_i = sum(Cv);
            nn_sum = Cv(i)*sum([Cv(Lind(i)), Cv(Rind(i)), Cv(Tind(i)), Cv(Bind(i))]);
            H_trial = -J*nn_sum - mu*B*Sig_i;

            delta_E = H_trial - H_qj; %delta_E of one spin flip
            
            if delta_E>0
                if rand() > exp(-beta*delta_E)
                    %revert change
                    Cv(i) = - Cv(i);
                    H = H_qj;
                    %Utotal doesn't change
                else %keep change
                    H = H_trial;
                    Utotal = Utotal + delta_E;
                end
            else
                H = H_trial;
                Utotal = Utotal + delta_E;
            end
        end
        C = reshape(Cv, [m m]);
        if mod(iter,10) == 0 %update figures every 10th iteration
            figure(1);
            imagesc(C); colormap(gray); xlim([1 100]); ylim([1 100]); title(['Temperature (K) = ', num2str(T)])
            drawnow;
        end

        iter = iter+1;
%         disp(iter);
        U(iter) = Utotal;
        
        if mod(iter,10) ==0
            figure(2);
            plot(iter, U(iter)/N, 'o'); title(['Energy']);
            hold on;
            drawnow;
        end
    end
    Uout = mean(U(end-100:end)); %take average once U has converged
    L = (sum(Cv==1)-sum(Cv==-1))/N; %order: (N+ - N-)/N
    
    %In case spontaneous magnetized regions: calculate L_avg
    a = 1; 
    b = length(Cv);
    r = randi([a b],1,1000); %pick 1000 random small regions for local order
    L_avg = zeros(1, length(r));
    for i =1:length(r)
        Ct = [Cv(Lind(r(i))), Cv(Rind(r(i))), Cv(Tind(r(i))), Cv(Bind(r(i))), Cv(r(i))];
        Nt = length(Ct);
        L_avg(i) = (sum(Ct==1)-sum(Ct==-1))/Nt;
    end
    L_avg = mean(abs(L_avg));
    
    %figure(2); plot(U, '.'); title('Energy as function of iteration')
    %figure(3); imagesc(C); colormap(gray); title('final config');
end