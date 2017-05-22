function Jortberg_Fourier(n_max)

    x= -5 : 0.0001 : 5;
    y= 1-x;
    %n_max=1000;

    FTcos=0.5;
    FTsin=0;
    for n = 1:n_max
        An=2*(1-(-1)^n)/(n^2*pi^2);
        FTcos = FTcos + An*cos(n*pi*x);
        Bn=2/(n*pi);
        FTsin=FTsin + Bn*sin(n*pi*x);
    end
    figure; plot(x, FTsin, 'g'); hold on, plot(x, FTcos, 'r'); hold on, plot(x,y, 'b');
    xlabel('x'); ylabel('y'); title([num2str(n_max), ' terms of Cosine and Sine Fourier Expansion']);
    legend('Sine Expansion', 'Cosine Expansion', 'Expanded function, y=1-x');

end