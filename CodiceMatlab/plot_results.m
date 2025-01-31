function [] = plot_results(f, x0, xseq, lowbound, upbound)
    
    [X, Y] = meshgrid(linspace(lowbound, upbound, 200), linspace(lowbound, upbound, 200));
    Z = 100*(Y - X.^2).^2 + (ones(size(X))-X).^2;     %Rosenbrock

    %PLOT-------------------------------------------------------
    fig = figure()
    low = lowbound-1;
    up = upbound+1;
    seq_contour = [low:0.135:up; low:0.135:up];
    contour(X, Y, Z, cellfun(f, num2cell(seq_contour, 1)));
    hold on
    plot([x0(1) xseq(1, :)], [x0(2) xseq(2, :)] , '--*')
    hold off
    grid on
   
    %fig2 = figure()
    %colormap(fig2, hot)
    %mesh(X,Y,Z)

    
end