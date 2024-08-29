function axchange(mode)
    fig = gcf;
    switch mode
        case 1
        case 2
            set(gca, 'XScale', 'log', 'YScale', 'log');
        case 3
            set(gca, 'XScale', 'log');
        case 4
            set(gca, 'YScale', 'log');
        otherwise
            disp('Invalid mode');
    end
end
