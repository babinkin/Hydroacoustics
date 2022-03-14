function semilogxOneNOctaveBar(y_center,f_low, f_high,width, color,faceAlpha)
    % ------------
    if gca == 0
        figure;
    end
    % ------------
    washold = ishold;
    if ~washold
        hold on;
    end
    % ------------
    for i = 1:length(y_center)
        a = 0.5*(f_low(i)*(1 + width) + f_high(i)*(1 - width) );
        b = 0.5*(f_low(i)*(1 - width) + f_high(i)*(1 + width) );
        rectangle('Position', [a 0 (b - a) y_center(i)]);
        bh = fill([a b b a], [0 0 y_center(i) y_center(i)], color);
        bh.FaceAlpha = faceAlpha;
    end
    % ------------
    if ~washold
        hold off;
    end
    % ------------
    set(gca,'XScale','log');
end