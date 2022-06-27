function CREx_scrollcallback( ax, dx, slider, varargin )
    % Scroller callback loops through the axes objects and updates the xlim
    val = slider.Value;
    for ii = 1:numel(ax)
        set( ax(ii), 'xlim', val + [0, dx] );
    end
end