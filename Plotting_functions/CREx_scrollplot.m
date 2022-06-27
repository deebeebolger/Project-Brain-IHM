function CREx_scrollplot( dx, x, ax )
    % Set appropriate axis limits
    for ii = 1:numel(ax)
        set( ax(ii), 'xlim', [0 dx] );
    end

    % Create Uicontrol slider
    % The callback is another local function, this gives us more
    % flexibility than a character array.
    uicontrol('style','slider',...
        'units', 'normalized', 'position', [0.1 0.01 0.8 0.05],...
        'callback', @(slider, ~) CREx_scrollcallback( ax, dx, slider ), ...
        'min', 0, 'max', max(x)-dx );
end