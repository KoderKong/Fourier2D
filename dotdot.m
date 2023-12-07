function dotdot(varargin)
persistent count
persistent width
switch nargin
    case 1
        assert(~isempty(count) && ~isempty(width))
        if varargin{1}
            count = count+1;
            fprintf('.')
            if count >= width
                count = 0;
                fprintf('\n')
            end
        else
            if count > 0
                count = 0;
                fprintf('\n')
            end
        end
    case 2
        fprintf('%s\n',varargin{1})
        count = 0;
        width = varargin{2};
end
