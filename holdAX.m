function holdAX(varargin)
            for n=2 : nargin
                if isa(varargin{n},'matlab.ui.control.UIAxes')
                    hold(varargin{n}, varargin{1})
                end
            end
end