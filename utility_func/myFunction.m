function myFunction(varargin)
    p = inputParser;

    % Define the expected parameters and their default values
    addParameter(p, 'myKey', '', @(x) ischar(x));
    addParameter(p, 'param1', 0, @(x) isnumeric(x));
    addParameter(p, 'param2', false, @(x) islogical(x));

    % Parse the input using name-value pairs
    parse(p, varargin{:});

    % Access the parsed parameters using p.Results
    myKeyValue = p.Results.myKey;
    param1Value = p.Results.param1;
    param2Value = p.Results.param2;

    fprintf('myKey: %s\n', myKeyValue);
    fprintf('param1: %d\n', param1Value);
    fprintf('param2: %d\n', param2Value);
end