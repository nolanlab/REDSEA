function TEXT = writeFCS(fname, DATA, TEXT, OTHER)
%writeFCS Save numeric data into FCS format used in flow cytometry.
%   writeFCS(FNAME, DATA) creates a ver 3.1 FCS file with file name FNAME
%       from the numerical matrix data. The program assings the parameters
%       by analyzing the properties of the numerical data.
%
%       FNAME               File name preferably with .fcs extension
%       DATA                NxM matrix with data. N matches the number of 
%                           events (cells) and M represents the number of
%                           channels. All negative values will be converted
%                           to zero. The datatype of the matric will
%                           determine the datatype in the file. If data is
%                           integer only, it will be saved as integer,
%                           otherwise single or double will be used,
%                           depending on the datatype.
%
%   writeFCS(FNAME, DATA, TEXT)  creates an FCS file with file name FNAME
%       from the numerical matrix data. The program takes parameters from
%       the struct TEXT or if they are missing, it guesses them from the
%       properties of the numerical data.
%
%       TEXT                is a struct with fields defined by the "Data 
%                           File Standard for Flow Cytometry Version FCS 
%                           3.1 Normative Reference". The fields can
%                           contain either text or numbers depending on
%                           their nature. If not all required fields are
%                           supplied, they will be guessed from the
%                           properties of the data.
%
%   writeFCS(FNAME, DATA, TEXT, OTHER)  creates an FCS file with file name 
%       FNAME from the numerical matrix data. The program takes parameters
%       from the struct TEXT or if they are missing, it guesses them from
%       the properties of the numerical data. OTHER is a struct that can
%       append any other custom data to the end of the FCS file
%
%       OTHER               is a struct with fields defined by the "Data 
%                           File Standard for Flow Cytometry Version FCS 
%                           3.1 Normative Reference". The fields can be
%                           chosen arbitrarily as long as they adhere to
%                           the format specification. 
%   Examples:
%
%       % This example creates an integer FCS file from random data points
%       DATA = round(horzcat(randn(1e4, 1)*20+120, randn(1e4, 1)*80+350));
%       writeFCS('integerFCS.fcs', DATA)
%
%       % This example creates an integer FCS file from random data points
%       DATA = single(horzcat(randn(1e4, 1)*20+120, randn(1e4, 1)*80+350));
%       writeFCS('singleFCS.fcs', DATA)
%
%
% Copyright 2013 Jakub Nedbal
% $Revision: 1.0 $  $Date: 2013/07/11 16:07:00 $ 

% check if TEXT is defined, else create a struct
if nargin < 3
    TEXT.BYTEORD = '1,2,3,4';
end

% If TEXT hasn't got BYTEORD specified, add it
if ~isfield(TEXT, 'BYTEORD')
    TEXT.BYTEORD = '1,2,3,4';
end

%% Parse the text, correct any errors if encountered and fix the DATA if
%  necessary according to the values in TEXT.
[DATA, TEXT] = parseTEXT(DATA, TEXT);


%% FCS standard version
fcsver = 3.1;


%% Convert filednames in TEXT into a consistent string
STR = fields2string(TEXT);

%% data offsets
Stext = 256;                            % Position of text start
Etext = Stext + numel(STR) - 1;         % Position of text end
Sdata = 2 ^ ceil(log2(Etext + 1));      % Position of data start
Edata = Sdata + size(DATA, 1) * sum(TEXT.PnB(:)) / 8 - 1;
                                        % Position of data end
Sgate = 0; Egate = 0;                   % Position of gates
if nargin > 3
    Sother = 256 * ceil((Edata + 1) / 256);
                                        % Position of other start
else
    Sother = [];                        % Leave empty if not in use
end

%% create text
outText = sprintf('FCS%3.01f    %8d%8d%8d%8d%8d%8d%8d', ...
                  fcsver, ...
                  Stext, ...
                  Etext, ...
                  Sdata, ...
                  Edata, ...
                  Sgate, ...
                  Egate, ...
                  Sother);
outText = sprintf('%s%s%s', outText, repmat(' ', 1, Stext - numel(outText)), STR);
outText = sprintf('%s%s', outText, repmat(' ', 1, Sdata - Etext - 1));

%% conform to the endianness of the output data
fields = fieldnames(TEXT);
if strcmp(TEXT.(fields{cellfun(@(x) strcmpi(x, 'byteord'), fields)}), ...
          '4,3,2,1')
    endian = true;
    machineformat = 'b';
else
    endian = false;
    machineformat = 'l';
end

%% Create data to match the data type and bites per number
datatype = fields{cellfun(@(x) strcmpi(x, 'datatype'), fields)};
switch TEXT.(datatype)
    case 'F'
        DT = 'single';
        DATAtmp = single(DATA');
    case 'D'
        DT = 'double';
        DATAtmp = double(DATA');
    case 'I'
        DT = 'uint8';
        % The DATA needs to be separated into 8-bit numbers.
        DATAtmp = zeros(sum(TEXT.PnB(:)) / 8, size(DATA, 1), 'uint8');
        u = 0;
        for i = 1 : size(DATA, 2)
            tmp = round(double(DATA(:, i)));
            % conform to endianness of data
            in = 1 : TEXT.PnB(i) / 8;       % big endian
            if endian
                in = fliplr(in);            % little endian
            end
            for j = 1 : numel(in)
                DATAtmp(u + in(j), :) = mod(floor(tmp / (256 ^ (j - 1))), 256);
            end
            u = u + j;
        end
       
end

fid = fopen(fname, 'w');
fprintf(fid, outText);
%ftell(fid)

fwrite(fid, DATAtmp, DT, 0, machineformat);
%fwrite(fid, DATA', 'uint16');
% if OTHER segment of the FCS file is specified
if nargin > 3
    outText = repmat(' ', 1, Sother - Edata - 1);
    STR = fields2string(OTHER);
    outText = sprintf('%s%s', outText, STR);
    fprintf(fid, outText);
end

fprintf(fid, '\n');
%ftell(fid)
fclose(fid);

end     % end of writeFCS function


function STR = fields2string(TEXT)
%% Convert fieldnames in TEXT into a consistent string
dlm = '/';      % delimiter

fields = fieldnames(TEXT);
STR = dlm;

for i = 1 : numel(fields)
    value = TEXT.(fields{i});
    field = upper(fields{i});
    if isnumeric(value) || iscell(value)
        for j = 1 : numel(value)
            f = field;
            if regexp(f, '^((PN)|(GN)|(RN))')
                f = [f(1), num2str(j), f(3 : end)];
            end
            if isnumeric(value)
                v = num2str(value(j));
            else
                v = value{j};
            end
            STR = sprintf('%s$%s%s%s%s', STR, f, dlm, v, dlm);
        end
    else
        STR = sprintf('%s$%s%s%s%s', STR, field, dlm, value, dlm);
    end
end

end     % End of fields2string function

function [DATA, TEXT] = parseTEXT(DATA, TEXT)

%% check if data is numeric
if ~isnumeric(DATA)
    error('DATA must be a numeric matrix');
end

%% check if text is a struct
if ~isstruct(TEXT)
    warning('TEXT must be a struct. It will be recreated automatically.')
    clear TEXT;
    TEXT.BYTEORD = '1,2,3,4';
end

%% Get the fieldnames of TEXT and process them
fns = fieldnames(TEXT);

% Process each required fieldname

%% $BYTEORD Byte order for data acquisition computer.
al = {'1,2,3,4', '4,3,2,1'};
fn = fns{cellfun(@(x) strcmpi(x, 'BYTEORD'), fns)};
if isempty(fn)
    TEXT.(fn) = '1,2,3,4';
    fprintf('Setting $BYTEORD to "1,2,3,4".\n');
else
    va = TEXT.(fn);
    if ~ischar(va)
        va = 'mistake';
    end
    if ~any(strcmpi(va, al))
        TEXT.(fn) = '1,2,3,4';
        warning('$BYTEORD must be either "1,2,3,4" or "4,3,2,1"');
    end
end


%% $DATATYPE Type of data in DATA segment (ASCII, integer, floating point).
al = {'I', 'F', 'D'}; nal = {'A'};
in = cellfun(@(x) strcmpi(x, 'DATATYPE'), fns);
if ~in
    fn = 'DATATYPE';
    TEXT = finddatatype(DATA, TEXT, fn);
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~ischar(va)
        va = 'mistake';
    end
    if ~any(strcmpi(va, al))
        warning(['$DATATYPE must be "I", "F" or "D"\nfor unsigned', ...
                 'integer, single or double float, respectively.']);
        if ~any(strcmpi(va, nal))
            warning('writeFCS does not support $DATATYPE "A".\n');
        end
        [DATA, TEXT] = finddatatype(DATA, TEXT, fn);
    end
end
TEXT.(fn) = upper(TEXT.(fn));


%% $MODE Data mode (list mode - preferred, histogram - deprecated).
al = {'L'}; nal = {'C', 'U'};
in = cellfun(@(x) strcmpi(x, 'MODE'), fns);
if ~in
    fn = 'MODE';
    TEXT.(fn) = 'L';
    fprintf('Setting $MODE to "L".\n');
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~ischar(va)
        va = 'mistake';
    end
    if ~any(strcmpi(va, al))
        warning('$MODE must be "L" for list mode.');
        if ~any(strcmpi(va, nal))
            warning('writeFCS does not support $MODE "C" or "U".\n');
        end
        TEXT.(fn) = 'L';
    end
end
TEXT.(fn) = upper(TEXT.(fn));


%% $NEXTDATA Byte offset to next data set in the file.
in = cellfun(@(x) strcmpi(x, 'NEXTDATA'), fns);
if ~in
    fn = 'NEXTDATA';
    TEXT.(fn) = 0;
    fprintf('Setting $NEXTDATA to 0.\n');
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = 0;
        warning('writeFCS only supports $NEXTDATA equal to 0 (zero).');
    elseif va ~= 0
        TEXT.(fn) = 0;
        warning('writeFCS only supports $NEXTDATA equal to 0 (zero).');
    end
end


%% $PAR Number of parameters in an event.
in = cellfun(@(x) strcmpi(x, 'PAR'), fns);
if ~in
    fn = 'PAR';
    TEXT.(fn) = size(DATA, 2);
    fprintf('Setting $PAR to %d.\n', size(DATA, 2));
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = size(DATA, 2);
        warning('writeFCS only supports numerical integer $PAR.');
        fprintf('Setting $PAR to %d.\n', size(DATA, 2));
    elseif va ~= size(DATA, 2);
        TEXT.(fn) = size(DATA, 2);
        warning('$PAR must be equal to the second dimension of DATA.');
        fprintf('Setting $PAR to %d.\n', size(DATA, 2));
    end
end
PAR = TEXT.(fn);


%% $TOT Total number of events in the data set.
in = cellfun(@(x) strcmpi(x, 'TOT'), fns);
if ~in
    fn = 'TOT';
    TEXT.(fn) = size(DATA, 1);
    fprintf('Setting $TOT to %d.\n', size(DATA, 1));
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = size(DATA, 1);
        warning('writeFCS only supports numerical integer $TOT.');
        fprintf('Setting $TOT to %d.\n', size(DATA, 1));
    elseif va ~= size(DATA, 1);
        TEXT.(fn) = size(DATA, 1);
        warning('$TOT must be equal to the first dimension of DATA.');
        fprintf('Setting $TOT to %d.\n', size(DATA, 1));
    end
end


%% $PnB Number of bits reserved for parameter number n.
in = cellfun(@(x) strcmpi(x, 'PnB'), fns);
if in
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = [];
        warning('writeFCS only supports numerical integer $PnB.');
    elseif numel(va) ~= PAR;
        TEXT.(fn) = [];
        warning(['Size of $PnB must be equal to the second ', ...
                 'dimension of DATA.']);
    end
end
TEXT = checkPnB(DATA, TEXT);


%% $PnE Amplification type for parameter n.
in = cellfun(@(x) strcmpi(x, 'PnE'), fns);
if in
    fn = fns{in};
    va = TEXT.(fn);
    if ~(iscell(va) || (ischar(va) && PAR == 1))
        TEXT.(fn) = {};
        warning('writeFCS only supports cell or character $PnE.');
    end
end
TEXT = checkPnE(DATA, TEXT);


%% $PnN Short name for parameter n.
in = cellfun(@(x) strcmpi(x, 'PnN'), fns);
if ~in
    fn = 'PnN';
    va = cellfun(@(x) sprintf('Ch %d', x), num2cell(1 : size(DATA, 2)), ...
                 'UniformOutput', false);
    err = true;
else
    fn = fns{in};
    va = TEXT.(fn);
    err = false;
end

% check if TEXT.PnN is defined correctly
if isempty(va) || (ischar(va) && size(DATA, 2) > 1) || ...
        (iscell(va) && size(DATA, 2) ~= numel(va))
    warning(['$PnE must be a cell. Its size must be equal to the ', ...
             'second dimension of DATA.']);
    va = cellfun(@(x) sprintf('Ch %d', x), num2cell(1 : size(DATA, 2)), ...
                 'UniformOutput', false);
    err = true;
end

% if TEXT.PnN is character, convert to a cell
if ischar(va)
    va = {va};
end

% Check that all values are a character, replace those which aren't
in = find(~cellfun(@ischar, va));
va(in) = cellfun(@(x) sprintf('Ch %d', x), num2cell(in), ...
                 'UniformOutput', false);
if ~isempty(in)
    warning(['$PnE must be a cell. Its size must be equal to the ', ...
             'second dimension of DATA.']);
    err = true;
end

% Make sure that commas are not used in the names, replace them with dashes
if any(~cellfun(@isempty, strfind(va, ',')))
    warning(['$PnE values must not include commas (,). Replacing by ', ...
         'dashes (-).']);
    err = true;
    va = cellfun(@(x) regexprep(x, ',', '-'), va, 'UniformOutput', false);
end

if any(err)
    txt = '';
    for i = 1 : numel(va)
        txt = sprintf('%s''%s'', ', txt, va{i});
    end
    fprintf('Setting $PnN to {%s}.\n', txt(1 : end - 2));
end
TEXT.(fn) = va;


%% $PnR Range for parameter number n.
in = cellfun(@(x) strcmpi(x, 'PnR'), fns);
if ~in
    fn = 'PnR';
    va = 2 .^ ceil(log2(max(DATA)));
    err = true;
else
    fn = fns{in};
    va = TEXT.(fn);
    err = false;
end

% check if TEXT.PnR is defined correctly
if isempty(va) || ~isnumeric(va)
    if ~isnumeric(va)
        warning(['$PnR must be a numeric array. Its size must be ', ...
                 'equal to the second dimension of DATA.']);
    end
    va = 2 .^ ceil(log2(max(DATA)));
    err = true;
end

% check if values of TEXT.PnR have integer values
in = va ~= round(va);
if any(in)
    err = true;
    va(in) = 2 .^ ceil(log2(max(DATA(in))));
    warning('$PnR must consist of integer values.');
end

% check if values of TEXT.PnR have positive values
in = va < 1;
if any(in)
    err = true;
    va(in) = 2 .^ ceil(log2(max(DATA(in))));
    warning('$PnR must consist of positive values.');
end


if err
    txt = sprintf('%d, ', va);
    fprintf('Setting $PnR to [%s].\n', txt(1 : end - 2));
end
TEXT.(fn) = va;

end     % end of parseTEXT function

function TEXT = finddatatype(DATA, TEXT, fn)
%% Work out the data type based on the properties of the DATA
% The functions checks whether data is integer, single or double. If single
% or double and consisting only of integer values, it converts to integer.

if isinteger(DATA)
    TEXT.(fn) = 'I';
    fprintf('Setting $DATATYPE to "I".\n');
elseif isfloat(DATA)
    % Check if all values are integers
    if ~any(DATA(:) - round(DATA(:)))
        TEXT.(fn) = 'I';
        fprintf('Setting $DATATYPE to "I".\n');
    else
        switch class(DATA)
            case 'single'
                TEXT.(fn) = 'F';
                fprintf('Setting $DATATYPE to "F".\n');
            case 'double'
                TEXT.(fn) = 'D';
                fprintf('Setting $DATATYPE to "D".\n');
            otherwise
                TEXT.(fn) = 'I';
                warning('Cannot determine $DATATYPE. Setting to "I".');
        end
    end
end

end     % end of finddatatype function

function TEXT = checkPnB(DATA, TEXT)
%% Check the values fo $PnB and make sure they comply with FCS format 
%  specification.
fns = fieldnames(TEXT);
in = cellfun(@(x) strcmpi(x, 'PnB'), fns);
if ~in
    fn = 'PnB';
    va = [];
    TEXT.(fn) = [];
    err = false;
else
    fn = fns{in};
    va = TEXT.(fn);
    err = true;
end

switch TEXT.(fns{cellfun(@(x) strcmpi(x, 'DATATYPE'), fns)})
    case 'I'        % integer
        uints = [8; 16; 32; 64];
        in = ~any(repmat(va(:)', numel(uints), 1) == ...
             repmat(uints, 1, numel(va)));
        if isempty(va)
            in = true(1, size(DATA, 2));
        elseif any(in) && err
            warning(['With integer $DATATYPE, $PnB can only be ', ...
                     '8, 16, 32 or 64']);
        end
        va(in) = 2 .^ sum(repmat(max(DATA(:, in)), 3, 1) > ...
                 repmat(2 .^ [8; 16; 32], 1, sum(in))) * 8;

    case 'F'        % single float
        in = ~any(va(:)' == 32 * ones(1, numel(va)));
        if (any(in) || isempty(in)) && err
            warning('With single float $DATATYPE, $PnB can only be 32.');
        end
        va = 32 * ones(1, size(DATA, 2));

    case 'D'        % double float
        in = ~any(va(:)' == 64 * ones(1, numel(va)));
        if (any(in) || isempty(in)) && err
            warning('With single float $DATATYPE, $PnB can only be 32.');
        end
        va = 64 * ones(1, size(DATA, 2));
end

if ~isequal(va, TEXT.(fn))
    txt = sprintf('%d, ', va);
    fprintf('Setting $PnB to [%s].\n', txt(1 : end - 2));
    TEXT.(fn) = va;
end

end     % end of checkPnB function

function TEXT = checkPnE(DATA, TEXT)
fns = fieldnames(TEXT);
in = cellfun(@(x) strcmpi(x, 'PnE'), fns);
if ~in
    fn = 'PnE';
    va = repmat({''}, 1, TEXT.(fns{cellfun(@(x) strcmpi(x, 'PAR'), fns)}));
    TEXT.(fn) = va;
    % error vector (zeros mean no error), last bit =0 means no warnings
    err = false(1, 5);
else
    fn = fns{in};
    va = TEXT.(fn);
    % error vector (zeros mean no error), last bit =1 means allow warnings
    err = [false(1, 4), true];
end

% check if it is a characted array with a single-dimensional data array
if ischar(va) && size(DATA, 2) == 1
    va = {va};
end

% check that it is a cell, else redefine it
if ~iscell(va)
    va = repmat({''}, 1, TEXT.(fns{cellfun(@(x) strcmpi(x, 'PAR'), fns)}));
    TEXT.(fn) = va;
    warning('$PnE must be a cell.')
end

switch TEXT.(fns{cellfun(@(x) strcmpi(x, 'DATATYPE'), fns)})
    case 'I'        % integer
        for i = 1 : size(DATA, 2)
            if numel(va) >= i
                if isempty(va{i})
                    cva = {[], []};
                else
                    cva = textscan(va{i}, '%f,%f');
                end
                if any(cellfun(@isempty, cva))
                    cva = {0, 0};
                    err(1) = true;
                end
            else
                cva = {0, 0};
                err(2) = true;
            end
            if any(cellfun(@(x) x < 0, cva))
                cva = {0, 0};
                err(3) = true;
            end
            if cva{1} > 0 && cva{2} == 0
                cva{2} = 10 .^ floor(log10(min(DATA(:, i))));
                err(3) = true;
            end
            va{i} = sprintf('%g,%g', cell2mat(cva));
        end
        if all(err([1 end]))
            warning(['$PnE must consist of two numbers separated by ', ...
                     'a comma such as "0,0" or "4.0,0.1".']);
        end
        if all(err([2 end]))
            warning(['Size of $PnE must be equal to the second ', ...
                     'dimension of DATA.']);
        end
        if all(err([3 end]))
            warning('$PnE must consist of non-negative numbers only.');
        end
        if all(err([4 end]))
            warning(['Unless being "0,0", $PnE must consist of ', ...
                     'positive numbers only.']);
        end

    case {'F', 'D'}     % single or double float
        if numel(va) ~= size(DATA, 2)
            err = true;
            warning(['Size of $PnE must be equal to the second ', ...
                     'dimension of DATA.']);
        end
        if ~iscell(va)
            warning('$PnE must be a cell');
            err = true;
        end
        try
            if any(~cellfun(@(x) isequal({0, 0}, textscan(x, '%f,%f')), va))
                err = true;
                warning(['With float $DATATYPE "F" or "D", $PnE must ', ...
                         'remain "0,0".']);
            end
        catch
            err(1) = true;
        end
        va = repmat({'0,0'}, 1, size(DATA, 2));
end

if any(err(1 : end - 1))
    txt = '';
    for i = 1 : numel(va)
        txt = sprintf('%s''%s'', ', txt, va{i});
    end
    fprintf('Setting $PnE to {%s}.\n', txt(1 : end - 2));
end
TEXT.(fn) = va;

end     % end of checkPnE function