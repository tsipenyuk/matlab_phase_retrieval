function cifData = cif2mat(readFile)
% cif2mat - read a CIF file and convert it to Matlab structure
%
% Synopsis ([]s are optional)
%   cifData = cif2mat(readFile)
%
% Caveat
%   Does not comply the formal requierements of the *.cif standard
%   (http://www.iucr.org/resources/cif/spec/version1.1). In
%   particular, does not support non-simple data values. (Data
%   values are split by whitespace. Adjust the 'strsplit' function in
%   the code to repair this behaviour). Use with caution, adjust to
%   your needs. See detailed parser description below. 
%
% Inputs ([]s are optional)
%   (string)     readFile   *.cif file to be parsed
%
% Outputs
%   (cell array) cifData    cell array containing all the
%                           information stored in the file in the
%                           form of structs. 
%
% Examples
%    Download http://files.rcsb.org/download/2OLX-sf.cif to your
%    Matlab folder. Then, call:
%    """
%    >> cifData = cif2mat('2olx-sf.cif')
%    """
%    Use the autocompletion feature or call 'fieldnames(cifData)'
%    to access the contents of the file.
%
%
% Description
%   Parses a CIF file to the Matlab cell array, converting numbers
%   to doubles where possible. To be more precise, this parser does 
%   the following:
%     1) Find lines in the *.cif file that start with the character
%   '#'. Split the file into corresponding blocks. (Each block is a
%   cell array containing lines between two consecutive #-starting
%   lines.) 
%     2) If the second line in the block (line immediately after
%   the #-line) does NOT start with the word 'loop_':
%         i) Create a struct named after the beginning of the
%         line. If the line has the form '_big.apple 33', the
%         first struct will be called 'big'. (Name starts at the
%         2nd line character and ends at the line character
%         preceeding the '.').
%
%         ii) Create a substructure with the name starting after
%         the '.' character and ending with any whitespace
%         character. If the line has the form '_big.apple 33', the
%         substructure will be called 'apple' (field name) and 
%         whatever follows the whitespace character will be
%         assigned as its value. The parser will attempt to use
%         str2num to convert the value to a number.
%
%         (Remark. Structure names are made valid using
%         matlab.lang.makeValidName. The program yields valid
%         results only if the lines in such blocks contain a '.',
%         a whitespace character, and a value. Otherwise, results
%         are unpredictable.)
%
%     3) If the second line in the block (line immediately after
%   the #-line) DOES start with the word 'loop_':
%        i) The parser counts all lines within the block that
%        start with the character '_'. The resulting number
%        numCols is declared as the number of columns in the table
%        stored in the block.
%
%        ii) For numCols lines following the 'loop_'-line, the
%        parser creates a structure name and a substructure name 
%        (From 2nd character in the line until the '.'-character
%        and from the '.'-character until the end of the line,
%        respectively). For example, for the line '_big.apple' a
%        structure 'big' (only one for the whole block' and a
%        substructure 'apple' is created. This substructure is
%        the cell array that contains the corresponding column
%        from the table.
%        
%        iii) For all lines following the numCols lines, the parser
%        splits each line at whitespaces and writes the first
%        numCols entries into the corresponding substructures.
%        CAVEAT. The parser does not respect the *cif-convention
%        that everything between apostrophes counts as a single
%        table entry. Example. The block
%        """
%        #
%        loop_
%        _symmetry_equiv.id
%        _symmetry_equiv.pos_as_xyz
%         1  'X,  Y,  Z'
%         2  '-X+1/2,  -Y,  Z+1/2'
%         3  'X+1/2,  -Y+1/2,  -Z'
%         #
%        """
%        will be parsed to cifData with entries
%        """
%        >> cifData.symmetry_equiv.id
%        
%        ans =
%        
%             1     2     3
%        
%        >> cifData.symmetry_equiv.pos_as_xyz
%         
%        ans = 
%        
%            ''X,'    ''-X+1/2,'    ''X+1/2,'
%             
%        """
%        instead of being parsed to 
%        """
%        >> cifData.symmetry_equiv.pos_as_xyz
%         
%        ans = 
%        
%            'X,  Y,  Z'  '-X+1/2,  -Y,  Z+1/2'  'X+1/2,  -Y+1/2,  -Z'
%        """
%        In other words, the apostrophes structure is not respected.
%
%        iv) The parser looks at the first 10 rows of the data (or at
%        all the data if smaller than 10 rows are provided). For each
%        given column, if at least one value V in the first 10 rows
%        has the 'true' status when '[res, status] = str2num(V)' is
%        called, all the values in this column are converted to
%        doubles using str2double.
%
%  
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   This code is a part of the Matlab Phase Retrieval Sandbox
%   software. The corresponding license may be found at
%   https://github.com/atsipenyuk/phase_retrieval/blob/master/LICENSE 
%
% Changes
%   2016-06-14  First Edition
    tic;
    
    cifData = struct;
    cifData.type = 'cifData';
    %cifData.outfile = readFile;
    
    % initialize file
    fileID = fopen(readFile, 'r');
    rawText = fread(fileID, inf, '*char');
    
    % parse lines by end-of-lines
    splitLines = strread(rawText, '%s', 'delimiter', '\n');
    numLines = length(splitLines);
    
    % and determine lines starting with '#'
    blockBounds = [];
    for iLine = 1:numLines
        if splitLines{iLine}(1) == '#'
            blockBounds = [blockBounds iLine];
        end
    end 
    
    % parse blocks
    blockData = {};
    for iBlock = 1:length(blockBounds)-1
        blockData{iBlock} = parseBlock({splitLines{blockBounds(iBlock): ...
                                              blockBounds(iBlock + 1)}});
        blockName = matlab.lang.makeValidName(cell2mat(...
            fieldnames(blockData{iBlock})));
        cifData.(blockName) = blockData{iBlock}.(blockName);
    end
    toc;
end


function blockData = parseBlock(blockLines)
% parseBlock - parses a block of .cif-file
%
% Input
%   (cell array) blockLines --- cell array of lines containing exactly one
%                               block of a *.cif file. In other words,
%                               the first and the last elements of
%                               blockLines are lines starting with '#'
%                               (other content of these lines is ignored).
% Output
%   (struct) blockData --- struct containing fields with values
%                          specified in the cif file
    if and(blockLines{1}(1) ~= '#', blockLines{end}(1) ~= '#')
        error(['mprs -- cif2mat: internal error #1. Please ' ...
               'report this incident.']);
    end


    % Distinguish between tables (arrays) and scalars
    if all(blockLines{2}(1:5) == 'loop_')
        % blockName is the string that is repeated before the
        % dot in a block of the *cif file. It does not contain the first '_'
        % character and the dot at the end.
        blockName = matlab.lang.makeValidName(blockLines{3}(2:strfind(blockLines{3}, '.')-1));
        blockData = struct(blockName, struct);
        numLines = length(blockLines);
    
        % Count how many columns are in the table
        numCols = 0;
        for iLine = 2:1:(numLines - 1)
            if blockLines{iLine}(1) == '_'
                numCols = numCols + 1;
            end
        end
        
        % number of rows = total number of lines - 2 ('#' lines) - numCols
        % -'loop' line
        numRows = numLines - 3 - numCols;
        if numRows <= 0
            error(['mprs - cif2mat: Internal error #2. Please report ' ...
                   'this incident.']);
        end
        
        
        % Declare a struct for each column. The name of the struct
        % is the string between '.' and ' ' (or EOL). 
        structName = {};
        for iCol = 3:1:(numCols + 2) % +2 because we start with 3rd
                                     % line (after '#' and 'loop_')
            % Determine the name
            dotPos = strfind(blockLines{iCol}, '.');
            blankPos = strfind(blockLines{iCol}, ' ');
            if blankPos ~= [] % Whitespace is found
                structName{iCol-2} = ...
                    matlab.lang.makeValidName(blockLines{iCol}(dotPos+1:blankPos-1));
            else
                structName{iCol-2} = ...
                    matlab.lang ...
                    .makeValidName(blockLines{iCol}(dotPos+1:end));
            end
            
            blockData.(blockName).(structName{iCol-2}) = {};
        end
        
        % Loop over the table and write down the contents
        for iRow = 3 + numCols:1:2 + numCols + numRows
            % Split at whitespace. This function must be adjusted
            % to be consistent with *CIF format
            currentLine = strsplit(blockLines{iRow}); 
                
            for iCol = 1:1:numCols
                blockData.(blockName).(structName{iCol}){iRow - numCols - 2} = ...
                    currentLine{iCol};
            end
        end

        % Also, check whether the entries of the column are numbers
        % or not (try to convert any of the first 10 values)
        structIsNum = zeros([1,numCols]);
        for iRow = 1:1:min(numRows,10)
            for iCol = 1:1:numCols
                [~, status] = ...
                    str2num(blockData.(blockName).(structName{iCol}){iRow});
                if status == true
                    structIsNum(iCol) = true;
                end
            end
        end

        % Convert the entries that have passed the test above
        for iCol = 1:1:numCols
            if structIsNum(iCol) == true
                blockData.(blockName).(structName{iCol}) = ...
                    str2double(blockData.(blockName).(structName{iCol}));
            end
        end
       
    else % scalar entries\
        % blockName is the string that is repeated before the
        % dot in a block of the *cif file. It does not contain the first '_'
        % character and the dot at the end.
        blockName = matlab.lang.makeValidName(blockLines{2}(2:strfind(blockLines{2}, '.')-1));
        blockData = struct(blockName, struct);
        numLines = length(blockLines);
        
        for iLine = 2:1:(numLines - 1)
            dotPos = strfind(blockLines{iLine}, '.');
            blankPos = strfind(blockLines{iLine}, ' ');
            structName = matlab.lang.makeValidName(blockLines{iLine}(dotPos+1:blankPos-1));
            blockData.(blockName).(structName) = blockLines{iLine}(blankPos+1:end);

            % Convert to number if possible
            [res, status] = str2num(blockData.(blockName).(structName));
            if status == true
                blockData.(blockName).(blockLines{iLine}(dotPos+1: ...
                                                  blankPos-1)) = res;
            end
        end
    end
end