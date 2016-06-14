function strOut = cifStrSplit(strIn)
% Split string respecting apostrophes; i.e., first split at
% apostrophes, then split the text outside apostrophes at
% whitespaces. E.g., string """x1 'x2 x3 x4' x5 x6 'x7 x8' x9"""
% is split into the cell array {'x1', 'x2 x3 x4', x5, x6, 'x7 x8', x9}.
% Efficiency is sacrificed for readability; regex masters
% are welcome to rewrite this into a 1-liner using their heavenly
% powers.
    
    
    % If there are not any apostrophes in the string, split at
    % whitespaces
    if length(strfind(strIn, '''')) == 0
        strOut = strsplit(strtrim(strIn)); % Split at whitespaces and be done
    else
        inFlag = false; % Determines whether we are inside the
                        % apostrophes
        strOut = {}; % Store substrings
        isIn = [];   % Remember whether substring was inside or
                     % outside
        iStart = 1; % Mark start and end of each substring
        for iChar = 1:length(strIn)
            if strIn(iChar) == ''''
                if inFlag == false % We were outside and now
                                   % passing to the inside
                    inFlag = true;
                    iEnd = iChar-1; % outside string ended
                    if (iEnd ~= 0)
                        strOut = {strOut{:} strIn(iStart:iEnd)};
                        inOrOut = [isIn false];
                    end
                    iStart = iChar+1; %inside string starting
                end
                
                if inFlag == true % We were inside and now passing
                                  % to the outside
                    inFlag = false;
                    iEnd = iChar - 1; % inside string ended
                    if (iEnd + 1 ~= length(strIn))
                        strOut = [strOut{:} strIn(iStart:iEnd)];
                        inOrOut = [isIn true];
                    end
                    iStart = iChar+1;
                end
            end
        end
        
        % Split all outside substrings at whitespaces
        for iSubstr = 1:length(strOut)
            if isIn(iSubstr) == false
                strOut{iSubstr} = strsplit(strtrim(strout{iSubstr}));
            end
        end
        strOut = flatten(strOut);
    end
end