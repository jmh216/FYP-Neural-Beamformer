function t = is_instring(longstr,shortstr)
%is_instring - determines whether shortstr is in longstr
%
% Syntax: t = is_instring(longstr,shortstr)
%
% Wei Xue, Imperial College, 20180712

if isempty(strfind(longstr,shortstr))
    t = false;
else
    t = true;    
end