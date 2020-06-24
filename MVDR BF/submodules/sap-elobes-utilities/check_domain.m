function [pp,fieldVal] = check_domain(pp,fieldName,fieldVal)
% check_domain(p,fieldName,fieldVal) 
% Check wether a struct contains a certain field name, and set a default value to it.
% 
% Wei Xue, Imperial College, July 13, 2018

if ~isfield(pp,fieldName)
    pp.(fieldName) = fieldVal;
else
    fieldVal = pp.(fieldName);
end