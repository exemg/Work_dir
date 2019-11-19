function errarrout = avg_err(errarrin,dim)
if ~isnumeric(errarrin)
    error('avg_err accepts only numeric arrays');
elseif ~any(size(errarrin))
    error('nonpositive array dimension found in avg_err array');
elseif prod(dim)==1
    warning('avg_err used on single value array. Verify the code does not add excessive complexity');
end

errarrout = sqrt(sum(errarrin.^2,dim))/size(errarrin,dim);
end