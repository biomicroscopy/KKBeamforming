function [pS] = convertParamsToSingle(p)


fields = fieldnames(p);
pS = p;
for i = 1:length(fields)
    if (isa(pS.(fields{i}),'double'))
        pS.(fields{i}) = single(pS.(fields{i}));
    end
end


end