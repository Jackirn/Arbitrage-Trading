function y = erfi(x)
    y = (2 / sqrt(pi)) * arrayfun(@(z) integral(@(t) exp(t.^2), 0, z, 'ArrayValued', true), x);
end