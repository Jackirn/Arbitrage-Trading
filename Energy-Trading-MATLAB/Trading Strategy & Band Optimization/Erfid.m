function val = Erfid(x, y)
    f = @(t) exp(t.^2/2);
    val = sqrt(2/pi) * integral(f, y, x);
end
