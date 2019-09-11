function [sol, its, hist] = myNewton(f, df, x0, TOL, maxit)
its = 0;
x0 = x0';
hist = x0;
x1 = x0;
for its = 1:maxit
    func_vals = num2cell(x0);
    if norm( f(func_vals{:}) )<TOL
        break
    end
    g = f(func_vals{:});
    H = df(func_vals{:});
    delta_x = -(H^-1)*g;
    x1 = x0 + delta_x;
    hist = [hist; x1];
    x0 = x1;
end
hist = hist;
sol = x1;
end