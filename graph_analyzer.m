%%
% <https://www.mathworks.com MathWorks> 

function results = graph_analyzer(f, xrange)
%GRAPH_ANALYZER  Automated calculus-based graph analyzer.
%   graph_analyzer(f) analyses the symbolic function f(x): domain, vertical
%   asymptotes, intercepts, horizontal/slant asymptotes, critical points,
%   monotonicity, local extrema, concavity and inflection points, prints a
%   report, and draws everything on one annotated graph.
%
%   f can be a symbolic expression, a string such as '(x^2-4)/(x^2-1)', or a
%   function handle such as @(x) x.^3 - 3*x. With no input, it asks you to
%   type f(x).
%
%   graph_analyzer(f, [a b]) uses [a b] as the plotting window and as the
%   window for numeric root finding (default search window is [-20 20]).
%
%   results = graph_analyzer(f) also returns the findings in a struct.
%
%   Requires the Symbolic Math Toolbox. Designed for polynomial and rational
%   functions (the standard curve-sketching cases).
%
%   Examples:
%       syms x
%       graph_analyzer(x^3 - 3*x);
%       graph_analyzer((x^2 - 4)/(x^2 - 1));
%       graph_analyzer((x^2 + 1)/x);

    syms x
    if nargin < 1 || isempty(f)
        f = input('Enter f(x) using x, e.g. (x^2-4)/(x^2-1) : ', 's');
    end
    f = toSym(f, x);
    userWindow = nargin >= 2 && ~isempty(xrange);
    if userWindow
        win = xrange;
    else
        win = [-20 20];
    end
    tol = 1e-8;

    fprintf('\n=========== CALCULUS GRAPH ANALYSIS ===========\n');
    fprintf('f(x) = %s\n', char(f));

    %% STEP 1: Domain and vertical asymptotes -----------------------------
    [~, den] = numden(f);
    cand = realRoots(den, x, win);            % points where denominator = 0
    dpts = sort(double(cand(:))');            % excluded from the domain
    VA = [];
    for k = 1:numel(cand)
        if isVerticalAsymptote(f, x, cand(k))
            VA(end+1) = double(cand(k)); %#ok<AGROW>
        end
    end
    VA = sort(VA);

    if isempty(dpts)
        fprintf('\n[1] Domain: all real numbers\n');
    else
        fprintf('\n[1] Domain: all real x except x = %s\n', numlist(dpts));
    end
    fprintf('    Vertical asymptotes: x = %s\n', numlist(VA));

    %% STEP 2: Intercepts --------------------------------------------------
    yint = [];
    if ~any(abs(dpts) < tol)
        y0 = double(subs(f, x, 0));
        if isfinite(y0) && isreal(y0)
            yint = [0 y0];
        end
    end
    numer = numden(f);
    xi = double(realRoots(numer, x, win));
    xi = excludePts(sort(xi(:)'), dpts);

    fprintf('\n[2] Intercepts\n');
    if isempty(yint)
        fprintf('    y-intercept: none\n');
    else
        fprintf('    y-intercept: (0, %.4g)\n', yint(2));
    end
    fprintf('    x-intercepts: x = %s\n', numlist(xi));

    %% STEP 3: Horizontal and slant asymptotes (end behaviour) -------------
    HA = zeros(0, 2);     % rows: [direction  y-value]
    SA = zeros(0, 3);     % rows: [direction  slope  intercept]
    for s = [1 -1]
        L  = tryLimit(f, x, s*Inf);
        Ld = double(L);
        if isfinite(Ld)
            HA(end+1, :) = [s Ld]; %#ok<AGROW>
        elseif isinf(Ld)
            m  = tryLimit(f/x, x, s*Inf);        % slope = lim f(x)/x
            md = double(m);
            if isfinite(md) && abs(md) > tol
                b  = tryLimit(f - m*x, x, s*Inf); % intercept = lim (f - m x)
                bd = double(b);
                if isfinite(bd)
                    SA(end+1, :) = [s md bd]; %#ok<AGROW>
                end
            end
        end
    end

    fprintf('\n[3] End behaviour\n');
    if isempty(HA) && isempty(SA)
        fprintf('    No horizontal or slant asymptotes\n');
    end
    for k = 1:size(HA, 1)
        fprintf('    Horizontal asymptote: y = %.4g  (x -> %s)\n', HA(k, 2), dirName(HA(k, 1)));
    end
    for k = 1:size(SA, 1)
        fprintf('    Slant asymptote: y = %.4g*x + %.4g  (x -> %s)\n', SA(k, 2), SA(k, 3), dirName(SA(k, 1)));
    end

    %% STEP 4: First derivative -> monotonicity and extrema ----------------
    f1 = simplify(diff(f, x));
    [n1, d1] = numden(f1);
    c1 = double(realRoots(n1, x, win));       % f'(x) = 0
    c2 = double(realRoots(d1, x, win));       % f'(x) undefined
    crit = sort([c1(:)' c2(:)']);
    crit = excludePts(crit, dpts);            % must lie in the domain of f

    bp1 = uniquetol(sort([crit dpts]), tol);
    [e1, s1] = signChart(f1, x, bp1);

    maxima = zeros(0, 2);
    minima = zeros(0, 2);
    for c = crit
        j  = find(abs(bp1 - c) < 1e-6, 1);
        sl = s1(j);                           % sign of f' just left of c
        sr = s1(j + 1);                       % sign of f' just right of c
        yc = double(subs(f, x, c));
        if sl > 0 && sr < 0
            maxima(end+1, :) = [c yc]; %#ok<AGROW>
        elseif sl < 0 && sr > 0
            minima(end+1, :) = [c yc]; %#ok<AGROW>
        end
    end

    fprintf('\n[4] First derivative: f''(x) = %s\n', char(f1));
    fprintf('    Critical points: x = %s\n', numlist(crit));
    fprintf('    Monotonicity:\n');
    for k = 1:numel(s1)
        fprintf('      %-24s %s\n', intervalStr(e1(k), e1(k+1)), ...
            trendWord(s1(k), 'increasing', 'decreasing', 'constant'));
    end
    for k = 1:size(maxima, 1)
        fprintf('    Local maximum at (%.4g, %.4g)\n', maxima(k, 1), maxima(k, 2));
    end
    for k = 1:size(minima, 1)
        fprintf('    Local minimum at (%.4g, %.4g)\n', minima(k, 1), minima(k, 2));
    end
    if isempty(maxima) && isempty(minima)
        fprintf('    No local extrema\n');
    end

    %% STEP 5: Second derivative -> concavity and inflection points --------
    f2 = simplify(diff(f1, x));
    n2 = numden(f2);
    ipc = double(realRoots(n2, x, win));      % f''(x) = 0
    ipc = excludePts(sort(ipc(:)'), dpts);

    bp2 = uniquetol(sort([ipc dpts]), tol);
    [e2, s2] = signChart(f2, x, bp2);

    infl = zeros(0, 2);
    for c = ipc
        j = find(abs(bp2 - c) < 1e-6, 1);
        if s2(j) * s2(j + 1) < 0              % concavity changes sign
            infl(end+1, :) = [c double(subs(f, x, c))]; %#ok<AGROW>
        end
    end

    fprintf('\n[5] Second derivative: f''''(x) = %s\n', char(f2));
    fprintf('    Concavity:\n');
    for k = 1:numel(s2)
        fprintf('      %-24s %s\n', intervalStr(e2(k), e2(k+1)), ...
            trendWord(s2(k), 'concave up', 'concave down', 'linear'));
    end
    for k = 1:size(infl, 1)
        fprintf('    Inflection point at (%.4g, %.4g)\n', infl(k, 1), infl(k, 2));
    end
    if isempty(infl)
        fprintf('    No inflection points\n');
    end

    %% STEP 6: Plot everything on one graph ---------------------------------
    keys = [xi crit ipc dpts 0];
    if userWindow
        pl = win(1);  ph = win(2);
    else
        pl = min(keys) - 2;  ph = max(keys) + 2;
        if ph - pl < 8
            c0 = (pl + ph)/2;  pl = c0 - 4;  ph = c0 + 4;
        end
    end

    g  = matlabFunction(f, 'Vars', x);
    xs = linspace(pl, ph, 6000);
    ys = evalVec(g, xs);
    ys(abs(imag(ys)) > 1e-12) = NaN;
    ys = real(ys);
    ys(abs(ys) > 1e6) = NaN;
    for a = VA                                % break the curve at each pole
        idx = find(xs(1:end-1) < a & xs(2:end) >= a);
        ys(idx) = NaN;
    end

    ykeys = [0 yint(2:end) maxima(:, 2)' minima(:, 2)' infl(:, 2)' HA(:, 2)'];
    ylo = min(ykeys);  yhi = max(ykeys);
    span = max(yhi - ylo, 2);
    ylo = ylo - 0.75*span;  yhi = yhi + 0.75*span;

    figure('Name', 'Calculus Graph Analyzer', 'Color', 'w');
    hold on; grid on; box on;
    hdl = gobjects(0);
    leg = {};

    hdl(end+1) = plot(xs, ys, 'b-', 'LineWidth', 1.8);
    leg{end+1} = 'f(x)';
    plot([pl ph], [0 0], 'k-', 'HandleVisibility', 'off');
    plot([0 0], [ylo yhi], 'k-', 'HandleVisibility', 'off');

    for a = VA
        hv = plot([a a], [ylo yhi], 'r--', 'LineWidth', 1.2);
    end
    if ~isempty(VA)
        hdl(end+1) = hv;
        leg{end+1} = 'vertical asymptote';
    end
    for k = 1:size(HA, 1)
        hh = plot([pl ph], HA(k, 2)*[1 1], 'm--', 'LineWidth', 1.2);
    end
    if ~isempty(HA)
        hdl(end+1) = hh;
        leg{end+1} = 'horizontal asymptote';
    end
    for k = 1:size(SA, 1)
        hs = plot([pl ph], SA(k, 2)*[pl ph] + SA(k, 3), 'g--', 'LineWidth', 1.2);
    end
    if ~isempty(SA)
        hdl(end+1) = hs;
        leg{end+1} = 'slant asymptote';
    end

    if ~isempty(xi)
        hdl(end+1) = plot(xi, zeros(size(xi)), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
        leg{end+1} = 'x-intercept';
        labelPts([xi(:) zeros(numel(xi), 1)], 'k');
    end
    if ~isempty(yint)
        hdl(end+1) = plot(yint(1), yint(2), 'kd', 'MarkerFaceColor', 'c', 'MarkerSize', 8);
        leg{end+1} = 'y-intercept';
        labelPts(yint, 'k');
    end
    if ~isempty(maxima)
        hdl(end+1) = plot(maxima(:, 1), maxima(:, 2), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
        leg{end+1} = 'local max';
        labelPts(maxima, 'r');
    end
    if ~isempty(minima)
        hdl(end+1) = plot(minima(:, 1), minima(:, 2), 'bv', 'MarkerFaceColor', 'b', 'MarkerSize', 9);
        leg{end+1} = 'local min';
        labelPts(minima, 'b');
    end
    if ~isempty(infl)
        hdl(end+1) = plot(infl(:, 1), infl(:, 2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
        leg{end+1} = 'inflection point';
        labelPts(infl, [0 0.5 0]);
    end

    xlim([pl ph]);
    ylim([ylo yhi]);
    xlabel('x'); ylabel('y');
    title(['f(x) = ' char(f)], 'Interpreter', 'none');
    legend(hdl, leg, 'Location', 'best');

    results = struct( ...
        'domainExclusions',     dpts, ...
        'verticalAsymptotes',   VA, ...
        'horizontalAsymptotes', HA, ...
        'slantAsymptotes',      SA, ...
        'xIntercepts',          xi, ...
        'yIntercept',           yint, ...
        'criticalPoints',       crit, ...
        'localMax',             maxima, ...
        'localMin',             minima, ...
        'inflectionPoints',     infl);
end


%% ======================= HELPER FUNCTIONS ================================

function f = toSym(f, x)
% Accepts a symbolic expression, a function handle, or text such as
% '(x^2-4)/(x^2-1)' and returns a symbolic expression in x.
    if isa(f, 'function_handle')
        f = f(x);
    elseif ischar(f) || isstring(f)
        f = char(f);
        try
            f = str2sym(f);
        catch
            f = sym(f);
        end
    end
    f = sym(f);
    v = symvar(f);
    for k = 1:numel(v)
        if ~strcmp(char(v(k)), 'x')
            error('graph_analyzer:badInput', ...
                'f must be a function of x only (found variable "%s").', char(v(k)));
        end
    end
end

function r = realRoots(expr, x, win)
% Real roots of expr = 0 as a symbolic column vector.
% Exact solve() first; then a numeric scan of the window adds any roots
% solve() missed (e.g. periodic roots of trig functions).
    r = sym(zeros(0, 1));
    if isAlways(expr == 0)
        return
    end
    w = warning('off', 'all');
    try
        s = solve(expr == 0, x, 'Real', true);
        s = s(:);
        if ~isempty(s) && isempty(symvar(s))
            r = s;
        end
    catch
    end
    warning(w);

    nr = numericRoots(expr, x, win);
    have = double(r);
    for k = 1:numel(nr)
        if isempty(have) || all(abs(have - nr(k)) > 1e-6)
            r = [r; sym(nr(k))]; %#ok<AGROW>
            have(end+1, 1) = nr(k); %#ok<AGROW>
        end
    end
end

function y = evalVec(g, xs)
% Evaluate g on the row vector xs and always return a row vector the same
% size as xs. Falls back to point-by-point evaluation if g is not vectorised
% or returns a constant.
    ok = false;
    try
        y = g(xs);
        if isscalar(y)
            y = y * ones(size(xs));
            ok = true;
        elseif isequal(size(y), size(xs))
            ok = true;
        end
    catch
    end
    if ~ok
        y = zeros(size(xs));
        for k = 1:numel(xs)
            try
                y(k) = g(xs(k));
            catch
                y(k) = NaN;
            end
        end
    end
end

function r = numericRoots(expr, x, win)
% Sign-change scan + fzero over the search window.
    r = [];
    g  = matlabFunction(expr, 'Vars', x);
    xs = linspace(win(1), win(2), 4001);
    ys = evalVec(g, xs);
    ys(~isfinite(ys) | abs(imag(ys)) > 0) = NaN;
    for k = 1:numel(xs) - 1
        if ys(k) == 0
            r(end+1) = xs(k); %#ok<AGROW>
        elseif ys(k) * ys(k+1) < 0
            try
                z = fzero(g, [xs(k) xs(k+1)]);
                if abs(g(z)) < 1e-6          % reject sign changes at poles
                    r(end+1) = z; %#ok<AGROW>
                end
            catch
            end
        end
    end
    if ~isempty(r)
        r = uniquetol(sort(r), 1e-6);
    end
end

function tf = isVerticalAsymptote(f, x, a)
% True if f blows up to +/-Inf from at least one side of x = a.
    tf = false;
    for side = {'left', 'right'}
        L = tryLimit(f, x, a, side{1});
        if isinf(double(L))
            tf = true;
            return
        end
    end
end

function L = tryLimit(expr, x, a, varargin)
% limit() that returns sym(NaN) instead of erroring when no limit exists.
    try
        L = limit(expr, x, a, varargin{:});
        d = double(L);
        if ~(isscalar(d) && isreal(d))
            L = sym(NaN);
        end
    catch
        L = sym(NaN);
    end
end

function [edges, sgn] = signChart(expr, x, bp)
% Sign of expr on every interval cut out by the breakpoints bp.
    g = matlabFunction(expr, 'Vars', x);
    edges = [-Inf, bp(:)', Inf];
    sgn = zeros(1, numel(edges) - 1);
    for k = 1:numel(edges) - 1
        a = edges(k);  b = edges(k+1);
        if isinf(a) && isinf(b)
            t = 0;
        elseif isinf(a)
            t = b - 1;
        elseif isinf(b)
            t = a + 1;
        else
            t = (a + b)/2;
        end
        v = g(t);
        if ~isreal(v) || ~isfinite(v)
            sgn(k) = NaN;
        elseif abs(v) < 1e-10
            sgn(k) = 0;
        else
            sgn(k) = sign(v);
        end
    end
end

function pts = excludePts(pts, bad)
% Remove from pts any value that is (numerically) in bad.
    keep = true(size(pts));
    for k = 1:numel(pts)
        if any(abs(pts(k) - bad) < 1e-6)
            keep(k) = false;
        end
    end
    pts = pts(keep);
end

function s = numlist(v)
    if isempty(v)
        s = 'none';
        return
    end
    s = strjoin(arrayfun(@(t) sprintf('%.4g', t), v, 'UniformOutput', false), ', ');
end

function s = intervalStr(a, b)
    s = sprintf('(%s, %s)', boundStr(a), boundStr(b));
end

function t = boundStr(v)
    if isinf(v) && v < 0
        t = '-Inf';
    elseif isinf(v)
        t = 'Inf';
    else
        t = sprintf('%.4g', v);
    end
end

function w = trendWord(s, pos, neg, zero)
    if isnan(s)
        w = 'undefined';
    elseif s > 0
        w = pos;
    elseif s < 0
        w = neg;
    else
        w = zero;
    end
end

function n = dirName(s)
    if s > 0
        n = '+Inf';
    else
        n = '-Inf';
    end
end

function labelPts(P, clr)
    for k = 1:size(P, 1)
        text(P(k, 1), P(k, 2), sprintf('  (%.3g, %.3g)', P(k, 1), P(k, 2)), ...
            'Color', clr, 'FontSize', 9, 'VerticalAlignment', 'bottom');
    end
end
