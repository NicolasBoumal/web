function results = comparison_plot()
% This script will run a number of methods on random synchronization
% tasks of increasing size, to compare running times and also record
% the singular values of the solutions found. You will need CVX,
% SeDuMi and SDPLR installed for the competitors, as well as Manopt
% for the staircase method (the proposed method).
% See www.manopt.org.
%
% Nicolas Boumal, UCLouvain, May 18, 2014.

    clear all;
    close all;
    clc;

    d = 3;
    sigma = .3;
    cvx_solver sdpt3;
    
    % dry run to get all codes into Matlab's memory.
    h = helper(10, d, sigma);  %#ok<NASGU>
    
    mm = union(round(linspace(10, 2000, 21)), round(logspace(log10(10), log10(2000), 11)));
    iters = length(mm);
    results = cell(iters, 1);
    
    for iter = 1 : iters
        results{iter} = helper(mm(iter), d, sigma);
    end
    
    figure;
    hold on;
    
    colors = { [ 27,158,119]/255, [217, 95,  2]/255, [117,112,179]/255, ...
               [231, 41,138]/255, [102,166, 30]/255, [230,171,  2]/255, ...
               [166,118, 29]/255, [102,102,102]/255, [77,  0,  75]/255, ...
               [  8, 64,129]/255 };

    
    fields = fieldnames(results{1});
    for k = 1 : numel(fields)
        field = fields{k};
        time = zeros(size(mm));
        for iter = 1 : iters
            res = results{iter};
            resf = res.(field);
            time(iter) = resf.time;
        end
        plot(mm, time, '.-', 'color', colors{k}, 'MarkerSize', 15, 'LineWidth', 2);
    end
    
    hold off;
    
    legend(fields);
    
end

function results = helper(m, d, sigma)

    results = struct();
    
    n = m*d;
    Rtrue = randortho(d, m);
    Rtrue_stacked = reshape(multitransp(Rtrue), [d, n])';
    W = sigma*randn(n);
    W = (W+W')/2;
    C = Rtrue_stacked * Rtrue_stacked' + W;

    sdpt3 = struct();
    if m <= 300
        t = tic;
        X = linearcost_cvx(-C/(n*m), m, d);
        sdpt3.time = toc(t);
        sdpt3.svd = svd(X);
    else
        sdpt3.time = NaN;
        sdpt3.svd = NaN(m, 1);
    end
    results.sdpt3 = sdpt3;

    sedumi = struct();
    if m <= 300
        t = tic;
        X = linearcost_sedumi(-C/(n*m), m, d);
        sedumi.time = toc(t);
        sedumi.svd = svd(X);
    else
        sedumi.time = NaN;
        sedumi.svd = NaN(m, 1);
    end
    results.sedumi = sedumi;

    sdplr_ = struct();
    if m <= 1000
        t = tic;
        Y = linearcost_sdplr(-C/(n*m), m, d);
        sdplr_.time = toc(t);
        sdplr_.svd = svd(Y);
    else
        sdplr_.time = NaN;
        sdplr_.svd = NaN(m, 1);
    end
    results.sdplr = sdplr_;

    sdplr_forcerank = struct();
    t = tic;
    Y = linearcost_sdplr(-C/(n*m), m, d, d+1);
    sdplr_forcerank.time = toc(t);
    sdplr_forcerank.svd = svd(Y);
    results.sdplr_forcerank = sdplr_forcerank;

    staircase = struct();
    t = tic;
    Y = linearcost_staircase(-C/(n*m), m, d);
    staircase.time = toc(t);
    staircase.svd = svd(Y);
    results.staircase = staircase;

    eigmethod = struct();
    t = tic;
    [V, D] = eigs(C, d);
    round2orthogonal(V*sqrt(D), d);
    eigmethod.time = toc(t);
    results.eig = eigmethod;
    
end
