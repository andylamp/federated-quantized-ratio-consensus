%% Federated quantized capacity ratio consensus /w finite time termination
%
% This is the main stub that performs the experiments in our paper. We test
% across a variety of parameters; concretely we define a parameter tuple as
% the node count, delay and run it across a number of trials - in the paper
% we used 10 trials.
%
% For more details, please see either the README.md or our paper which can
% be found here: https://arxiv.org/abs/2104.03126
%
% Authors:
%
%  - Andreas A. Grammenos (ag926@cl.cam.ac.uk)
%  - Apostolos Rikos (rikos@kth.se)
%  - Themistoklis Charalambous (themistoklis.charalambous@aalto.fi)
%
% License: GPLv3
%

%% Initialisation
%
clear; clc; close all;

% for reproducibility, fix the seed
rng("default")

% the maximum iterations to go through
max_iter = 2000;
% attempts to try with different graphs
max_attempts = 10;
% the graph connectivity target
graph_connectivity = 0.6;
% the number of nodes
nodes_to_test = [50, 100, 200, 300, 600, 1000, 2000, 5000, 10000]; % , 100, 150, 200, 400, 600
% nodes for "large scale" testing
% nodes_to_test = [20, 200, 500, 1000, 5000, 10000];
nodes_to_test_len = length(nodes_to_test);
% set the node length as a convenience variable
node_len_array = 1:nodes_to_test_len;

% trials 50 regular testing
trials = 5;
% trials for large scale testing
% trials = 50;
trials_arr = 1:trials;

% epsilon to stop
epsilon = 1;

% the quantisation step to use
quantisation_step = 100;

% converge statistics
cov_min_global = zeros(nodes_to_test_len, trials);
cov_max_global = zeros(nodes_to_test_len, trials);
cov_mean_global = zeros(nodes_to_test_len, trials);
cov_win_global = zeros(nodes_to_test_len, trials);

% execution time
total_time_global = zeros(nodes_to_test_len, trials);
total_trial_time = zeros(trials, 1);

% setup variables
params.type = "quant-normal";   % normal async
params.pflag = 0;               % enable printing
params = setup_vars(params);    % setup environment variables

for t=trials_arr
    fprintf("\n** Running trial %d\n", t);
    trial_tic = tic;
    for n=node_len_array
      completed = 0;
      attempt = 1;
      while completed == 0
        nodes = nodes_to_test(n);
        fprintf(" -- Running for node size: %d (attempt %d)\n", ...
          nodes, attempt);
        
        % generate the capacity for z0
        z0 = gen_workload(nodes);
        % now generate the capacity for y0
        y0_opt.multiplier = quantisation_step;
        y0_opt.distribution_type = "randomised";
        y0 = gen_workload(nodes, y0_opt);
        sum_y0 = sum(y0);
        sum_z0 = sum(z0);
        
        % set the initialisation values
        z = z0;
        y = y0;
        init_avg = sum_y0 / sum_z0;
        
        % sanity check to ensure that the invariant holds
        assert(sum(y0) > sum(z0));
        
        fprintf("\t== DEBUG INFO: \n\t Initial Average: %d,\n\t Sum(y0): %d,\n\t Sum(z0): %d.\n", ...
          init_avg, sum_y0, sum_z0);
        
        % firstly, lets generate the graph
        [~, diameter, nodes, adjMatrix] = gen_graph(nodes);
        
        % the node states
        z_states = z;
        y_states = y;    
  
        % plot variables
        max_plot = NaN(max_iter, 1);
        min_plot = NaN(max_iter, 1);
        
        node_stats(:, 1) = max_iter;
        
        % start ticking the global one
        global_tic = tic;

        % assign parameters
        con_params.z = z;
        con_params.y = y;
        con_params.can_stop = 0;
        con_params.epsilon = epsilon;
        con_params.diameter = diameter;
        con_params.adj_matrix = adjMatrix;
        

        % node converge statistics, the structure is as follows:
        %
        % row count: nodes
        %
        % description of each column:
        %
        % 1. min iteration that satisfied the constraint
        % 2. max iteration that satisfied the constraint (and then not
        % changed)
        % 3. number of flips between min, max (not used currently)
        %
        con_params.node_stats = zeros(nodes, 2);
        con_params.node_stats(:, 1) = max_iter;
        
        % now run for all the iterations
        k = 1;
        while con_params.can_stop == 0 && k < max_iter
          % perform the consensus iteration
          [con_params] = consensus_iter(k, con_params);
          % increment the counter
          k = k + 1;
        end
        
        % compute the stats
        cov_min = min(con_params.node_stats(:, 1));        % min converge.
        cov_max = max(con_params.node_stats(:, 2));        % max converge.
        cov_mean = mean(con_params.node_stats(:, 1));      % mean converge for max
        
        % update the variables.
        cov_min_global(n, t) = cov_min;
        cov_max_global(n, t) = cov_max;
        cov_mean_global(n, t) = cov_mean;
        cov_win_global(n, t) = cov_max - cov_min; 

        % assert(cov_win_global(n, t) >= 0)      

        % record the time it took to converge
        total_time_global(n, t) = toc(global_tic);
        
        fprintf("\t== DEBUG INFO: \n\t Converged at iteration:  %d (ouf ot max: %d),\n\t Initial avg:  %d, \n\t Final avg:    %d,\n\t Time elapsed: %d seconds.\n", ...
          k, max_iter, init_avg, sum_y0/sum_z0, total_time_global(n, t));
        if k >= max_iter
          fprintf("\t^^ ERROR: Exhausted max iterations (max=%d) for converging -- stopping\n", k);
        else
          fprintf("\t^^ INFO: Converged after %d out of %d iterations for %d nodes\n", k, max_iter, nodes);
        end

        % check what happened
        if cov_win_global(n, t) < 0 && attempt < max_attempts
          attempt = attempt + 1;
        elseif attempt >= max_attempts
          error("Failed to converge and max attempts (%d) exhausted", ...
            max_attempts)
        else
          % notify we completed
          completed = 1;
        end

        % while end
      end
      
      fprintf(" -- Finished for node size: %d (attempts: %d)\n", ...
        nodes, attempt);
    end
    trial_toc = toc(trial_tic);
    fprintf("** Finished trial %d, elapsed time: %d seconds\n", t, trial_toc);
    total_trial_time(t) = trial_toc;
end

fprintf("== Finished running %d trials, total time elapsed: %d seconds\n", ...
  t, sum(total_trial_time));

% -- plot execution time -- %

fig = figure;
hold on; box on;
  cur = total_time_global;
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);
hold off;

ylabel("Time (s)", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
xlabel("Nodes", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
title("Total time to converge", ...
  "Interpreter", "Latex", "FontName", "Times New Roman");
xticks(node_len_array)
xticklabels(num2cell(nodes_to_test))
legend("Seconds");

% print the figure 
st = sprintf("nmax_%d_trials_%d_total_exec_time", ...
  nodes_to_test(end), trials);
print_fig(fig, st, params);

% -- finished plotting execution time -- %

% -- plot convergence -- %

msg = "Converge iterations";
% create the figure
fig = figure;

% enable hold and boxing
hold on; box on;
  cur = cov_max_global;
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);
hold off;

% handle figure presentation
ylabel("Iterations", "FontName", "Times New Roman");
xlabel("Nodes", "Interpreter", "Latex", "FontName", "Times New Roman");
st = sprintf(msg);
title(st, "Interpreter", "Latex", "FontName", "Times New Roman");
xticks(node_len_array)
xticklabels(num2cell(nodes_to_test))
legend("Iterations");

% print the figure 
st = sprintf("nmax_%d_trials_%d_converge_iterations", ...
  nodes_to_test(end), trials);
print_fig(fig, st, params);

% -- finished convergence plot -- %

% -- plot convergence statistics -- %

wlegs = ["min", "max", "mean", "window"];
msg = "Converge statistics";
% create the figure
fig = figure;

% enable hold and boxing
hold on; box on;

  % min cov
  cur = cov_min_global;
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

  % max cov
  cur = cov_max_global;
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

  % mean
  cur = cov_mean_global;
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

  % window cov
  cur = cov_win_global;
  err_std = std(cur, 0, 2)' / sqrt(trials);
  errorbar(node_len_array, mean(cur, 2)', err_std, "LineWidth", 2);

hold off;

% handle figure presentation
ylabel("Iterations", "FontName", "Times New Roman");
xlabel("Nodes", "Interpreter", "Latex", "FontName", "Times New Roman");
st = sprintf(msg);
title(st, "Interpreter", "Latex", "FontName", "Times New Roman");
xticks(node_len_array)
xticklabels(num2cell(nodes_to_test))
legend(wlegs);

% print the figure 
st = sprintf("nmax_%d_trials_%d_converge_statistics", ...
  nodes_to_test(end), trials);
print_fig(fig, st, params);

% -- finished plotting convergence statistics -- %




