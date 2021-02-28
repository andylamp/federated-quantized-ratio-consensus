%% Federated quantized capacity ratio consensus /w finite time termination
%
% This is the main stub that performs the experiments in our paper. We test
% across a variety of parameters; concretely we define a parameter tuple as
% the node count, delay and run it across a number of trials - in the paper
% we used 10 trials.
%
% For more details, please see either the README.md or our paper which can
% be found here: https://arxiv.org/abs/xxxx.xxxx
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
max_iter = 100;
% the graph connectivity target
graph_connectivity = 0.5;
% the number of nodes
nodes_to_test = [20, 50]; % , 100, 150, 200, 400, 600
% nodes for "large scale" testing
% nodes_to_test = [20, 200, 500, 1000, 5000, 10000];
nodes_to_test_len = length(nodes_to_test);
% set the node length as a convenience variable
node_len_array = 1:nodes_to_test_len;

% trials for regular testing
trials = 2;
% trials for large scale testing
% trials = 5;
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

% setup variables
params.type = "quant-normal";   % normal async
params.pflag = 0;               % enable printing
params = setup_vars(params);    % setup environment variables

for t=trials_arr
    fprintf("\n** Running trial %d\n", t);
    
    for n=nodes_to_test_len
      nodes = nodes_to_test(n);
      fprintf(" -- Running for node size: %d\n", nodes);
      
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
      
      % our termination flag
      can_terminate = 0;
      
      % start ticking the global one
      global_tic = tic;
      
      % now run for all the iterations
      for k=1:max_iter
        
        % transmission buffers for the nodes
        transm_z = zeros(nodes, nodes);
        transm_y = zeros(nodes, nodes);
        
        % get the indices that are greater than zero
        z_idcs = find(z > 0);
        % now get the ceil for the y-states based on the predicate above
        y_states(z_idcs) = ceil(y(z_idcs)./z(z_idcs));
        
        % get current y states of the nodes
        y_states = y_states .* z0;
        
        % check if the diameter is right for us to check
        if mod(k, diameter) == 1
          vote_states = y ./ z;
          max_votes = ceil(vote_states);
          min_votes = floor(vote_states); 
        end
        
        % get the nodes that need to be adjusted
        n_idcs = find(z > 1);
        
        
        % -- transmission process -- %
        
        % outgoing transmission
        for idc=1:length(n_idcs)
          out_edge = n_idcs(idc);
          while z(out_edge) > 1
            % get a random node out of the available pool
            node = randi(nodes);
            
            % now adjust the transmitted weights
            if adjMatrix(node, out_edge) > 0
              flr = floor(y(out_edge) / z(out_edge));
              transm_y(node, out_edge) = transm_y(node, out_edge) + flr;
              transm_z(node, out_edge) = transm_z(node, out_edge) + 1;
              
              y(out_edge) = y(out_edge) - flr;
              z(out_edge) = z(out_edge) - 1;
            end
            
          end
          
        end
        
        % sum with the incoming transmitted values
        y = y + sum(transm_y, 2);
        z = z + sum(transm_z, 2);
        
        % -- end transmission process -- %
        
        % find the min/max votes after values settle
        [v_row, v_col] = find(adjMatrix == 1);
        max_votes = max(max_votes(v_row), max_votes(v_col));
        min_votes = min(min_votes(v_row), min_votes(v_col));
        
        % check if the termination condition is met
        if mod(k, diameter) == 0
          if max_votes - min_votes <= epsilon
            can_terminate = 1;
          end
        end
        
        % check if we can terminate
        if can_terminate == 1
          break
        end
        
      end
      
      % record the time it took to converge
      total_time_global(n, t) = toc(global_tic);
      
      fprintf("\t== DEBUG INFO: \n\t Converged at iteration:  %d (ouf ot max: %d),\n\t Initial avg:  %d, \n\t Final avg:    %d,\n\t Time elapsed: %d seconds.\n", ...
        k, max_iter, init_avg, sum_y0/sum_z0, total_time_global(n, t));
      if can_terminate == 0
        fprintf("\t^^ ERROR: Exhausted max iterations (max=%d) for converging -- stopping\n", k);
      else
        fprintf("\t^^ INFO: Converged after %d out of %d iterations for %d nodes\n", k, max_iter, nodes);
      end
      
      fprintf(" -- Finished for node size: %d\n", nodes);
    end
    fprintf("** Finished trial %d\n", t);
end

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
%legend("secs");

% print the figure 
st = sprintf("nmax_%d_trials_%d_total_exec_time", ...
  nodes_to_test(end), trials);
print_fig(fig, st, params);

% -- finished plotting execution time -- %




