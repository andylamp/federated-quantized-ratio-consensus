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
node_to_test_len = length(nodes_to_test);
% set the node length as a convenience variable
node_len_array = 1:node_to_test_len;

% trials for regular testing
trials = 2;
% trials for large scale testing
% trials = 5;
trials_arr = 1:trials;

% the quantisation step to use
quantization_step = 100;

for t=trials_arr
    fprintf("\n** Running trial %d\n", t);
    
    for n=node_to_test_len
      nodes = nodes_to_test(n);
      fprintf(" -- Running for node size: %d", nodes);
      
      % generate the capacity for z0
      z0 = gen_workload(nodes);
      % now generate the capacity for y0
      y0_opt.multiplier = quantization_step;
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
      
      fprintf(" == DEBUG INFO: Initial Average %d, sum(y0): %d, sum(z0): %d", ...
        init_avg, sum_y0, sum_z0);
      
      % firstly, lets generate the graph
      [~, diameter, nodes, adjMatrix] = gen_graph(nodes);
      
      % transmit the nodes
      transm_z = zeros(nodes);
      transm_y = zeros(nodes);
      
      % the node states
      z_states = z;
      y_states = y;
      
      % variables to hold the max/nim
      first_max = NaN(max_iter, 1);
      first_min = NaN(max_iter, 1);
      
      % now run for all the iterations
      for k=1:max_iter
        % get the indices that are greater than zero
        z_idcs = find(z > 0);
        % now get the ceil for the y-states based on the predicate above
        y_states(z_idcs) = ceil(y(z_idcs)./z(z_idcs));
      end
      
      fprintf(" -- Finished for node size: %d", nodes);
    end
    fprintf("** Finished trial %d\n", t);
end

