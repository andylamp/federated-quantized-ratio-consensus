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
rng("default");

%largeconstant=10^10                %%% WORK ALSO FOR RING
connectivity=0.5;                   
nodes = 20;
max_iter = 100;

% epsilon to stop
epsilon = 1;

% % generate a square matrix for the adjacency matrix.
% adjMatrix = rand(nodes);
% % set the 1's first
% adjMatrix(adjMatrix > connectivity) = 1;
% % now set the zeroes
% adjMatrix(adjMatrix <= connectivity) = 0;
% % now set the diagonal values
% adjMatrix(1:nodes+1:end) = 1:nodes;
% 
% diameter=1;
% AAA = adjMatrix;
% while (find(AAA==0)>0)
%     AAA = AAA*adjMatrix;
%     diameter=diameter+1;
% end

% rng('default');


Ad = zeros(nodes);
for i=1:nodes
    for j=1:nodes
        if rand()>connectivity
            Ad(i,j)=1;              % AD - Adjacency Martix
        end
    end
end


for i=1:nodes
    Ad(i,i)=1;
end


diam=1;
AAA = Ad;
while (find(AAA==0)>0)
    AAA = AAA*Ad;
    diam=diam+1;
end


quant_step = 1000;   % we can multiply by max of z0
y0 = quant_step * randi(100, nodes, 1);
z0 = 2 * ones(nodes, 1);
for j=1:nodes
    if (mod(j,2)==1)
        z0(j) = 100;
    else
        z0(j) = 300;
    end
end

% display the vectors
y=y0
z=z0

sum(y)
sum(z)

adjMatrix = Ad;

% quant_step = 1000;   % we can multiply by max of z0
% y0 = quant_step * randi(100, nodes, 1);
% z0 = 2 * ones(nodes, 1);
% % even cap
% z0(2:2:end) = 100;
% % odd cap
% z0(1:2:end) = 300;
% 
% % copy the vectors
% y = y0;
% z = z0;

% % display the vectors
% y=y0
% z=z0




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   - RING
% Ad = zeros(nodes,nodes);
% 
% for i=1:nodes-1
%     Ad(i+1,i)=1;
%     Ad(i,i)=1;
% end
% Ad(1,nodes)=1;
% 
% y0 = [4; 8; 10; 3; 5; 1];
% z0 = [1; 1; 1; 1; 1; 1];
% y=y0;
% z=z0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Ad = [0 1 1 0 0 0 0; 1 0 0 1 0 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1; 1 1 1 1 0 0 0; 0 0 0 0 1 0 1; 0 0 0 0 1 0 0];   % CDC paper plot
%  nodes = 7;
%  y0 = [15; 5; 11; 4; 3; 13; 9]
%  z0 = [1; 1; 1; 1; 1; 1; 1]
%  y=y0;
%  z=z0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:nodes
%     Ad(i,i)=1;
% end


% 
% transm_z = zeros(nodes);
% transm_y = zeros(nodes);

initial_aver = sum(y0)/sum(z0);

first = zeros(max_iter, nodes);
first2 = zeros(max_iter, nodes);

y_states = y0;
% nodes_states_z = z0;

node_stats = zeros(nodes, 2);
node_stats(:, 1) = max_iter;


y
z
initial_aver

% first_max = zeros(max_iter, 1);
% first_min = zeros(max_iter, 1);

max_plot = zeros(max_iter, 1);
min_plot = zeros(max_iter, 1);
% 
% nonz = zeros(1,nodes);


for k=1:max_iter
  % transmission buffers for the nodes
  transm_z = zeros(nodes, nodes);
  transm_y = zeros(nodes, nodes);

  % get the indices that are greater than zero
  z_idcs = find(z > 0);
  % now get the ceil for the y-states based on the predicate above
  y_states(z_idcs) = ceil(y(z_idcs)./z(z_idcs));
  
  % update the first2
  first2(k, :) = y_states';
  
  % get current y states of the nodes
  y_states = y_states .* z0;
  
  % update the first
  first(k, :) = y_states';
  
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
  for h=1:nodes
    [~, cols] = find(adjMatrix(h, :) == 1);
    max_votes(h) = max(max_votes(h), max(max_votes(cols)));
    min_votes(h) = min(min_votes(h), min(min_votes(cols)));
  end

  % check if the termination condition is met
  if mod(k, diameter) == 0
    pad_idc = (k/diameter);
    
    min_idc = (pad_idc-1)*diameter + 1;
    max_idc = (pad_idc)*diameter;
    max_plot(min_idc:max_idc) = max_votes(1);
    min_plot(min_idc:max_idc) = min_votes(1);
    
    vote_index = max_votes - min_votes <= epsilon;
    nodes_converged = find(vote_index == 1);
    nodes_diverged = find(vote_index ~= 0);

    node_stats(nodes_converged, 1) = min(node_stats(nodes_converged, 1), k);
    node_stats(nodes_converged, 2) = max(node_stats(nodes_converged, 2), k);
    %min(node_stats(:, 1))

    % check if we can terminate, which is when all votes are raised
    if vote_index
      break;
    end
  end
end

y
z

% normalise the values
idc = find(z > 0);
y_states(idc) = ceil(y(idc)./z(idc));

%     for j=1:nodes
%         if (z(j) >0)
%             y_states1(j) = ceil(y(j)/z(j));
%         end
%     end

% figure
% hold
% title('Quantized Average Consensus via Mass Equalization - Envelope','FontWeight','Normal')
% plot(first_max)
% plot(first_min)
% hold

%max_iter


%first = (first'*z0)'

figure
hold
stairs(ceil(first(1:k, :)/quant_step))
title('Load per Node According to Processing Capacity', 'FontWeight','Normal')
ylabel('Node State Variables ($q_j^s[k]$)', 'interpreter', "latex", "fontsize", 18)
xlabel('Number of Iterations ($k$)', 'interpreter', "latex", "fontsize", 18)

C_average = y_states

initial_aver 


figure
hold
stairs(first2(1:k, :))
title('Load per Processing Cycle + Max/Min plots (Dashed)')
ylabel('Node State Variables ($q_j^s[k]$)', 'interpreter', "latex", "fontsize", 18)
xlabel('Number of Iterations ($k$)', 'interpreter', "latex", "fontsize", 18)
stairs(max_plot(1:k), '--')
stairs(min_plot(1:k), '--')


