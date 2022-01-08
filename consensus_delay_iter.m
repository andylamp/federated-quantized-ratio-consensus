function [params] = consensus_delay_iter(k, params)
%% Delay tolerant consensus iteration
%
  nodes = size(params.adj_matrix, 1);
  % transmission buffers for the nodes
  transm_z = zeros(nodes, nodes);
  transm_y = zeros(nodes, nodes);
  
  % get the indices that are greater than zero
  %z_idcs = find(z > 0);
  % now get the ceil for the y-states based on the predicate above
  %y_states(z_idcs) = ceil(y(z_idcs)./z(z_idcs));
  
  % get current y states of the nodes
  %y_states = y_states .* z0;
  
  % check if the diameter is right for us to check
  if mod(k, params.diameter*params.max_delay) == 1
    vote_states = params.y ./ params.z;
    params.max_votes = ceil(vote_states);
    params.min_votes = floor(vote_states); 
  end
  
  % get the nodes that need to be adjusted
  n_idcs = find(params.z > 1);
  
  
  % -- transmission process -- %
  
  % outgoing transmission
  for idc=1:length(n_idcs)
    out_edge = n_idcs(idc);
    
    % check if the process time is zero, only then execute
    if params.max_process_times(idc) == 0

      while params.z(out_edge) > 1
        % get a random node out of the available pool
        node = randi(nodes);
        
        % now adjust the transmitted weights
        if params.adj_matrix(node, out_edge) > 0
          flr = floor(params.y(out_edge) / params.z(out_edge));
          transm_y(node, out_edge) = transm_y(node, out_edge) + flr;
          transm_z(node, out_edge) = transm_z(node, out_edge) + 1;
          
          params.y(out_edge) = params.y(out_edge) - flr;
          params.z(out_edge) = params.z(out_edge) - 1;
        end
        
      end
    end
  end
        
  % sum with the incoming transmitted values
  params.y = params.y + sum(transm_y, 2);
  params.z = params.z + sum(transm_z, 2);
  % this is for delays
  idc = find(sum(transm_z, 2) > 0);
  params.max_process_times(idc) = randi(params.max_delay, size(idc, 1), 1);

  % sum incoming
  params.max_process_times = params.max_process_times - 1;
  
  % -- end transmission process -- %
  
  % find the min/max votes after values settle
  for h=1:nodes
    [~, cols] = find(params.adj_matrix(h, :) == 1);
    params.max_votes(h) = max(params.max_votes(h), max(params.max_votes(cols)));
    params.min_votes(h) = min(params.min_votes(h), min(params.min_votes(cols)));
  end
  
  % check if the termination condition is met
  if mod(k, params.diameter) == 0
    vote_index = params.max_votes - params.min_votes <= params.epsilon;
    nodes_converged = find(vote_index == 1);
    % nodes_diverged = find(vote_index ~= 0);
     
    params.node_stats(nodes_converged, 1) = min(params.node_stats(nodes_converged, 1), k);
    params.node_stats(nodes_converged, 2) = max(params.node_stats(nodes_converged, 2), k);
    %min(node_stats(:, 1))
    
    % check if we can terminate, which is when all votes are raised
    if vote_index
      params.can_stop = 1;
    end
  end
end