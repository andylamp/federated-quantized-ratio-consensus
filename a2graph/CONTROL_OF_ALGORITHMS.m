close all
clear all
clc

connectivity=0.5;                   
nodes=20;
iter=100;

% MAXIMUM PROCESSS TIMES
Max_Process_B = 5; 
processs_times = randi(Max_Process_B,nodes,1);


Ad = zeros(nodes,nodes);

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


quant_step=100;   % we can multiply by max of z0
y0=10 + quant_step*randi(90,nodes,1);
% z0=2*ones(nodes,1);
% for j=1:nodes
%     if (mod(j,2)==1)
%         z0(j) = 10;
%     else
%         z0(j) = 30;
%     end
% end
z0 = 10 + randi(90,nodes,1);
y=y0
z=z0

y0 = y0.*z0;



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


%first1 = zeros(iter,nodes);

%first2 = zeros(iter,nodes);



first1 = Quantized_Mass_Splitting_NO_OSCILLATIONS_STOPPING_S_V1_FED(diam,nodes,Ad,y,z,y0,z0,iter); 

first2 = Quantized_Mass_Splitting_NO_OSCILLATIONS_STOPPING_S_V1_PR_FED(diam,nodes,Ad,y,z,y0,z0,iter,quant_step,Max_Process_B,processs_times); 



figure
hold
subplot(2,1,1)
%stairs(ceil(first1/quant_step))
stairs(first1)
title('(A)','Interpreter','latex')
ylabel('Node State Variables ($q_j^s$[k])','Interpreter','latex')
xlabel('Number of Iterations (k)','Interpreter','latex')


subplot(2,1,2)
%stairs(ceil(first2/quant_step))
stairs(first2)
title('(B)','Interpreter','latex')
ylabel('Node State Variables ($q_j^s[k]$)','Interpreter','latex')
xlabel('Number of Iterations (k)','Interpreter','latex')


initial_aver = sum(y0)/sum(z0);

initial_aver 


