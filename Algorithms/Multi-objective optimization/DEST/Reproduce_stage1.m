function [Population] = Reproduce_stage1(Problem,Population,W,type,Z)
% The environmental selection of SPEA2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

   delta = 0.9;
   nr = 2;
   N = size(W,1);
   T = ceil(N/10);
   B = pdist2(W,W);
   [~,B] = sort(B,2);
   B = B(:,1:T);
   for i = 1 : N
        % Choose the parents
        if rand < delta
              P = B(i,randperm(size(B,2)));
        else
              P = randperm(N);
        end

        % Generate an offspring
        if type == 1          
           Offspring = OperatorGAhalf(Problem,Population(P(1:2)));
        else 
           Offspring  = OperatorDE(Problem,Population(i),Population(P(1)),Population(P(2)));
        end
                
        theta = rad2deg(acos(1-pdist2((Population(P).objs-repmat(Z,length(P),1)),W(P,:),'cosine')));
        g_old = [] ;
        Dominate = true(1,size(P,2));
        for j = 1 :size(P,2)
            g_old(j) = theta(j,j);
            k = any(Population(P(j)).objs<Offspring.objs) - any(Population(P(j)).objs>Offspring.objs);
                if k == 1
                    Dominate(j) = false;
                end
        end    
        g_new = rad2deg(acos(1-pdist2((Offspring.objs-Z),W(P,:),'cosine')));

        Population(P(find((g_old>=g_new)&Dominate,nr))) = Offspring;
     
    end
   
   
end

