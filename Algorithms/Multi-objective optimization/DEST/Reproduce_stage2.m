function [Population] = Reproduce_stage2(Problem,Population,W,type,Z)
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
   T = ceil(N/20);
   B = pdist2(W,W);
   [~,B] = sort(B,2);
   B = B(:,1:T);
   Angle_P = acos(1-pdist2((Population.objs-repmat(Z,N,1)),W,'cosine'));
   [~,associate_P] = min(Angle_P,[],2);
   CV = sum(max(0,Population.cons),2);
  
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
        CV_O = sum(max(0,Offspring.cons),2);
        Angle_O = acos(1-pdist2((Offspring.objs-Z),W,'cosine'));
        [~,associate_O] = min(Angle_O,[],2);
        

        
        count = 0;
        for k = 1 : size(P,2) %当前选择P(k)向量进行比较 
            y = any(Offspring.objs<Population(P(k)).objs) - any(Offspring.objs>Population(P(k)).objs);
            if (associate_P(P(k)) == P(k)) && (associate_O == P(k))
               if CV_O < CV(P(k))
                   Population(P(k)) = Offspring;
                   count = count + 1;
               else
                   if CV_O == CV(P(k)) && (y == 1)
                      Population(P(k)) = Offspring;
                      count = count + 1;
                   elseif CV_O == CV(P(k)) && (y == 0)
                        newpop = [Population,Offspring];
                        CrowdDis1 = CrowdingDistance(newpop.objs);
                        if CrowdDis1(P(k)) < CrowdDis1(end)
                            Population(P(k)) = Offspring;
                            count = count + 1;
                        end
                   end
               end    
            elseif associate_O == P(k)
                Population(P(k)) = Offspring;
                count = count + 1;
            elseif associate_P(P(k)) == P(k) 
                v=1;
            else
                theta_P = rad2deg(acos(1-pdist2((Population(P(k)).objs-Z),W(P(k),:),'cosine')));
                theta_O = rad2deg(acos(1-pdist2((Offspring.objs-Z),W(P(k),:),'cosine')));
                if (theta_O < theta_P)&& (CV_O < CV(P(k)))
                    Population(P(k)) = Offspring;
                    count = count + 1;
                else
                    if (theta_O < theta_P)&& (CV_O == CV(P(k)))&& (y==1)
                         Population(P(k)) = Offspring;
                         count = count + 1;
                    elseif (theta_O < theta_P)&& (CV_O == CV(P(k)))&& (y==0)
                        newpop = [Population,Offspring];
                        CrowdDis1 = CrowdingDistance(newpop.objs);
                        if CrowdDis1(P(k)) < CrowdDis1(end)
                            Population(P(k)) = Offspring;
                            count = count + 1;
                        end
                    end
                end
            end
            if count == nr
               break 
            end
        end
            
    end
   
   
end

 function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end  

