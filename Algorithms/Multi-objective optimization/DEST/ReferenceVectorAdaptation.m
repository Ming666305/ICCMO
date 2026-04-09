function [Population,W] = ReferenceVectorAdaptation(Problem,Population,W)
% Reference vector adaption strategy

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%     V = V.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V,1),1);
    N = size(W,1);
    lim = 0.4*size(W,1);
    [FrontNo,~] = NDSort(Population.objs,Population.cons,inf);
    flag = false(N,1);
    ra = zeros(N,1);
    for i = 1 : N
        if FrontNo(i) == 1 
            ra(i) = 1;
              if i == 1
                 if  ra(i+1) == 0
                     ra(i+1) = 2;
                 end
                 if  ra(i+2) == 0
                     ra(i+2) = 2;
                 end
                 if  ra(i+3) == 0
                     ra(i+3) = 2;
                 end
             elseif i ==2
                 if  ra(i+1) == 0
                     ra(i+1) = 2;
                 end
                 if ra(i+2) == 0
                     ra(i+2) = 2;
                 end
                 if  ra(i+3) == 0
                     ra(i+3) = 2;
                 end
                 if ra(i-1) == 0
                     ra(i-1) = 2;
                 end
              elseif i ==3
                  if  ra(i+1) == 0
                     ra(i+1) = 2;
                 end
                 if ra(i+2) == 0
                     ra(i+2) = 2;
                 end
                 if  ra(i+3) == 0
                     ra(i+3) = 2;
                 end
                 if ra(i-1) == 0
                     ra(i-1) = 2;
                 end
                 if ra(i-2) == 0
                     ra(i-2) = 2;
                 end
              elseif i == N-1
                  if  ra(i+1) == 0
                     ra(i+1) = 2;
                  end
                  if ra(i-3) == 0
                     ra(i-3) = 2;
                 end
                 if ra(i-2) == 0
                     ra(i-2) = 2;
                 end
                 if ra(i-1) == 0
                     ra(i-1) = 2;
                  end
              elseif i== N
                  if  ra(i-1) == 0
                     ra(i-1) = 2;
                  end
                 if ra(i-2) == 0
                     ra(i-2) = 2;
                 end
                 if ra(i-3) == 0
                     ra(i-3) = 2;
                 end
              elseif i== N-2
                  if  ra(i-1) == 0
                     ra(i-1) = 2;
                  end
                 if ra(i-2) == 0
                     ra(i-2) = 2;
                 end
                 if ra(i-3) == 0
                     ra(i-3) = 2;
                 end
                 if  ra(i+1) == 0
                     ra(i+1) = 2;
                 end
                 if  ra(i+2) == 0
                     ra(i+2) = 2;
                 end
             else
                  if  ra(i+1) == 0
                     ra(i+1) = 2;
                  end
                 if ra(i-2) == 0
                     ra(i-2) = 2;
                 end
                 if ra(i-1) == 0
                     ra(i-1) = 2;
                 end
                  if ra(i+2) == 0
                     ra(i+2) = 2;
                  end
                  if ra(i-3) == 0
                     ra(i-3) = 2;
                 end
                  if ra(i+3) == 0
                     ra(i+3) = 2;
                  end
              end
        end
    end
    for i = 1 : N
        if (ra(i) == 1) || (ra(i) == 2)
            flag(i) = true;
        end
    end
    if sum(flag==true) >= lim
        W=W(flag,:);
        Population = Population(flag);
    end
     W = W.*repmat(max(Population.objs,[],1)-min(Population.objs,[],1),size(W,1),1);
end  