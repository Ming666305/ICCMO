classdef DEST < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <constrained/large>
% Coevolutionary constrained multi-objective optimization framework加入向量删除和自适应
% type --- 1 --- Type of operator (1. GA 2. DE) 
            %% Initialization
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            type = Algorithm.ParameterSet(1);
            [W,~] = UniformPoint(2*Problem.N,Problem.M);

            %% Generate random population
            Population = Problem.Initialization();
            Fitness    = CalFitness(Population.objs);
            gen = 0;
            Z          = min(Population.objs,[],1);
            h = 1 ;
            last_gen         = 0.4 * ceil(Problem.maxFE/Problem.N);
            stage = 1;
            archive = Archive(Population,2*Problem.N);
            archive1 = Archive(Population,Problem.N,0);
            stage3_lim = 0.15*ceil(Problem.maxFE/(Problem.N));
            stage3_flag = 0;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if stage == 1
                  gen = gen + 1 ;
                  if gen >=  last_gen 
                        stage = 2;
                  end
                  if type == 1
                      MatingPool = TournamentSelection(2,2*Problem.N,Fitness);             
                      Offspring  = OperatorGAhalf(Problem,Population(MatingPool));
                      archive = Archive([Offspring,archive],2*Problem.N,10000);
                      archive1 = Archive([Offspring,archive1],Problem.N,0);
                      Z = min([Z;Offspring.objs],[],1);
                  elseif type == 2
                      MatingPool1 = TournamentSelection(2,Problem.N,Fitness); 
                      MatingPool2 = TournamentSelection(2,Problem.N,Fitness); 
                      MatingPool3 = TournamentSelection(2,Problem.N,Fitness); 
                      Offspring  = OperatorDE(Problem,Population(MatingPool1),Population(MatingPool2),Population(MatingPool3));
                     
%                         Offspring  = OperatorDE(Problem,Population,Population(randperm(end)),Population(randperm(end)),Population(randperm(end)));
                      archive = Archive([Offspring,archive],2*Problem.N,10000);
                      archive1 = Archive([Offspring,archive1],Problem.N,0);
                      Z = min([Z;Offspring.objs],[],1);
                 
                  

                  end
                  [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N,false);  

                elseif stage == 2
                    Z = min([Z;Population.objs],[],1);
                    if h == 1
                      Population = archive;
                      Angle = acos(1-pdist2(W,(Population.objs-repmat(Z,2*Problem.N,1)),'cosine'));
                      [~,associate] = min(Angle,[],2);
                      Population = Population(associate);
                      h = 2;
                    end
                    Population =  Reproduce_stage1(Problem,Population,W,type,Z);
                    archive1 = Archive([Population,archive1],Problem.N,0);
                    if theta_maxchange(Population,W,Z) || Problem.FE > 0.75*Problem.maxFE
                        stage = 3;
                    end
                else
                    Z = min([Z;Population.objs],[],1);
                    Population =  Reproduce_stage2(Problem,Population,W,type,Z);
                    archive1 = Archive([Population,archive1],Problem.N,0);
                    stage3_flag = stage3_flag+1;
                    if stage3_flag == stage3_lim
                         [Population,W] = ReferenceVectorAdaptation(Problem,Population,W);
                    end
                end
                if Problem.FE >= Problem.maxFE
                    Population = archive1;
                end
            end
        end
    end
end

function flage = theta_maxchange(Population,W,Z)
   PopObj = Population.objs;    
   [N,~]  = size(PopObj);
    Angle = acos(1-pdist2((Population.objs-repmat(Z,N,1)),W,'cosine'));
    [~,associate] = min(Angle,[],2);
    flage = false;
    a = 0;
    for i = 1 : N
        if associate(i) == i
           a = a+1;
        end
    end
    if(a/N)>0.8
       flage = true;  
    end
end

