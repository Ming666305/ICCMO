classdef ICCMO < ALGORITHM
    % <2026> <multi> <real> <constrained>
    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population{1} = Problem.Initialization(); % Main task
            Fitness{1}    = CalFitness(Population{1}.objs,Population{1}.cons);
            Population{2} = Problem.Initialization(); % Global auxiliary task
            Fitness{2}   = CalFitness(Population{2}.objs);
            Population{3} = Problem.Initialization(); % Local auxiliary task
            Fitness{3}   = CalFitness(Population{3}.objs,Population{3}.cons);
            %% Evaluate the Population
            Nt = 2;      % The number of tasks
            Beta = 0.5;   % Learning phase
            alpha = 30000;
            Pmin = 0.1;
            % Pmins = ((Problem.maxFE/2-Problem.FE)/Problem.maxFE)^2;% Minimum selection probability
            Pmins = 0.1;
            A = Population{1};
            rwd   =  0.5.*ones(2, 1);
            % counter = 0;  % 用于记录连续满足条件的代数
            % nn = 5;       % 设定连续代数阈值
            rwds  =  0.5*ones(4, 1);
            cons = Population{1}.cons;
            cons(cons<0) = 0;
            VAR0 = max(sum(cons,2));
            if VAR0 == 0
                VAR0 = 1;
            end
            X=0;

            %% Optimization
            while Algorithm.NotTerminated(A)
                Q = [];
                cp  = (-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*(1-X)^cp;
                pro   =   Pmin / Nt + (1 - Pmin) * rwd ./ sum(rwd);
                pro   =   pro ./ sum(pro);

                pros   =   Pmins / 4 + (1 - Pmins) * rwds ./ sum(rwds);
                pros   =   pros ./ sum(pros);
                r = rand;
                for i = 1:4
                    if r <= sum(pros(1:i))
                        j = i;
                        break;
                    end
                end
                MatingPool1 = TournamentSelection(2,Problem.N,Fitness{1});
                MatingPool2 = TournamentSelection(2,Problem.N,Fitness{2});
                MatingPool3 = TournamentSelection(2,Problem.N,Fitness{3});
                if  alpha > Problem.FE
                    Offspring1(1:Problem.N/2) = GenerateOffspring(Problem,Population{1},MatingPool1,'GA');
                    Offspring1(Problem.N/2+1:Problem.N) = GenerateOffspring(Problem,Population{1},MatingPool1,'DE');
                    Offspring2(1:Problem.N/2) = GenerateOffspring(Problem,Population{2},MatingPool2,'GA');
                    Offspring2(Problem.N/2+1:Problem.N) = GenerateOffspring(Problem,Population{2},MatingPool2,'DE');
                    Offspring3(1:Problem.N/2) = GenerateOffspring(Problem,Population{3},MatingPool3,'GA');
                    Offspring3(Problem.N/2+1:Problem.N) = GenerateOffspring(Problem,Population{3},MatingPool3,'DE');
                    Offspring = [Offspring1,Offspring2,Offspring3];
                    [Population{1},Fitness{1}] = EnvironmentalSelection([Population{1},Offspring],Problem.N,true);
                    [Population{2},Fitness{2}] = EnvironmentalSelection([Population{2},Offspring],Problem.N,false);
                    [Population{3},Fitness{3}] = EnvironmentalSelection([Population{3},Offspring],Problem.N,true);
                    [A,~] =  ArchiveUpdate([A Offspring],Problem.N);
                    % MatingPool3 = TournamentSelection(2,min(length(Population{3}),Problem.N/2),Fitness{3});
                elseif alpha <= Problem.FE && Problem.FE <= Beta * Problem.maxFE
                    if j == 1
                        Offspring1 = GenerateOffspring(Problem,Population{1},MatingPool1,'GA');
                        Offspring3 = GenerateOffspring(Problem,Population{3},MatingPool3,'GA');
                        % EliteFrom3 = GetElite(Offspring3, Problem.N/4);
                        % DiverseFrom2 = diversity(Problem, Offspring3, Problem.N/2);
                        Offspring = [Offspring1, Offspring3];
                    elseif j == 2
                        MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness{1});
                        MatingPool3 = TournamentSelection(2,2*Problem.N,Fitness{3});
                        Offspring1 = GenerateOffspring(Problem,Population{1},MatingPool1,'DE');
                        %Offspring2 = GenerateOffspring(Problem,Population{2},MatingPool2,'DE');
                        Offspring3 = GenerateOffspring(Problem,Population{3},MatingPool3,'DE');
                        % EliteFrom3 = GetElite(Offspring3, Problem.N/4);
                        %DiverseFrom2 = diversity_generation(Problem, Offspring2, Problem.N/2);
                        Offspring = [Offspring1, Offspring3];
                    elseif j==3
                        Offspring2 = GenerateOffspring(Problem,Population{2},MatingPool2,'GA');
                        Offspring1 = GenerateOffspring(Problem,Population{3},MatingPool3,'GA');
                        %EliteFrom3 = GetElite(Offspring3, Problem.N/4);
                        % DiverseFrom2 = diversity_generation(Problem, Offspring2, Problem.N/2);
                        Offspring = [Offspring1, Offspring2];
                    elseif j==4
                        MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness{2});
                        Offspring1  = OperatorDE(Problem,Population{2},Population{2}(MatingPool1(1:end/2)),Population{2}(MatingPool1(end/2+1:end)));
                        %Offspring2 = GenerateOffspring(Problem,Population{2},MatingPool2,'DE');
                        MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness{3});
                        Offspring2  = OperatorDE(Problem,Population{3},Population{3}(MatingPool1(1:end/2)),Population{3}(MatingPool1(end/2+1:end)));
                        %EliteFrom3 = GetElite(Offspring2, Problem.N/4);
                        % DiverseFrom2 = diversity_generation(Problem, Offspring2, Problem.N/2);
                        Offspring = [Offspring1, Offspring2];

                    end
                    DiverseFrom = diversity(Problem, Offspring, Problem.N);
                    Q = [Q DiverseFrom];
                    if size(Q,2) > 0
                        s = size(A,2);
                        [A,Next] =  ArchiveUpdate([A,Q],Problem.N);
                        if s >= Problem.N
                            rwds(j) = rwds(j) + sum(Next(Problem.N+1:end))/Problem.N;
                        end
                    end
                    disp(['Current rwds: ', num2str(rwds')]);

                    [Population{1},Fitness{1}] = EnvironmentalSelection([Population{1},Offspring1],Problem.N,true);
                    [Population{2},Fitness{2}] = EnvironmentalSelection([Population{2},Offspring],Problem.N,false);
                    [Population{3},Fitness{3}] = EnvironmentalSelection([Population{3},Offspring],Problem.N,true);
                else
                    r = rand;
                    for t = 1:Nt
                        if r <= sum(pro(1:t))
                            k = t;
                            break;
                        end
                    end
                    a = pros(1)+pros(3);
                    b = pros(2)+pros(4);
                    if a > b
                        Offspring1 = GenerateOffspring(Problem,Population{1},MatingPool1,'GA');
                        if k == 1
                            Offspring2 = GenerateOffspring(Problem,Population{2},MatingPool2,'GA');
                            Offspring=[Offspring1,Offspring2];
                        end
                        if k == 2
                            % MatingPool3 = TournamentSelection(2,min(length(Population{3}),Problem.N/2),Fitness{3});
                            Offspring3 = GenerateOffspring(Problem,Population{3},MatingPool3,'GA');
                            Offspring=[Offspring1,Offspring3];
                        end
                    else
                        Offspring1 = GenerateOffspring(Problem,Population{1},MatingPool1,'DE');
                        if k == 1
                            Offspring2 = GenerateOffspring(Problem,Population{2},MatingPool2,'DE');
                            Offspring=[Offspring1,Offspring2];
                        end
                        if k == 2
                            % MatingPool3 = TournamentSelection(2,min(length(Population{3}),Problem.N/2),Fitness{3});
                            Offspring3 = GenerateOffspring(Problem,Population{3},MatingPool3,'DE');
                            Offspring=[Offspring1,Offspring3];
                        end
                    end

                    [Population{1},Fitness{1}] = EnvironmentalSelection( [Population{1},Offspring],Problem.N,true);
                    [Population{2},Fitness{2}] = EnvironmentalSelection( [Population{2},Offspring],Problem.N,false);
                    [Population{3},Fitness{3}] = Auxiliray_task_EnvironmentalSelection([Population{3},Offspring],Problem.N,VAR);

                    num = size(Offspring1,2);
                    Q = [Q Offspring];
                    s = size(A,2);
                    [A,Next] =  ArchiveUpdate([A Q],Problem.N);

                    if size(Q,2) > 0 && Problem.FE > 0.1 * Problem.maxFE
                        if s >= Problem.N
                            rwd(k)   =  rwd(k) + sum(Next(Problem.N+num+1:end))/Problem.N;   % update the reward of the main task
                        end
                    end
                end
                X=Problem.FE/Problem.maxFE;
            end
        end
    end
end

function Offspring = GenerateOffspring(Problem, Population, MatingPool, operator)
if strcmp(operator, 'GA')
    Offspring = OperatorGAhalf(Problem, Population(MatingPool));
else
    Offspring = OperatorDE(Problem, Population(1:Problem.N/2), Population(MatingPool(1:end/2)), Population(MatingPool(end/2+1:end)));
end
end


function Offspring = diversity(Problem, Population, N)
% 计算理想点
Z = [0,0];
M = size(Population(1).objs, 2); % 目标数
W = UniformPoint(N, M);
Selected = environmental_selection(Population, N, W, Z);
Fitness = CalFitness(Selected.objs,Selected.cons);
MatingPool = TournamentSelection(2, N, Fitness);
Offspring =  GenerateOffspring(Problem, Population, MatingPool, 'DE');
end

function Selected = environmental_selection(Population, N, W, Z)
PopObj = Population.objs;
Angle = acos(1 - pdist2(PopObj, W, 'cosine')); % 计算角度
[~, associate] = min(Angle, [], 2); % 关联到最近的参考向量

% 选择每个参考向量中最优的个体
Selected = [];
uniqueRefs = unique(associate);
for i = 1:length(uniqueRefs)
    currentRef = uniqueRefs(i);
    candidates = find(associate == currentRef);
    if ~isempty(candidates)
        % 选择距离最近的个体
        [~, idx] = min(sqrt(sum((PopObj(candidates, :) - W(currentRef, :)).^2, 2)));
        Selected = [Selected, Population(candidates(idx))];
    end
    if length(Selected) >= N
        break;
    end
end
% 如果不足N个，随机补充
if length(Selected) < N
    remaining = N - length(Selected);
    extra = randi(length(Population), 1, remaining);
    Selected = [Selected, Population(extra)];
end
Selected = Selected(1:N); % 确保精确返回N个个体
end