function flag = Analysis_of_relationship(temp1,temp2,c1,c2,target1,target2,beta)

pop1 = temp1.decs;
conv1 = temp1.cons;
conv1(conv1<=0)=0;
obj1 = temp1.objs;

pop2 = temp2.decs;
conv2 = temp2.cons;
conv2(conv2<=0)=0;
obj2 = temp2.objs;

[FrontNo,MaxFNo] = NDSort(obj1,conv1(:,target1),inf);
x1 = find(FrontNo==1);
[FrontNo,MaxFNo] = NDSort(obj2,conv2(:,target2),inf);
x2 = find(FrontNo==1);

obj1 = obj1(x1,:);
obj2 = obj2(x2,:);
[FrontNo,MaxFNo] = NDSort([obj1;obj2],inf);


alpha1 = length(find(FrontNo(1:length(x1))==1))/length(x1);
alpha2 = length(find(FrontNo(length(x1)+1:end)==1))/length(x2);

if alpha1 >= beta && alpha2 < beta
    flag = 0;
elseif  alpha1 < beta && alpha2 >= beta
    flag = 1;  %3
elseif  alpha1 < beta && alpha2 < beta
    flag = 2; %0
elseif  alpha1 >= beta && alpha2 >= beta
    if c1 == 1
        flag = 1;
    elseif c2 == 1
        flag = 0;
    else
        if alpha1 >= alpha2
            flag = 1;
        else
            flag = 0;
        end
    end
end

end


