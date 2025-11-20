clear
clc
close all

con=configure();
repeatMax=con.repeat;
functions=con.TestFunctions;
CreatTruePOF()
Alfa=0.7;
mum = 20;
for group=1:size(con.T_parameter,1)
    diary on;  T_parameter=con.T_parameter(group,:)
    diary off;
    for testfunc=1:size(functions,2)
        Problem=TestFunctions(functions{testfunc});
        diary on;   fprintf('\n');   fprintf(Problem.Name);  diary off;
        mkdir( ['.\IGDResults\' '\parameter' num2str(T_parameter) '\' Problem.Name '\']);
        %         mkdir( ['.\PFResults\' '\T_parameter' num2str( T_parameter) '\' Problem.Name '\']);
        
        for rep=1:repeatMax   %repeatMax Maximum number of repeated runs
            t = 0;         % the initial moment
            Problem=TestFunctions(functions{testfunc});
            CostFunction=Problem.FObj;    % Cost Function
            gen = T_parameter(1,2);     % max generations of evolution
            M = Problem.NObj;         % number of objects
            V = size(Problem.XLow,1);          % number of decision variables
            pop = 100;      % poplation number
            cp=1;
            if M==3
                pop=150;
            end
            min_range = Problem.XLow;
            max_range = Problem.XUpp;
            TW=20;     % The number of weight vectors in the neighborhood
            [Popsize,weights,neighbour]=Init_weights(pop,M,TW);    % Initialize weight vector matrix
            chromosome = initialize_variables(Popsize, M, V, min_range, max_range,Problem,t); % Initialize the population randomly
            [z,~]=min(chromosome(:,V+1:end));
            for T = 1:T_parameter(1,3)/T_parameter(1,2)
                t = 1/T_parameter(1,1)*(T);
                if T==1
                    %The initial two generations of the algorithm do not undergo mutation because the second generation cannot calculate the centroids of the first two generations.
                    %The 'if' here is only for programming convenience, and the algorithm does not know when the environment changes during runtime
                    for i = 1:2
                        chrom = chromosome(:,1:V); obj = chromosome(:,V+1:M + V);
                        for j=1:Popsize 
                            newpoint=evolution(chrom,neighbour,j,V,TW);
                            newpoint=fixnew(newpoint,min_range',max_range');
                            newpoint=gaussian_mutate(newpoint,V,min_range',max_range');
                            CostFunction=Problem.FObj;
                            zz = CostFunction(newpoint,t);
                            zz = zz';
                            z=min(z,zz);
                            nei=neighbour(j,:);
                            for t1=1:TW
                                k=nei(t1);
                                [chrom(k,:),obj(k,:)]=updates(weights(k,:),chrom(k,:),obj(k,:),newpoint,zz,z);
                            end
                        end
                        chromosome = [chrom,obj];                       
                        centerPoints(cp,:)=mean(chromosome(:,1:V),1);
                        cp=cp+1;
                    end
                    for i = 3:gen
                        previous_index = 0;
                        chromosome = non_domination_sort_mod(chromosome, M, V);
                        current_index = max(find(chromosome(:,M + V + 1) == 1));
                        distance=abs(centerPoints(cp-1,:)-centerPoints(cp-2,:));
                        Y=sign(distance);
                        p=(size(chromosome,1)-current_index)/size(chromosome,1);
                        A = Alfa*(1-p)*exp(p);
                        if current_index<=pop-round(A*pop)
                            index=round(A*pop);
                        else
                            index=pop-current_index;
                            temp_pop = ...
                                chromosome(previous_index + 1 : current_index, :);
                            [temp_sort,temp_sort_index] = ...
                                sort(temp_pop(:, M + V + 2),'ascend');
                            randomI=norm(distance,2).*rand(1,V);
                            for j = 1 : current_index-pop+round(A*pop)
                                chromosome(temp_sort_index(j),1: V)=chromosome(temp_sort_index(j),1: V)+randomI;
                                for k=1:V
                                    if chromosome(temp_sort_index(j),k)>max_range(k)
                                        chromosome(temp_sort_index(j),k)=max_range(k);
                                    else if chromosome(temp_sort_index(j),k)<min_range(k)
                                            chromosome(temp_sort_index(j),k)=min_range(k);
                                        end
                                    end
                                end
                            end
                        end
                        for n = pop+1-index:pop
                            child_3 = chromosome(n,:);
                            for j = 1 : V
                                r(j) = rand(1);
                                if r(j) < 0.5
                                    delta(j) = (2*r(j))^(1/(mum+1)) - 1;
                                else
                                    delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
                                end
                                chromosome(n,j) = child_3(j) + delta(j);
                                if chromosome(n,j) > max_range(j)
                                    chromosome(n,j) = max_range(j);
                                elseif chromosome(n,j) < min_range(j)
                                    chromosome(n,j) = min_range(j);
                                end
                            end
                        end
                        
                        chrom = chromosome(:,1:V); obj = chromosome(:,V+1:M + V);
                        for j=1:Popsize
                            newpoint=evolution(chrom,neighbour,j,V,TW);
                            newpoint=fixnew(newpoint,min_range',max_range');
                            newpoint=gaussian_mutate(newpoint,V,min_range',max_range');
                            CostFunction=Problem.FObj;
                            zz = CostFunction(newpoint,t);
                            zz = zz';
                            z=min(z,zz);
                            nei=neighbour(j,:);
                            for t1=1:TW
                                k=nei(t1);
                                [chrom(k,:),obj(k,:)]=updates(weights(k,:),chrom(k,:),obj(k,:),newpoint,zz,z);
                            end
                        end
                        chromosome = [chrom,obj];
                        centerPoints(cp,:)=mean(chromosome(:,1:V),1);
                        cp=cp+1;
                    end
                else
                    for i = 1 : gen
                        previous_index = 0;
                        chromosome = non_domination_sort_mod(chromosome, M, V);
                        current_index = max(find(chromosome(:,M + V + 1) == 1));
                        distance=abs(centerPoints(cp-1,:)-centerPoints(cp-2,:));
                        Y=sign(distance);
                        p=(size(chromosome,1)-current_index)/size(chromosome,1);
                        A = Alfa*(1-p)*exp(p);
                        if current_index<=pop-round(A*pop)
                            index=round(A*pop);
                        else
                            index=pop-current_index;
                            temp_pop = ...
                                chromosome(previous_index + 1 : current_index, :);
                            [temp_sort,temp_sort_index] = ...
                                sort(temp_pop(:, M + V + 2),'ascend');
                            randomI=norm(distance,2).*rand(1,V);
                            for j = 1 : current_index-pop+round(A*pop)
                                chromosome(temp_sort_index(j),1: V)=chromosome(temp_sort_index(j),1: V)+randomI;
                                for k=1:V
                                    if chromosome(temp_sort_index(j),k)>max_range(k)
                                        chromosome(temp_sort_index(j),k)=max_range(k);
                                    else if chromosome(temp_sort_index(j),k)<min_range(k)
                                            chromosome(temp_sort_index(j),k)=min_range(k);
                                        end
                                    end
                                end
                            end
                        end
                        for n = pop+1-index:pop
                            child_3 = chromosome(n,:);
                            for j = 1 : V
                                r(j) = rand(1);
                                if r(j) < 0.5
                                    delta(j) = (2*r(j))^(1/(mum+1)) - 1;
                                else
                                    delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
                                end
                                chromosome(n,j) = child_3(j) + delta(j);
                                if chromosome(n,j) > max_range(j)
                                    chromosome(n,j) = max_range(j);
                                elseif chromosome(n,j) < min_range(j)
                                    chromosome(n,j) = min_range(j);
                                end
                            end
                        end
                                               chrom = chromosome(:,1:V); obj = chromosome(:,V+1:M + V);
                        for j=1:Popsize
                            newpoint=evolution(chrom,neighbour,j,V,TW);
                            newpoint=fixnew(newpoint,min_range',max_range');
                            newpoint=gaussian_mutate(newpoint,V,min_range',max_range');
                            CostFunction=Problem.FObj;
                            zz = CostFunction(newpoint,t);
                            zz = zz';
                            z=min(z,zz);
                            nei=neighbour(j,:);
                            for t1=1:TW
                                k=nei(t1);
                                [chrom(k,:),obj(k,:)]=updates(weights(k,:),chrom(k,:),obj(k,:),newpoint,zz,z);
                            end
                        end
                        chromosome = [chrom,obj];
                        centerPoints(cp,:)=mean(chromosome(:,1:V),1);
                        cp=cp+1;
                    end
                end
                chromosome = non_domination_sort_mod(chromosome, M, V);
                current_index = max(find(chromosome(:,M + V + 1) == 1));
                POF_iter = chromosome(1:current_index,V+1:V+M);
                TruePOF = getBenchmarkPOF(Problem.Name,1,T,T_parameter);
                iterIGD(T) = IGD(POF_iter,TruePOF);
            end
            resIGD(rep)=mean(iterIGD);
            diary on;   fprintf('\n %.3d',resIGD(rep));  diary off;
            %% Save Metrics to flies
            FILEPATH = ['.\IGDResults\' '\parameter' num2str(T_parameter) '\' Problem.Name '\'];
            filename = ['IGD' num2str(rep) '.mat'];
            save([FILEPATH,filename],'iterIGD','-mat');
        end
        rIGD=mean(resIGD); sIGD=std(resIGD);
        diary on;   fprintf('\n %.3d',rIGD);  fprintf('\t %.3d',sIGD);    diary off;
    end
end 
