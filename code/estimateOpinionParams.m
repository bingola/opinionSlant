function graph=estimateOpinionParams(graphType,avg)
tic
avgm=avg;
if avg<1
    avgm=ceil(1/avg)+1;
end
 
FILE=strcat('../data/',graphType,'512Graphout',num2str(avgm),'XContainedOpinion');
load(FILE);
graph=graphData;
USS=graph.USSx{1};
obj1=0;
obj2=0;
A=graph.Adj;
N=max(size(A));
regVal=1;
clear Ahat;
clear muhat;
Ahat(N,N)=0;
muhat(N)=0;
for user=1:max(size(A))
    user
if ~isempty(USS{user})    
 Nuser=find(A(:,user));
 n=length(Nuser);
 ahat=inferOpParams(user,graph,regVal);
 Ahat(Nuser,user)=ahat(2:end);
 muhat(user)=ahat(1);
end  
end

graph.AhatOpinion=Ahat;
graph.muhatOpinion=muhat;

FILE=strcat('../data/',graphType,'512Graphout',num2str(avgm),'XContainedOpinionX');

eval(['save -v7.3 ',FILE,' graph;']);
toc
end
%%
%evaluation

