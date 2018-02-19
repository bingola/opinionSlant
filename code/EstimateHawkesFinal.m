function graph=EstimateHawkesFinal(graphType,avg)
tic
avgm=avg;
if avg<1
    avgm=ceil(1/avg)+1;
end

FILE=strcat('../data/',graphType,'512Graphout',num2str(avgm));
load(FILE);
graph=graphout;
USS=graph.USSx{1};
A=graph.Adj;


for user=1:max(size(A))
%     user
if ~isempty(USS)    
 Nuser=find(A(:,user));
 n=length(Nuser);
 ahat=inference(user,graph);
 Ahat(Nuser,user)=ahat(2:end);
 muhat(user)=ahat(1);
end  
end

graph.Ahat=Ahat;
graph.muhat=muhat;
FILE=strcat('../data/',graphType,'512Graphout',num2str(avgm),'X');
graphout=graph;
eval(['save -v7.3 ',FILE,' graphout;']);
toc
end