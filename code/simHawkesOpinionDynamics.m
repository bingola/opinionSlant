function graphData=simHawkesOpinionDynamics(graphType,avg)
avgm=avg;
 
 
FILE=strcat('../data/',graphType,'512Graphout',num2str(avgm),'X');
 load(FILE);% 
graph=graphout;
 N=length(graph.mu);
a=graph.A;
alpha=graph.alpha;
mu=graph.mu;
maxUses=avg*N;
USS=graph.USSx{1};
eventVec=[];
userVec=[];
for n=1:N
eventVec=[eventVec;USS{n}];
userVec=[userVec;n*ones(length(USS{n}),1)];
end

[x1,y1]=sort([eventVec]);
tu_events(:,2)=x1;
tu_events(:,1)=userVec(y1);
% Initialization
OSS=cell(N,1);
for n=1:N
    OSS{n}=zeros(maxUses,1);    
end
tu=zeros(1,N);%number of events per user
opinionTransient=zeros(1,N);
opinionIni=alpha;
TotalEvents=length(x1);
t=0;
om=4;
%% general subroutine
 while t<TotalEvents
    t=t+1;
    t
    user=tu_events(t,1); 
    if t==1
        s=tu_events(t,2);
    end
    if t>1
    s=tu_events(t,2)-tu_events(t-1,2);
    end    
    opinionTransient=opinionTransient*exp(-om*s);
    currentOpinion=opinionIni+opinionTransient;
    %%%%%change here%%%%%
    sentimentValue=normrnd(currentOpinion(user),1); 
    %%%%%%%%%%%%%%%%%%%%
    opinionTransient=opinionTransient+a(user,:)*sentimentValue;
    
    sentiment_m(1,t)=sentimentValue;
    sentiment_m(2,t)=user;   
    tu(user)=tu(user)+1;
    OSS{user}(tu(user),1)=sentimentValue;
 
end


% events per user
for n=1:N
    OSS{n}=OSS{n}(1:tu(n),:);
end
%%
%
clear g
clear G
g{N}=[]; 
G{N}=[];
clc
A=graph.Adj;
for user=1:N
    Nuser=find(A(:,user)); 
   for j=1:length(Nuser) 
      Nuser(j)
    g1=(ones(length(USS{Nuser(j)}),1)*USS{user}'-USS{Nuser(j)}*ones(1,length(USS{user})));
    g1(g1<0)=Inf;
    SentimentMatrix=(OSS{Nuser(j)}*ones(1,length(USS{user})));
    g{user}(j,:)=sum(exp(-om*g1).*SentimentMatrix,1);

   end
end
 
graphData=graph;
graphData.gOpinion=g;
graphData.alpha=alpha;
graphData.OSS=OSS;
graphData.sentiment=sentiment_m;

USS=graph.USSx{1};
for user=1:N
  if ~isempty({user})  
   Nuser=find(A(:,user));
   T=max(USS{user});
   for j=1:length(Nuser)
      if ~isempty(USS{Nuser(j)})
       T=max(T,max(USS{Nuser(j)}));
     end
   end
  end
end
% a=AIn;
dt=(1e-3)*T;
ts=0:dt:T;
% mm=[user length(ts)]
xstar(user,length(ts))=0;

 for user=1:N
%      user
    Nuser=find(graph.Adj(:,user));
    clear G
    for j=1:length(Nuser)
        USSX=graph.USSx{1};
        tv=USSX{Nuser(j)};

        gMat=(ones(length(USS{Nuser(j)}),1)*ts-USS{Nuser(j)}*ones(1,length(ts)));
        gMat(gMat<0)=Inf;
        SM=ones(length(ts),1)*OSS{Nuser(j)}';
%         size(SM)
%         size(gMat)
        g=exp(-om*gMat).*SM';
        v=g';
       G(j,:)=sum(v,2);
         
    end
    
    if isempty(Nuser)
    xstar(user,:)=alpha(user)*ones(1,length(ts));
    else
     xstar(user,:)=alpha(user)+a(Nuser,user)'*G;  
    end
 end

 graphData.xstar=xstar;
 graphData.ts=ts;
%%
 
FILE=strcat('../data/',graphType,'512Graphout',num2str(avgm),'XContainedOpinion');

eval(['save -v7.3 ',FILE,' graphData;']);
end