clc
clear all
 graphTypeSet={'KroneckerCP'};
Avg=[25];
for ix=1
for avg=Avg
    graphType=graphTypeSet{ix};
        FILE=strcat('../data/',graphType,'512Graph');
        load(FILE)
        graphx=graph;    
        simulateHawkesEventsOnly(graphx,graphType,avg);
        EstimateHawkesFinal(graphType,avg);
        simHawkesOpinionDynamics(graphType,avg);
        estimateOpinionParams(graphType,avg);
end
end