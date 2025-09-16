%% Phase 1

close all;
clear;
clc;

functions = {'F-1', 'F-2', 'F-3', 'F-4', 'F-5', 'F-6', 'F-7', 'F-8', 'F-9', 'F-10', 'F-11', 'F-12', 'F-13', 'F-14', 'F-15'};
numFunctions = length(functions);

downstreamFunctions = {
{'F-8'},                                                             % F-1
{'F-3', 'F-4', 'F-5', 'F-6', 'F-7', 'F-8', 'F-10', 'F-11', 'F-15'},  % F-2
{'F-5', 'F-6', 'F-7', 'F-10', 'F-11', 'F-15'},                       % F-3
{'F-5', 'F-6', 'F-7', 'F-10', 'F-15'},                               % F-4
{'F-2', 'F-8'},                                                      % F-5
{'F-8'},                                                             % F-6
{'F-8'},                                                             % F-7
{'F-9', 'F-10', 'F-11', 'F-15'},                                     % F-8
{'F-1', 'F-3', 'F-6', 'F-8'},                                        % F-9
{'F-9'},                                                             % F-10
{'F-12'},                                                            % F-11
{},                                                                  % F-12
{'F-1'},                                                             % F-13
{'F-2', 'F-9'},                                                      % F-14
{'F-1', 'F-3', 'F-4', 'F-6','F-7', 'F-8', 'F-9', 'F-10'}             % F-15
};


edges = [];
for f = 1:numFunctions
for d = 1:length(downstreamFunctions{f})
downstreamIdx = find(strcmp(functions, downstreamFunctions{f}{d}));
if ~isempty(downstreamIdx)
edges = [edges; f, downstreamIdx];
end
end
end

G = digraph(edges(:,1), edges(:,2), [], functions);

nodeColors = lines(numFunctions); 
edgeColors = [0.3 0.3 0.3]; 

figure;
h1 = plot(G, 'Layout', 'layered', 'LineWidth', 1.5, 'ArrowSize', 12);
title('FRAM Function Relationship Graph', 'FontSize', 14);
h1.MarkerSize = 10;
h1.NodeColor = nodeColors;
h1.EdgeColor = [0.3 0.3 0.6];
EdgeAlpha = 0.7;  
h1.NodeFontSize = 11;
h1.NodeLabel = {};
labelOffset = 0.2; 
for i = 1:numFunctions
text(h1.XData(i) + labelOffset, h1.YData(i) + labelOffset, functions{i}, ...
'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
end

% Time (V_t): [Too Early, On-Time, Too Late, Not at all] -> Values: [1,2,3,4]
timeVals = [1, 2, 3, 4];
timeProbs = [0.08, 0.80, 0.09, 0.03];
if abs(sum(timeProbs) - 1) > 1e-6, error('Time probabilities do not sum to 1. Sum = %.6f', sum(timeProbs)); end

% Precision (V_p): [Precise, Acceptable, Imprecise, Wrong] -> Values: [1,2,3,4]
precVals = [1, 2, 3, 4];
precProbs = [0.15, 0.75, 0.08, 0.02];
if abs(sum(precProbs) - 1) > 1e-6, error('Precision probabilities do not sum to 1. Sum = %.6f', sum(precProbs)); end

% Force (V_f): [Too Low, On Target, Too High, Not at all] -> Values: [1,2,3,4]
forceVals = [1, 2, 3, 4];
forceProbs = [0.20, 0.65, 0.14, 0.01];
if abs(sum(forceProbs) - 1) > 1e-6, error('Force probabilities do not sum to 1. Sum = %.6f', sum(forceProbs)); end

% Sequence (V_s): [Correct Order, Wrong Order] -> Values: [1,2]
seqVals = [1, 2];
seqProbs = [0.90, 0.10];
if abs(sum(seqProbs) - 1) > 1e-6, error('Sequence probabilities do not sum to 1. Sum = %.6f', sum(seqProbs)); end

numSim = 1000; % number of Monte Carlo iterations
OV = zeros(numFunctions, numSim); 
for f = 1:numFunctions
    vt = randsample(timeVals, numSim, true, timeProbs);
    vp = randsample(precVals, numSim, true, precProbs);
    vf = randsample(forceVals, numSim, true, forceProbs);
    vs = randsample(seqVals, numSim, true, seqProbs);
    OV(f, :) = vt .* vp .* vf .* vs;
end

Final_OV = zeros(numFunctions, 1);
for f = 1:numFunctions
    threshold = prctile(OV(f,:), 90);
    worst_outcomes = OV(f, OV(f,:) >= threshold); % Use >= for CVaR
    if isempty(worst_outcomes)
        Final_OV(f) = threshold;
    else
        Final_OV(f) = mean(worst_outcomes);
    end
end

CI_expert_judgment = [1; 2; 1; 0.5; 2; 2; 1; 2; 0.5; 1; 1; 1; 0.5; 1; 2];
if length(CI_expert_judgment) ~= numFunctions, error('Expert CI judgments do not match the number of functions.'); end
CI = CI_expert_judgment;

b = [0,1,0; 0.5,0.5,1; 1,0,1; 1,0,1; 0.5,1,1; 1,0.5,0; 1,1,1; 1,1,1; 0.5,1,0; 1,0.5,0; 0,0,0; 0,0,0; 0,1,0; 0,0.5,1; 1,0.5,0];
m = size(b, 2);

S = [ 2, 2, 2; 4, 2, 2; 2, 4, 2; 2, 2, 4; ];
numScenarios = size(S, 1);

e = zeros(numFunctions, numScenarios);
for s = 1:numScenarios
    for f = 1:numFunctions
        e(f,s) = sum(S(s,:) .* b(f,:)) / m;
        if e(f, s) == 0, e(f, s) = 1; end
    end
end

VPN = Final_OV .* CI .* e;
initial_VPN = Final_OV .* CI .* e; 

initial_pathVPN_matrix = zeros(numFunctions, numScenarios);
for s = 1:numScenarios
    for f = 1:numFunctions
        downstreamFuncs = downstreamFunctions{f};
        sumVPN_val = 0;
        for i = 1:length(downstreamFuncs)
            idx = find(strcmp(functions, downstreamFuncs{i}));
            if ~isempty(idx), sumVPN_val = sumVPN_val + VPN(idx, s); end
        end
        initial_pathVPN_matrix(f, s) = sumVPN_val;
    end
end
totalPathVPN = sum(initial_pathVPN_matrix, 2);

VPN = Final_OV .* CI .* e; 

initial_pathVPN_matrix = zeros(numFunctions, numScenarios);
basePsoParams_temp = struct('numFunctions', numFunctions, 'numScenarios', numScenarios, 'downstreamFunctions', {downstreamFunctions}, 'functions', {functions});
initial_downstreamSumVPN = calculateTotalPathVPN(initial_VPN, basePsoParams_temp); 

initial_intrinsicVPN = sum(initial_VPN, 2);

initial_rankingMetric = initial_intrinsicVPN + initial_downstreamSumVPN;


figure;

h_graph_vpn = plot(G, 'Layout', 'layered', 'LineWidth', 1.5, 'ArrowSize', 14);
title('Path Graph with Propagated Risk (Downstream Sum VPN) Weighted Connections', 'FontSize', 12);

h_graph_vpn.MarkerSize = 12;
h_graph_vpn.NodeColor = nodeColors;
h_graph_vpn.EdgeAlpha = 0.7;
h_graph_vpn.NodeLabel = {}; 

labelOffset = 0.2;
for i = 1:numFunctions
    text(h_graph_vpn.XData(i) + labelOffset, h_graph_vpn.YData(i) + labelOffset, functions{i}, ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
end

edges_graph = G.Edges; 
edgeWeights = zeros(height(edges_graph), 1);
for i = 1:height(edges_graph)
    sourceNodeName = edges_graph.EndNodes{i, 1};
    sourceIdx = find(strcmp(functions, sourceNodeName));
    if ~isempty(sourceIdx)
        edgeWeights(i) = initial_downstreamSumVPN(sourceIdx); 
    end
end

minWeight = min(edgeWeights);
maxWeight = max(edgeWeights);
if abs(maxWeight - minWeight) < 1e-6 
    edgeThickness = ones(length(edgeWeights), 1) * 2; 
else
    edgeThickness = 1 + (edgeWeights - minWeight) / (maxWeight - minWeight) * 7; 
end

edgeColors_vis = [
    60, 179, 113;  
    147, 112, 219;  
    255, 140, 0;    
    0, 0, 128;      
    139, 0, 0      
]/255; 
numVisColors = size(edgeColors_vis, 1);

scaledEdgeWeights = (edgeWeights - minWeight); 
maxScaled = max(scaledEdgeWeights);
if maxScaled > 1e-6 
    scaledEdgeWeights = scaledEdgeWeights / maxScaled; 
else
    scaledEdgeWeights = zeros(size(scaledEdgeWeights)); 
end
edgeColorIdx = round(scaledEdgeWeights * (numVisColors - 1)) + 1; 
edgeColorIdx = min(max(edgeColorIdx, 1), numVisColors); 
edgeColorMap = edgeColors_vis(edgeColorIdx, :);

h_graph_vpn.EdgeColor = edgeColorMap;
h_graph_vpn.LineWidth = edgeThickness;

colorLegend = axes('Position', [0.832, 0.13, 0.05, 0.8]); 
set(colorLegend, 'Color', 'none', 'XColor', 'none', 'YColor', 'none'); 
colormap(colorLegend, edgeColors_vis); 
caxis(colorLegend, [0, 1]); 

cb = colorbar('peer', colorLegend, 'Ticks', [0, 1], 'TickLabels', {'Low VPN', 'High VPN'}, 'FontSize', 10);

set(gcf, 'Position', [100, 100, 1000, 700]);
colorbarPosition = get(cb, 'Position');
colorbarPosition(3) = 0.02; 
colorbarPosition(4) = 0.3; 
set(cb, 'Position', colorbarPosition);


initial_b = b;
initial_Final_OV = Final_OV;
initial_CI = CI;
initial_e = e;
initial_downstreamSumVPN_stored = initial_downstreamSumVPN;
initial_rankingMetric_stored = initial_rankingMetric; 
fprintf('--------------- Phase 1 ---------------\n');
fprintf('Initial Max Ranking Metric (Intrinsic + Downstream): %.3f\n', max(initial_rankingMetric));
fprintf('Initial Sum Ranking Metric (Intrinsic + Downstream): %.3f\n', sum(initial_rankingMetric));
fprintf('(For reference) Initial Max Propagated Risk (Downstream Sum VPN): %.3f\n', max(initial_downstreamSumVPN));
fprintf('(For reference) Initial Sum Propagated Risk (Downstream Sum VPN): %.3f\n', sum(initial_downstreamSumVPN));

%% Phase 2

fprintf('\n--------------- Phase 2 ---------------');
% Define Mitigation Strategies

strategies = struct();
numStrategies = 7;
feasibilityMap = containers.Map({'Easy', 'Medium', 'Difficult'}, {1.0, 0.7, 0.3}); 
weightMap = containers.Map({'High', 'Medium'}, {1.0, 0.7}); 
susImpactMap = containers.Map({'Low impact', 'Medium impact', 'High impact'}, {0.5, 1.0, 1.5}); 

strategies(1).Name = 'Worker-centric wearable design';
strategies(1).RelatedSPC = {'WC', 'RA'};
strategies(1).RelatedFunctions = {'F-1', 'F-3', 'F-7', 'F-8', 'F-10', 'F-13'};
strategies(1).FeasibilityScore = feasibilityMap('Difficult');
strategies(1).WeightScore = weightMap('High');
strategies(1).SustainabilityScores = [susImpactMap('Low impact'), susImpactMap('Medium impact'), susImpactMap('High impact')];
strategies(1).Effect.Type = 'Improve Worker Satisfaction';
strategies(1).Effect.TargetParameter = 'Final_OV';
strategies(1).Effect.Magnitude = 0.80;
strategies(1).Applied = false;

strategies(2).Name = 'Low-power sensors';
strategies(2).RelatedSPC = {'WP'};
strategies(2).RelatedFunctions = {'F-3'};
strategies(2).FeasibilityScore = feasibilityMap('Difficult');
strategies(2).WeightScore = weightMap('Medium');
strategies(2).SustainabilityScores = [susImpactMap('High impact'), susImpactMap('Medium impact'), susImpactMap('Low impact')];
strategies(2).Effect.Type = 'Improve Environmental Aspects';
strategies(2).Effect.TargetParameter = 'b';
strategies(2).Effect.RelatedSPCIndices = [1, 3];
strategies(2).Effect.Magnitude = 0.90;
strategies(2).Applied = false;

strategies(3).Name = 'Wearable maintenance program';
strategies(3).RelatedSPC = {'WP', 'RA'};
strategies(3).RelatedFunctions = {'F-3', 'F-4', 'F-6', 'F-9', 'F-11', 'F-15'};
strategies(3).FeasibilityScore = feasibilityMap('Medium');
strategies(3).WeightScore = weightMap('High');
strategies(3).SustainabilityScores = [susImpactMap('High impact'), susImpactMap('Medium impact'), susImpactMap('Low impact')];
strategies(3).Effect.Type = 'Improve the Devices Lifetime';
strategies(3).Effect.TargetParameter = 'b';
strategies(3).Effect.RelatedSPCIndices = [1, 3];
strategies(3).Effect.Magnitude = 0.75;
strategies(3).Applied = false;

strategies(4).Name = 'Wellness monitoring and breaks';
strategies(4).RelatedSPC = {'WC'};
strategies(4).RelatedFunctions = {'F-8', 'F-9', 'F-13'};
strategies(4).FeasibilityScore = feasibilityMap('Easy');
strategies(4).WeightScore = weightMap('High');
strategies(4).SustainabilityScores = [susImpactMap('Low impact'), susImpactMap('Medium impact'), susImpactMap('High impact')];
strategies(4).Effect.Type = 'Improve Worker Satisfaction';
strategies(4).Effect.TargetParameter = 'Final_OV';
strategies(4).Effect.Magnitude = 0.75;
strategies(4).Applied = false;

strategies(5).Name = 'Tool sharing optimization';
strategies(5).RelatedSPC = {'WP', 'RA'};
strategies(5).RelatedFunctions = {'F-2', 'F-9', 'F-14'};
strategies(5).FeasibilityScore = feasibilityMap('Easy');
strategies(5).WeightScore = weightMap('Medium');
strategies(5).SustainabilityScores = [susImpactMap('High impact'), susImpactMap('Medium impact'), susImpactMap('Low impact')];
strategies(5).Effect.Type = ', Reducing Idle Resources';
strategies(5).Effect.TargetParameter = 'b';
strategies(5).Effect.RelatedSPCIndices = [1, 3];
strategies(5).Effect.Magnitude = 0.85;
strategies(5).Applied = false;

strategies(6).Name = 'Reusable tool systems';
strategies(6).RelatedSPC = {'RA', 'WP'};
strategies(6).RelatedFunctions = {'F-2', 'F-14'};
strategies(6).FeasibilityScore = feasibilityMap('Medium');
strategies(6).WeightScore = weightMap('Medium');
strategies(6).SustainabilityScores = [susImpactMap('High impact'), susImpactMap('High impact'), susImpactMap('Low impact')];
strategies(6).Effect.Type = 'Improve Resource Management';
strategies(6).Effect.TargetParameter = 'b';
strategies(6).Effect.RelatedSPCIndices = [3, 1];
strategies(6).Effect.Magnitude = 0.85;
strategies(6).Applied = false;

strategies(7).Name = 'Digital inventory management';
strategies(7).RelatedSPC = {'RA'};
strategies(7).RelatedFunctions = {'F-2', 'F-5', 'F-11', 'F-14'};
strategies(7).FeasibilityScore = feasibilityMap('Medium');
strategies(7).WeightScore = weightMap('High');
strategies(7).SustainabilityScores = [susImpactMap('Medium impact'), susImpactMap('High impact'), susImpactMap('Medium impact')];
strategies(7).Effect.Type = 'Improve Resource Management';
strategies(7).Effect.TargetParameter = 'Final_OV';
strategies(7).Effect.Magnitude = 0.70;
strategies(7).Applied = false;

fprintf('\n--- Defining Mitigation Strategies ---\n');

for s_idx = 1:numStrategies
    if ~isfield(strategies(s_idx),'Effect') || ~isfield(strategies(s_idx).Effect, 'TargetParameter') || ~isfield(strategies(s_idx).Effect, 'Magnitude') 
        error('Incomplete effect definition for Strategy %d: Missing TargetParameter or Magnitude.', s_idx);
    end
    if strcmp(strategies(s_idx).Effect.TargetParameter, 'b') && (~isfield(strategies(s_idx).Effect, 'RelatedSPCIndices') || isempty(strategies(s_idx).Effect.RelatedSPCIndices))
         error('Strategy %d targeting "b" is missing required RelatedSPCIndices.', s_idx);
    end
end

%% Phase 3

fprintf('\n--------------- Phase 3 ---------------\n');


maxMitigationSteps = 4; 
% --- Decision Maker Weights ---
sustainabilityWeights = [0.334, 0.333, 0.333];
if abs(sum(sustainabilityWeights) - 1) > 1e-6, error('Sustainability weights sum');
end
fprintf('Sustainability Weights: Environmental=%.2f, Economic=%.2f, Social=%.2f\n', sustainabilityWeights(1), sustainabilityWeights(2), sustainabilityWeights(3));

% --- PSO Settings for runPSO Function ---
basePsoParams = struct();
basePsoParams.numFunctions = numFunctions;
basePsoParams.numScenarios = numScenarios;
basePsoParams.downstreamFunctions = downstreamFunctions;
basePsoParams.functions = functions; 

basePsoParams.pso.nParticles = 50;
basePsoParams.pso.nIterations = 500; 
basePsoParams.pso.w_start = 0.9;
basePsoParams.pso.w_end = 0.4;
basePsoParams.pso.c1 = 0.9;
basePsoParams.pso.c2 = 1.5;
basePsoParams.pso.mutation_prob = 0.05;
basePsoParams.pso.maxVelFactor = 0.4;
basePsoParams.pso.initSpreadFactor = 0.6;
basePsoParams.pso.lowerBoundFactor = 0.5;
basePsoParams.pso.upperBoundFactor = 1.5;
basePsoParams.pso.fitnessType = 'max'; 

mitigationHistory_ActualMaxPropRisk = zeros(maxMitigationSteps + 1, 1);
mitigationHistory_ActualSumPropRisk = zeros(maxMitigationSteps + 1, 1);
appliedStrategyHistory = cell(maxMitigationSteps, 1);
costHistory = zeros(maxMitigationSteps + 1, 1);

current_b = initial_b;
current_Final_OV = initial_Final_OV;
current_CI = initial_CI;
current_e = initial_e;
current_VPN = initial_VPN; 
current_downstreamSumVPN = initial_downstreamSumVPN_stored; 
current_rankingMetric = initial_rankingMetric_stored;

mitigationHistory_ActualMaxPropRisk(1) = max(current_downstreamSumVPN); 
mitigationHistory_ActualSumPropRisk(1) = sum(current_downstreamSumVPN); 
costHistory(1) = 0;
cumulativeMitigationCost = 0;

targetSelectionCount = zeros(numFunctions, 1); 

for s_idx_reset = 1:numStrategies, strategies(s_idx_reset).Applied = false;
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN MITIGATION LOOP ~~~~~~~~~~~~~~~~~~~~~~~~~~
for step = 1:maxMitigationSteps
    fprintf('\n======= Mitigation Step %d =======\n', step);

    strategyAppliedThisStep = false;
    triedTargetsInThisStep = false(numFunctions, 1); 
    selectedTargetForStep = -1; 
    appliedStrategyNameForStep = 'None'; 

    while true 

        targetFuncIndex = -1; 
        maxRankingMetric_candidate = -Inf;

        [sortedRankMetric, sortedIndices] = sort(current_rankingMetric, 'descend');

        for i = 1:numFunctions
            candidateIndex = sortedIndices(i);
            candidateRankMetric = sortedRankMetric(i); 

            isGloballyEligible = candidateRankMetric > 1e-6 && targetSelectionCount(candidateIndex) < 2;
            isEligibleForThisTry = isGloballyEligible && ~triedTargetsInThisStep(candidateIndex);

            if isEligibleForThisTry
                hasAvailableStrategy = false;
                for s_idx_check = 1:numStrategies
                    if any(strcmp(functions{candidateIndex}, strategies(s_idx_check).RelatedFunctions)) && ~strategies(s_idx_check).Applied
                        hasAvailableStrategy = true;
                        break; 
                    end
                end

                if hasAvailableStrategy
                    targetFuncIndex = candidateIndex;
                    maxRankingMetric_candidate = candidateRankMetric;
                    fprintf('  Considering target candidate for Step %d: %s (Ranking Metric: %.3f, Global Selections: %d)\n', ...
                            step, functions{targetFuncIndex}, maxRankingMetric_candidate, targetSelectionCount(targetFuncIndex));
                    break; 
                else
                     triedTargetsInThisStep(candidateIndex) = true;
                end
            end
        end 

        if targetFuncIndex == -1
            fprintf('No further eligible target functions with available strategies found in Step %d.\n', step);
            break; 
        end

        triedTargetsInThisStep(targetFuncIndex) = true;

        targetSelectionCount(targetFuncIndex) = targetSelectionCount(targetFuncIndex) + 1;
        fprintf('>>> Targeting function "%s" for mitigation attempt (Current Global Selection Count: %d)\n', ...
                 functions{targetFuncIndex}, targetSelectionCount(targetFuncIndex));


        possibleStrategies_for_target = []; 
        fprintf('  Evaluating potential strategies for current target "%s":\n', functions{targetFuncIndex});

        for s_idx = 1:numStrategies
            isApplicable = any(strcmp(functions{targetFuncIndex}, strategies(s_idx).RelatedFunctions));
            isApplied = strategies(s_idx).Applied;

            if isApplicable && ~isApplied
                 fprintf('    Simulating Strategy %d: %s...\n', s_idx, strategies(s_idx).Name);
                 strategyInfo = struct(); strategyInfo.StrategyIndex = s_idx; strategyInfo.Name = strategies(s_idx).Name;
                 susCost = dot(sustainabilityWeights, strategies(s_idx).SustainabilityScores);
                 strategyInfo.TotalCost = susCost + strategies(s_idx).FeasibilityScore;
                 temp_b = current_b; temp_Final_OV = current_Final_OV; temp_CI = current_CI;
                 affectedFuncIndices = find(ismember(functions, strategies(s_idx).RelatedFunctions));
                 effect = strategies(s_idx).Effect; validEffect = true;
                 try
                      switch effect.TargetParameter
                          case 'b', temp_b(affectedFuncIndices, effect.RelatedSPCIndices) = max(0, temp_b(affectedFuncIndices, effect.RelatedSPCIndices) * effect.Magnitude);
                          case 'Final_OV', temp_Final_OV(affectedFuncIndices) = max(0, temp_Final_OV(affectedFuncIndices) * effect.Magnitude);
                          case 'CI', temp_CI(affectedFuncIndices) = max(0, temp_CI(affectedFuncIndices) * effect.Magnitude);
                          otherwise, warning('Unknown TargetParameter: %s', effect.TargetParameter); validEffect = false;
                      end
                 catch ME, validEffect = false; warning('Sim Error strat %d: %s', s_idx, ME.message); end
                 if ~validEffect, continue; end
                 temp_e = zeros(numFunctions, numScenarios);
                 for scenario_idx = 1:numScenarios, for f_idx = 1:numFunctions
                     temp_e(f_idx,scenario_idx) = sum(S(scenario_idx,:) .* temp_b(f_idx,:)) / m;
                     if temp_e(f_idx, scenario_idx) == 0, temp_e(f_idx, scenario_idx) = 1; end
                 end; end
                 temp_VPN = temp_Final_OV .* temp_CI .* temp_e;
                 strategyInfo.temp_VPN_state = temp_VPN;
                 fprintf('      Running PSO simulation for Strategy %d...\n', s_idx); tic;
                 [psoFitnessAfterStrategy, ~] = runPSO(temp_VPN, basePsoParams); psoTime = toc;
                 fprintf('      PSO finished (%.2f sec). Predicted Best Overall Propagated Risk (%s): %.4f\n', psoTime, basePsoParams.pso.fitnessType, psoFitnessAfterStrategy);
                 strategyInfo.Predicted_PSO_Fitness = psoFitnessAfterStrategy;
                 if strcmp(basePsoParams.pso.fitnessType, 'max')
                     currentFitness_for_Benefit = max(current_downstreamSumVPN);
                 else, currentFitness_for_Benefit = sum(current_downstreamSumVPN); end
                 benefit = currentFitness_for_Benefit - psoFitnessAfterStrategy;
                 costForScore = max(strategyInfo.TotalCost, 1e-6);
                 strategyInfo.Score = (max(0, benefit) * strategies(s_idx).WeightScore) / costForScore;
                 possibleStrategies_for_target = [possibleStrategies_for_target, strategyInfo];
            end
         end 


        bestStratInfoForTarget = struct();
        foundSuitableStrategyForTarget = false;
        if isempty(possibleStrategies_for_target)
            fprintf('  No applicable/unapplied strategies evaluated for target "%s". Moving to check next target.\n', functions{targetFuncIndex});
        else
            fprintf('  Selecting best strategy for "%s" based on score AND >=5%% target''s OWN AVG VPN improvement:\n', functions{targetFuncIndex});
            [~, sortedScoreIndices] = sort([possibleStrategies_for_target.Score], 'descend');

            for i = 1:length(sortedScoreIndices) 
                 currentLocalIndex = sortedScoreIndices(i); tempStratInfo = possibleStrategies_for_target(currentLocalIndex);
                 score = tempStratInfo.Score; name = tempStratInfo.Name;
                 fprintf('    Checking Candidate Strategy: "%s" (Score=%.4f)\n', name, score);
                 if score <= 0, fprintf('      Skipped: Score <= 0.\n'); continue; end

                 baseline_target_avg_VPN = mean(current_VPN(targetFuncIndex, :));
                 predicted_target_avg_VPN = mean(tempStratInfo.temp_VPN_state(targetFuncIndex, :));
                 meetsThreshold = false; percent_improvement = 0;
                 if baseline_target_avg_VPN > 1e-6
                     improvement = baseline_target_avg_VPN - predicted_target_avg_VPN;
                     percent_improvement = (improvement / baseline_target_avg_VPN) * 100;
                     if percent_improvement >= 15.0 
                         meetsThreshold = true;
                     end
                 elseif predicted_target_avg_VPN < (baseline_target_avg_VPN - 1e-6)
                      meetsThreshold = true; percent_improvement = 100.0;
                 end
                 fprintf('      Target OWN AVG VPN: Before=%.4f, After(Predicted)=%.4f. Improvement=%.2f%%\n', baseline_target_avg_VPN, predicted_target_avg_VPN, percent_improvement);

                 if meetsThreshold
                     fprintf('      PASSED 15%% threshold check (Own Avg VPN).\n');
                     bestStratInfoForTarget = tempStratInfo;
                     foundSuitableStrategyForTarget = true;
                     fprintf('    >>> Strategy "%s" selected for application this step.\n', name);
                     break; 
                 else
                     fprintf('      FAILED 5%% threshold check (Own Avg VPN).\n');
                 end
            end 

            if foundSuitableStrategyForTarget 
                selectedStratGlobalIndex = bestStratInfoForTarget.StrategyIndex;
                fprintf('  Applying Strategy: "%s" (Target: %s)\n', bestStratInfoForTarget.Name, functions{targetFuncIndex});
                affectedFuncIndices_perm = find(ismember(functions, strategies(selectedStratGlobalIndex).RelatedFunctions));
                effect_perm = strategies(selectedStratGlobalIndex).Effect;
                try 
                     switch effect_perm.TargetParameter
                        case 'b', current_b(affectedFuncIndices_perm, effect_perm.RelatedSPCIndices) = max(0, current_b(affectedFuncIndices_perm, effect_perm.RelatedSPCIndices) * effect_perm.Magnitude);
                        case 'Final_OV', current_Final_OV(affectedFuncIndices_perm) = max(0, current_Final_OV(affectedFuncIndices_perm) * effect_perm.Magnitude);
                        case 'CI', current_CI(affectedFuncIndices_perm) = max(0, current_CI(affectedFuncIndices_perm) * effect_perm.Magnitude);
                        otherwise, warning('PERSISTENT Error applying effect: Unknown TargetParameter: %s', effect_perm.TargetParameter);
                     end

                     strategies(selectedStratGlobalIndex).Applied = true; 
                     strategyAppliedThisStep = true;
                     selectedTargetForStep = targetFuncIndex;
                     appliedStrategyNameForStep = bestStratInfoForTarget.Name;
                     fprintf('    Applied effect permanently.\n');

                     for scenario_idx = 1:numScenarios, for f_idx = 1:numFunctions
                         current_e(f_idx,scenario_idx) = sum(S(scenario_idx,:) .* current_b(f_idx,:)) / m;
                         if current_e(f_idx, scenario_idx) == 0, current_e(f_idx, scenario_idx) = 1; end
                     end; end
                     current_VPN = current_Final_OV .* current_CI .* current_e;
                     current_downstreamSumVPN = calculateTotalPathVPN(current_VPN, basePsoParams);
                     current_intrinsicVPN = sum(current_VPN, 2); 
                     current_rankingMetric = current_intrinsicVPN + current_downstreamSumVPN; 

                     cumulativeMitigationCost = cumulativeMitigationCost + bestStratInfoForTarget.TotalCost;
                     break; 

                catch ME 
                    warning('Error applying permanent effect strat %d: %s. Trying next target.', selectedStratGlobalIndex, ME.message);
                    strategyAppliedThisStep = false; 
                end 
            else 
                fprintf('  No strategy met 15%% threshold for target "%s". Moving to check next target function.\n', functions{targetFuncIndex});
            end 

        end 
    end 

    if strategyAppliedThisStep
        mitigationHistory_ActualMaxPropRisk(step + 1) = max(current_downstreamSumVPN);
        mitigationHistory_ActualSumPropRisk(step + 1) = sum(current_downstreamSumVPN);
        costHistory(step + 1) = cumulativeMitigationCost;
        appliedStrategyHistory{step} = sprintf('%s (Target: %s)', appliedStrategyNameForStep, functions{selectedTargetForStep});
         fprintf('--- State After Step %d (Strategy Applied) ---\n', step);
         fprintf('  Applied: %s\n', appliedStrategyHistory{step});
         fprintf('  Max Propagated Risk : %.4f\n', mitigationHistory_ActualMaxPropRisk(step+1));
         fprintf('  Sum Propagated Risk : %.4f\n', mitigationHistory_ActualSumPropRisk(step+1));
         fprintf('  Cumulative Cost: %.3f\n', costHistory(step+1));
    else 
        mitigationHistory_ActualMaxPropRisk(step + 1) = mitigationHistory_ActualMaxPropRisk(step);
        mitigationHistory_ActualSumPropRisk(step + 1) = mitigationHistory_ActualSumPropRisk(step);
        costHistory(step + 1) = costHistory(step);
        appliedStrategyHistory{step} = 'None (No target/strategy met criteria this step)';
        fprintf('--- State After Step %d (No Strategy Applied) ---\n', step);
        fprintf('  Max Propagated Risk : %.4f\n', mitigationHistory_ActualMaxPropRisk(step+1));
        fprintf('  Sum Propagated Risk : %.4f\n', mitigationHistory_ActualSumPropRisk(step+1));
        fprintf('  Cumulative Cost: %.3f\n', costHistory(step+1));
    end

     if sum([strategies.Applied]) == numStrategies, fprintf('Stopping: All strategies applied.\n'); break; end
     if max(current_downstreamSumVPN) < 1e-2, fprintf('Stopping: Max Propagated Risk very low (%.4f).\n', max(current_downstreamSumVPN)); break; end

end 
    
final_downstreamSumVPN = current_downstreamSumVPN; 
final_intrinsicVPN = current_intrinsicVPN;        
final_rankingMetric = final_intrinsicVPN + final_downstreamSumVPN;

actual_steps_run = step; 

last_successful_apply_step = find(cellfun(@(x) ~isempty(x) && ~startsWith(x, 'None'), appliedStrategyHistory(1:actual_steps_run)), 1, 'last');
if isempty(last_successful_apply_step)
    finalMitigationStep = 0; 
else
    finalMitigationStep = last_successful_apply_step;
end

actual_steps_run = step; 

last_successful_apply_step = find(cellfun(@(x) ~isempty(x) && ~startsWith(x, 'None'), appliedStrategyHistory(1:actual_steps_run)), 1, 'last');
if isempty(last_successful_apply_step)
    finalMitigationStep = 0; 
else
    finalMitigationStep = last_successful_apply_step;
end


mitigationHistory_ActualMaxPropRisk = mitigationHistory_ActualMaxPropRisk(1:actual_steps_run+1);
mitigationHistory_ActualSumPropRisk = mitigationHistory_ActualSumPropRisk(1:actual_steps_run+1);
costHistory = costHistory(1:actual_steps_run+1);
appliedStrategyHistory = appliedStrategyHistory(1:actual_steps_run);


final_downstreamSumVPN = current_downstreamSumVPN; 

fprintf('\n======= Final Results =======\n');
fprintf('Completed %d mitigation loop iterations.\n', actual_steps_run);
fprintf('Last step where a strategy was successfully applied: %d\n', finalMitigationStep);
fprintf('Initial Max Propagated Risk : %.4f\n', mitigationHistory_ActualMaxPropRisk(1));
fprintf('Final Max Propagated Risk   : %.4f\n', mitigationHistory_ActualMaxPropRisk(end));
fprintf('Initial Sum Propagated Risk : %.4f\n', mitigationHistory_ActualSumPropRisk(1));
fprintf('Final Sum Propagated Risk   : %.4f\n', mitigationHistory_ActualSumPropRisk(end));
fprintf('Sequence of Applied/Attempted Strategies:\n');
for i = 1:actual_steps_run
    fprintf('  Step %d: %s\n', i, appliedStrategyHistory{i});
end
fprintf('Final Cumulative Cost: %.3f\n', costHistory(end));

%% Final Plotting


pct_reduction_ranking = (initial_rankingMetric_stored - final_rankingMetric) ./ initial_rankingMetric_stored * 100;
pct_reduction_ranking(initial_rankingMetric_stored == 0 | isnan(initial_rankingMetric_stored)) = 0;

final_rankingMetric_safe = final_rankingMetric;
final_rankingMetric_safe(isnan(final_rankingMetric_safe) | isinf(final_rankingMetric_safe)) = 0;
initial_rankingMetric_safe = initial_rankingMetric_stored;
initial_rankingMetric_safe(isnan(initial_rankingMetric_safe) | isinf(initial_rankingMetric_safe)) = 0;


figure;
barData = [initial_rankingMetric_safe(:), final_rankingMetric_safe(:)];
hBar = bar(barData); 
ax1 = gca; 

set(ax1, 'XTick', 1:numFunctions, 'XTickLabel', functions, 'XTickLabelRotation', 45);
ylabel('Ranking Metric (Intrinsic + Downstream)', 'FontSize', 14); 
title('Ranking Metric Comparison: Initial vs. Final Mitigated', 'FontSize', 15);
legend('Initial Ranking Metric', 'Final Ranking Metric', 'Location', 'northwest', 'FontSize', 15); 
grid off; 
set(ax1, 'FontSize', 15); 

hold(ax1, 'on'); 
y_max_overall_with_text = max(ylim(ax1));

for i = 1:numFunctions
    x_pos = i;

    y_pos_bar = max(initial_rankingMetric_safe(i), final_rankingMetric_safe(i));

    if isnan(y_pos_bar) || isinf(y_pos_bar) || y_pos_bar <= 1e-6
         current_ylim = ylim(ax1);
         y_pos_text = current_ylim(1) + 0.03 * (current_ylim(2) - current_ylim(1));
    else
         current_ylim = ylim(ax1);
         y_offset = 0.03 * (current_ylim(2) - current_ylim(1));
         y_pos_text = y_pos_bar + y_offset;
    end

    reduction_val = pct_reduction_ranking(i);

    if isnan(reduction_val) || isinf(reduction_val)
        text_str = 'N/A';
    else
        text_str = sprintf('%.1f%%', reduction_val); 
    end

    text(ax1, x_pos, y_pos_text, text_str, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 15, 'Color', 'k');

    y_max_overall_with_text = max(y_max_overall_with_text, y_pos_text);
end
hold(ax1, 'off'); 

current_ylim = ylim(ax1);
new_ylim_upper = y_max_overall_with_text + 0.05 * (y_max_overall_with_text - current_ylim(1));
new_ylim_upper = max(current_ylim(2), new_ylim_upper);
new_ylim_lower = max(0, current_ylim(1));
ylim(ax1, [new_ylim_lower, new_ylim_upper]); 



figure; 
ax_prog = gca; 

h_max = plot(ax_prog, 0:actual_steps_run, mitigationHistory_ActualMaxPropRisk, '-^k', 'LineWidth', 5, 'MarkerSize', 9, 'DisplayName', 'Max Path VPN');
hold(ax_prog, 'on');
h_sum = plot(ax_prog, 0:actual_steps_run, mitigationHistory_ActualSumPropRisk, '-ob', 'LineWidth', 5, 'MarkerSize', 9, 'DisplayName', 'Sum Path VPN');

y_lims = ylim(ax_prog);
y_range = y_lims(2) - y_lims(1);
y_offset_factor = 0.04; 
y_offset_val = y_range * y_offset_factor;

for step_num = 1:actual_steps_run 
    prev_max = mitigationHistory_ActualMaxPropRisk(step_num);
    curr_max = mitigationHistory_ActualMaxPropRisk(step_num + 1);
    prev_sum = mitigationHistory_ActualSumPropRisk(step_num);
    curr_sum = mitigationHistory_ActualSumPropRisk(step_num + 1);

    perc_impr_max = 0;
    if prev_max > 1e-9 
        improvement_max = prev_max - curr_max;
        perc_impr_max = (improvement_max / prev_max) * 100;
    end

    perc_impr_sum = 0;
    if prev_sum > 1e-9
        improvement_sum = prev_sum - curr_sum;
        perc_impr_sum = (improvement_sum / prev_sum) * 100;
    end

    if abs(perc_impr_max) > 0.1
        text_str_max = sprintf('%.1f%%', perc_impr_max);
        text(ax_prog, step_num, curr_max + y_offset_val, text_str_max, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 15, 'Color', h_max.Color); 
    end

    if abs(perc_impr_sum) > 0.1
        text_str_sum = sprintf('%.1f%%', perc_impr_sum);
        text(ax_prog, step_num, curr_sum - y_offset_val, text_str_sum, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
             'FontSize', 15, 'Color', h_sum.Color);
    end
end

hold(ax_prog, 'off'); 

title('Mitigation Progress: Propagated Risk Reduction Over Steps', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Mitigation Step Applied', 'FontSize', 18);
ylabel('Path VPN Value', 'FontSize', 18);
legend('Location', 'best', 'FontSize', 16); grid on;
xlim(ax_prog, [0, max(1, actual_steps_run)]);
xticks(ax_prog, 0:1:actual_steps_run);
set(ax_prog, 'FontSize', 16);

final_y_lims = ylim(ax_prog);
ylim(ax_prog, [final_y_lims(1) - y_offset_val*0.5, final_y_lims(2) + y_offset_val*0.5]); 


if finalMitigationStep > 0 
    figure('Name', 'Strategy Affect Imagesc');

    successful_indices = find(cellfun(@(x) ~isempty(x) && ~startsWith(x, 'None'), appliedStrategyHistory(1:finalMitigationStep)));
    num_successful_steps = length(successful_indices);
    
    if num_successful_steps > 0
        stratVsFuncMatrix = zeros(num_successful_steps, numFunctions);
        appliedStratNames = cell(num_successful_steps, 1);
        
        for plot_idx = 1:num_successful_steps
            step_idx = successful_indices(plot_idx); 
            full_hist_entry = appliedStrategyHistory{step_idx};
            
            name_parts = regexp(full_hist_entry, '^(.*) \(Target: F-\d+\)$', 'tokens');
            if ~isempty(name_parts)
                strat_name = strtrim(name_parts{1}{1});
            else
                strat_name = full_hist_entry; 
            end

            stratIdx = find(strcmp({strategies.Name}, strat_name), 1);

            if ~isempty(stratIdx)
                 affectedFuncNames = strategies(stratIdx).RelatedFunctions;
                 affectedIndices = find(ismember(functions, affectedFuncNames));
                 stratVsFuncMatrix(plot_idx, affectedIndices) = 1; 
                 appliedStratNames{plot_idx} = sprintf('S%d:%s', step_idx, strat_name(1:min(15, length(strat_name))));
             else
                 appliedStratNames{plot_idx} = sprintf('S%d: Unknown', step_idx);
             end
        end

        imagesc(stratVsFuncMatrix);
        colormap([1 1 1; 0.1 0.1 0.8]);

        ax = gca;
        xLim = xlim(ax); yLim = ylim(ax);
        for k = 1.5 : 1 : (numFunctions - 0.5), line([k k], yLim, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); end
        for k = 1.5 : 1 : (num_successful_steps - 0.5), line(xLim, [k k], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); end

        ax.XTick = 1:numFunctions;
        ax.XTickLabel = functions;
        ax.XTickLabelRotation = 90;
        ax.YTick = 1:num_successful_steps;
        ax.YTickLabel = appliedStratNames;
        ax.XAxis.TickLength = [0 0];
        ax.YAxis.TickLength = [0 0];
        xlim(ax, [0.5, numFunctions + 0.5]);
        ylim(ax, [0.5, num_successful_steps + 0.5]); 
        title('Functions Affected by Each Successfully Applied Strategy', 'FontSize', 14);
        xlabel('Functions', 'FontSize', 15);
        ylabel('Applied Strategy', 'FontSize', 15);
        colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Not Affected', 'Affected'});
        box on;
        set(ax, 'FontSize', 10);
    else
         disp('Plot 6 skipped: Although steps were run, no strategies were successfully applied.');
    end
else
    disp('Plot 6 skipped: No mitigation steps completed or no strategies successfully applied.');
end


%% Helper Function 
function totalPathVPN = calculateTotalPathVPN(vpnMatrix, baseParams)
    numF = baseParams.numFunctions;
    numS = baseParams.numScenarios;
    dsFuncs = baseParams.downstreamFunctions;
    funcs = baseParams.functions;

    pathVPN_matrix = zeros(numF, numS);
    for s = 1:numS
        for f = 1:numF
            downstreamNames = dsFuncs{f}; sumVPN_p = 0;
            for j = 1:length(downstreamNames)
                idx = find(strcmp(funcs, downstreamNames{j}));
                if ~isempty(idx) && idx <= size(vpnMatrix, 1) && s <= size(vpnMatrix, 2)
                    sumVPN_p = sumVPN_p + vpnMatrix(idx, s);
                end
            end
            pathVPN_matrix(f, s) = sumVPN_p;
        end
    end
    totalPathVPN = sum(pathVPN_matrix, 2);
end

function [bestFitness, bestVPNVector] = runPSO(initialVPNGuess, baseParameters)
    numFunctions = baseParameters.numFunctions;
    numScenarios = baseParameters.numScenarios;
    nParticles = baseParameters.pso.nParticles;
    nIterations = baseParameters.pso.nIterations;
    w_start = baseParameters.pso.w_start;
    w_end = baseParameters.pso.w_end;
    c1 = baseParameters.pso.c1;
    c2 = baseParameters.pso.c2;
    mutation_prob = baseParameters.pso.mutation_prob;
    maxVelFactor = baseParameters.pso.maxVelFactor;
    initSpreadFactor = baseParameters.pso.initSpreadFactor;
    lowerBoundFactor = baseParameters.pso.lowerBoundFactor;
    upperBoundFactor = baseParameters.pso.upperBoundFactor;
    fitnessType = baseParameters.pso.fitnessType;

    nDims = numFunctions * numScenarios;
    if ~isequal(size(initialVPNGuess), [numFunctions, numScenarios])
        error('runPSO: initialVPNGuess size mismatch. Expected [%d, %d], got [%d, %d]', ...
              numFunctions, numScenarios, size(initialVPNGuess,1), size(initialVPNGuess,2));
    end
    baseVPN_vector = reshape(initialVPNGuess, 1, []);

    if any(isnan(baseVPN_vector)) || any(isinf(baseVPN_vector))
       warning('runPSO: NaN or Inf found in initialVPNGuess. Replacing with zeros.');
       baseVPN_vector(isnan(baseVPN_vector) | isinf(baseVPN_vector)) = 0;
    end

    particles = zeros(nParticles, nDims);
    velocity = zeros(nParticles, nDims); 

    lowerBound = lowerBoundFactor * abs(baseVPN_vector); 
    upperBound = upperBoundFactor * abs(baseVPN_vector);
    upperBound = max(upperBound, lowerBound + 1e-6); 

    for i = 1:nParticles
        noise = initSpreadFactor * (2 * rand(1, nDims) - 1) .* abs(baseVPN_vector); 
        particles(i, :) = baseVPN_vector + noise;
        particles(i, :) = max(min(particles(i, :), upperBound), lowerBound);
        particles(i, particles(i,:) < 0) = 0; 
    end
    particles(1, :) = max(min(baseVPN_vector, upperBound), lowerBound);
    particles(1, particles(1,:) < 0) = 0;

    range = upperBound - lowerBound;
    velocity = 0.1 * (rand(nParticles, nDims) - 0.5) .* range; 

    pBestPosition = particles;
    pBestScore = inf(nParticles, 1); 
    gBestScore = inf;
    gBestPosition = particles(1, :); 

    for i = 1:nParticles
        particleVPN = reshape(particles(i, :), numFunctions, numScenarios);
        if size(particleVPN, 1) ~= numFunctions || size(particleVPN, 2) ~= numScenarios
             warning('runPSO Initial Score Calc: Particle VPN size mismatch for particle %d. Skipping.', i);
             continue; 
        end
        particleTotalPathVPN = calculateTotalPathVPN(particleVPN, baseParameters);

        if any(isnan(particleTotalPathVPN)) || any(isinf(particleTotalPathVPN))
            warning('runPSO Initial Score Calc: NaN/Inf in particleTotalPathVPN for particle %d. Assigning Inf fitness.', i);
            fitness = inf;
        else
            if strcmp(fitnessType, 'max')
                fitness = max(particleTotalPathVPN);
            else 
                fitness = sum(particleTotalPathVPN);
            end
        end

        pBestScore(i) = fitness;
        if fitness < gBestScore
            gBestScore = fitness;
            gBestPosition = particles(i, :);
        end
    end
    if isinf(gBestScore) && nParticles > 0
       warning('runPSO Initial Score Calc: gBestScore is Inf. Check initial state or fitness calculation.');
       finite_scores = pBestScore(isfinite(pBestScore));
       if ~isempty(finite_scores)
           [min_finite_score, min_idx] = min(finite_scores);
           original_indices = find(isfinite(pBestScore));
           gBestScore = min_finite_score;
           gBestPosition = particles(original_indices(min_idx), :);
       else
            gBestScore = pBestScore(1);
            gBestPosition = particles(1,:);
       end

    end

    w = w_start; 
    for iter = 1:nIterations
        for i = 1:nParticles
        
             particleVPN = reshape(particles(i, :), numFunctions, numScenarios);
             if size(particleVPN, 1) ~= numFunctions || size(particleVPN, 2) ~= numScenarios
                 warning('runPSO Loop: Particle VPN size mismatch for particle %d, iter %d. Skipping fitness.', i, iter);
                 continue; 
             end
             particleTotalPathVPN = calculateTotalPathVPN(particleVPN, baseParameters);

             if any(isnan(particleTotalPathVPN)) || any(isinf(particleTotalPathVPN))
                 fitness = inf; 
             else
                 if strcmp(fitnessType, 'max')
                    fitness = max(particleTotalPathVPN);
                 else
                     fitness = sum(particleTotalPathVPN);
                 end
             end

            if rand < mutation_prob
                 mutationIdx = randi(nDims);
                 noise_scale = max(abs(baseVPN_vector(mutationIdx)), range(mutationIdx) * 0.1);
                 noise = 0.1 * noise_scale * randn;
                 particles(i, mutationIdx) = particles(i, mutationIdx) + noise;
                 particles(i, mutationIdx) = max(min(particles(i, mutationIdx), upperBound(mutationIdx)), lowerBound(mutationIdx));
                 particles(i, particles(i,:) < 0) = 0;
            end

            if fitness < pBestScore(i)
                pBestScore(i) = fitness;
                pBestPosition(i, :) = particles(i, :);
            end
            if fitness < gBestScore
                gBestScore = fitness;
                gBestPosition = particles(i, :);
            end
        end

        maxVelocity = maxVelFactor * range;
        minVelocity = -maxVelocity;

        for i = 1:nParticles
            r1 = rand(1, nDims); r2 = rand(1, nDims);
            cognitive = c1 * r1 .* (pBestPosition(i, :) - particles(i, :));
            social = c2 * r2 .* (gBestPosition - particles(i, :)); 

            velocity(i, :) = w * velocity(i, :) + cognitive + social;
            velocity(i, :) = max(min(velocity(i, :), maxVelocity), minVelocity);

            particles(i, :) = particles(i, :) + velocity(i, :);

            particles(i, :) = max(min(particles(i, :), upperBound), lowerBound);
            particles(i, particles(i,:) < 0) = 0; 
        end

        w = w_start - (w_start - w_end) * (iter / nIterations);

    end 

    bestFitness = gBestScore;
    bestVPNVector = gBestPosition; 
    if isinf(bestFitness)
       warning('runPSO finished with Inf bestFitness. PSO might not have converged or initial state was problematic.');
    end
end