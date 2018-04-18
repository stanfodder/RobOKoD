function robustKnockSolutionRxnList = runRobustKnock(model_filename, biomass_reaction_id, target_reaction_id, ...
    constraint_reaction_list, constraint_reaction_values, max_knockouts)

    %results = robustKnock(chemicalInd, objectiveInd, knockoutNum, maxW, optKnockFlag, robustKnockFlag, constraintsList)
    %chemicalInd - chemical to produce
    %objectiveInd - organizms objective
    %knockoutNum - number of knockouts
    %maxW - maximal value of dual variables (higher number will be more accurate but takes more calculation time)
    %optKnockFlag - indicates if optKnock will be calculated 
    %robustFlag - indicates if robustKnock will be calculated 
    %fixedGlu - size of fixed glucose reaction
    
    % Read model:
    cobra_model = readCbModel(model_filename);
    
    % Set parameters:
    robustFlag = 1;
    knockoutNum = max_knockouts;
    optKnockFlag = 0;
    
    % Constants:
    maxWopt = 1000;
    maxWrobust = 1e7;
    
    % Specify reaction indices:
    chemicalInd = find(ismember(cobra_model.rxns, target_reaction_id));
    biomassInd = find(ismember(cobra_model.rxns, biomass_reaction_id));
    organismObjectiveInd = find(ismember(cobra_model.rxns, biomass_reaction_id));

    % Specify constraints:
    constraintsList.reactions = constraint_reaction_list;
    constraintsList.lb = constraint_reaction_values;
    constraintsList.ub = constraint_reaction_values;
    
    if robustFlag
        maxW = maxWrobust;
    else
        maxW = maxWopt;
    end

    robustKnockSolution = robustKnock(cobra_model, chemicalInd, biomassInd, organismObjectiveInd, ...
        knockoutNum, maxW, optKnockFlag, robustFlag, constraintsList);
    
    % TODO: check whether predictions are stored in knockedYindRound or
    % knockedVind:
    robustKnockSolutionRxnList = cobra_model.rxns(robustKnockSolution.Result_OptKnock2.knockedVind);