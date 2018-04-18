function optKnockSolution = runOptKnock(model_filename, biomass_reaction_id, target_reaction_id, ...
    constraint_reaction_list, constraint_reaction_values, max_knockouts, exclude_exchange_reactions)

    cobra_model = readCbModel(model_filename);
    cobra_model = changeObjective(cobra_model, biomass_reaction_id);

    options.targetRxn = target_reaction_id; 
    options.vMax = 1000; 
    options.numDel = max_knockouts; 
    options.numDelSense = 'L'; 
    constrOpt.rxnList = constraint_reaction_list; 
    constrOpt.values = constraint_reaction_values; 
    constrOpt.sense = 'GE';

    reaction_indices = 1:length(cobra_model.rxns);
    
    % Ignore biomass reaction:
    excluded_reactions = ismember(cobra_model.rxns, biomass_reaction_id);
    
    % Ignore target reaction:
    excluded_reactions = excluded_reactions + ismember(cobra_model.rxns, target_reaction_id);
    
    % Ignore reaction list:
    excluded_reactions = excluded_reactions + ismember(cobra_model.rxns, constrOpt.rxnList);
    
    % Ignore exchange reactions:
    if exclude_exchange_reactions
        excluded_reactions = excluded_reactions + findExcRxns(cobra_model,1);
    end
    
    reaction_indices = setdiff(reaction_indices, find(excluded_reactions));

    optKnockSolution = OptKnock(cobra_model, {cobra_model.rxns{reaction_indices}}, options, constrOpt);
end