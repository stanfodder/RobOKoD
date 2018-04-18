function [mins, maxs] = multiFluxVariability(cobra_model, reaction_ids, defined_reaction_id, objective_reaction_id, range)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    num_reactions = length(reaction_ids);
    range_length = size(range, 2);
    mins = zeros(num_reactions, range_length);
    maxs = zeros(num_reactions, range_length);
    
    % BUG? v0.9 had biomass set as objective for both biomass and target cases.
    % Optimal biomass flux is 0.8739, optimal butanol flux is 10,
    % which is the same as glucose uptake.
    %cobra_model = changeObjective(cobra_model, 'Biomass_Ecoli_core_w_GAM');
    cobra_model = changeObjective(cobra_model, defined_reaction_id);
    optimal_defined_solution = optimizeCbModel(cobra_model, 'max');
    optimal_defined_flux = optimal_defined_solution.f;
    
    cobra_model = changeObjective(cobra_model, objective_reaction_id);
    
    for k = 1:range_length
        percentage = range(k);
        cobra_model = changeRxnBounds(cobra_model, defined_reaction_id, (optimal_defined_flux * (percentage / 100)), 'b');
        [min,max] = fluxVariability(cobra_model, 95, 'max', cobra_model.rxns(reaction_ids), [], '0');
        min((abs(min)) < 1e-6) = 0;
        max((abs(max)) < 1e-6) = 0;

        mins(:,k) = min;
        maxs(:,k) = max;
    end
end