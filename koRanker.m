function [ranked_reactions, range, OFDB_min, OFDB_max] ...
    = koRanker(cobra_model, biomass_reaction_id, target_reaction_id, important_reactions)

cobra_model = changeObjective(cobra_model,biomass_reaction_id);

% Ignore exchange reactions:
exchange_reaction_ids = find(findExcRxns(cobra_model,1));
non_exchange_reaction_ids = setdiff(1:length(cobra_model.rxns), exchange_reaction_ids);

% Optimise target while making biomass at different percentages (from 0 to
% 99%) of maximum...
range = 0:5:99;
[OFDB_min, OFDB_max] = multiFluxVariability(cobra_model, non_exchange_reaction_ids, biomass_reaction_id, target_reaction_id, range);

ranked_reactions = [];


% %Reaction ranking

ranked_reactions = [];

number_fva_points = size(OFDB_min, 2);
percentile_5 = round(number_fva_points * 5/100);
percentile_10 = round(number_fva_points * 10/100);
percentile_25 = round(number_fva_points * 25/100);
percentile_75 = round(number_fva_points * 75/100);

for k = 1:length(non_exchange_reaction_ids)
    % Exclude reactions vital for target production (i.e. those that have
    % a non-zero max or min flux at high target production):
    target_zero = [OFDB_max(k, percentile_5:percentile_10); OFDB_min(k, percentile_5:percentile_10)];

    if any(abs(target_zero) < 1e-6)
        OFDB_min_percentile_25 = OFDB_min(k, 1:percentile_25);
        OFDB_max_percentile_25 = OFDB_max(k, 1:percentile_25);
        OFDB_min_percentile_75 = OFDB_min(k, percentile_75+1:end);
        OFDB_max_percentile_75 = OFDB_max(k, percentile_75+1:end);
        
        if sum(OFDB_min_percentile_25) < 0 && sum(OFDB_max_percentile_25) <= 0
            target_ex_ub_diff = sum(abs(OFDB_min_percentile_75) - abs(OFDB_min_percentile_25));
            area_ub = sum(abs(OFDB_min_percentile_75) - abs(OFDB_max_percentile_75));
        elseif sum(OFDB_min_percentile_25) >= 0 && sum(OFDB_max_percentile_25) > 0
            target_ex_ub_diff = sum(OFDB_max_percentile_75 - OFDB_max_percentile_25);
            area_ub = sum(OFDB_max_percentile_75 - OFDB_min_percentile_75);
        else
            continue
        end
        
        computeTotal = target_ex_ub_diff / area_ub;
        
        if ~isinf(computeTotal)
        reaction_index = non_exchange_reaction_ids(k);
        ranked_reactions = [ranked_reactions; reaction_index, any(important_reactions == reaction_index), computeTotal];
        else
            continue
        end
    end
end

% Sort:
if ~isempty(ranked_reactions)
    ranked_reactions = sortrows(ranked_reactions, [-3 -2]);
    
    % Only those FVA values for ranked_reactions need to be returned.
    % Map indexes of ranked_reactions to non_exchange_reaction_ids:
    fva_reaction_indices = arrayfun(@(x) find(non_exchange_reaction_ids == x, 1), ranked_reactions(:,1));
    OFDB_min = OFDB_min(fva_reaction_indices, :);
    OFDB_max = OFDB_max(fva_reaction_indices, :);
else
    OFDB_min = [];
    OFDB_max = [];
end