function [ranked_reactions, range, OFDB_min, OFDB_max, OBDF_min, OBDF_max] ...
    = overexpressRanker(cobra_model, biomass_reaction_id, target_reaction_id)

% Ignore exchange reactions:
exchange_reaction_ids = find(findExcRxns(cobra_model,1));
non_exchange_reaction_ids = setdiff(1:length(cobra_model.rxns), exchange_reaction_ids);

% Optimise target while making biomass at different percentages (from 0 to
% 99%) of maximum...
range = 0:5:99;
[OFDB_min, OFDB_max] = multiFluxVariability(cobra_model, non_exchange_reaction_ids, biomass_reaction_id, target_reaction_id, range);

% Optimise biomass while making target at different percentages (from 0 to
% 99%) of maximum...
[OBDF_min, OBDF_max] = multiFluxVariability(cobra_model, non_exchange_reaction_ids, target_reaction_id, biomass_reaction_id, range);
OBDF_min = fliplr(OBDF_min);
OBDF_max = fliplr(OBDF_max);

ranked_reactions = [];

number_fva_points = size(OFDB_min,2);
percentile_25 = round(number_fva_points * 25/100);

for k = 1:length(non_exchange_reaction_ids)
    %find lowest variability for optimizing target from 5% target.
    
    if sum(OBDF_max(k,1:percentile_25)')>0
    %target variability when making biomass
    max_OBDF = OBDF_max(k,1:percentile_25)';
    max_OFDB = OFDB_max(k,1:percentile_25)';
    min_OFDB = OFDB_min(k,1:percentile_25)';
    min_max_total = sum(max_OFDB - min_OFDB);
    elseif sum(OBDF_min(k,1:percentile_25)')<0
    max_OBDF = OBDF_min(k,1:percentile_25)';
    max_OFDB = OFDB_min(k,1:percentile_25)';
    min_OFDB = OFDB_max(k,1:percentile_25)';
    min_max_total = sum(abs(min_OFDB) - abs(max_OFDB)); 
    else
    min_max_total = 0;
    end
    
    
    if abs(min_max_total) > 1e-6
        target_ub_diff = sum(abs(max_OFDB) - abs(max_OBDF));
        ranked_reactions = [ranked_reactions; non_exchange_reaction_ids(k), target_ub_diff / min_max_total];
    end
end

% Sort:
if ~isempty(ranked_reactions)
    ranked_reactions = sortrows(ranked_reactions, -2);
    
    % Only those FVA values for ranked_reactions need to be returned.
    % Map indexes of ranked_reactions to non_exchange_reaction_ids:
    fva_reaction_indices = arrayfun(@(x) find(non_exchange_reaction_ids == x, 1), ranked_reactions(:,1));
    OFDB_min = OFDB_min(fva_reaction_indices, :);
    OFDB_max = OFDB_max(fva_reaction_indices, :);
    OBDF_min = OBDF_min(fva_reaction_indices, :);
    OBDF_max = OBDF_max(fva_reaction_indices, :);
else
    OFDB_min = [];
    OFDB_max = [];
    OBDF_min = [];
    OBDF_max = [];
end