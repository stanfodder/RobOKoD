function robokod(model_filename, biomass_reaction_id, target_reaction_id, max_knockouts, results_directory, geneKO, MCT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Make results directory:
if ~exist(fullfile(cd, results_directory), 'file')
    mkdir(results_directory);
end

% Read model:
cobra_model = readCbModel(model_filename);

noKO = 1;
lethal_KO_list = [];
%% Iterative knockout analysis:
%for i = 1:max_knockouts
while noKO <= max_knockouts
    i = noKO;
    % MCT:
    important_reactions = metaboliteCentralityTest(cobra_model, biomass_reaction_id, target_reaction_id);
    
    % FVAp:
    [ranked_reactions, range, optimal_target_defined_biomass_fva_min, optimal_target_defined_biomass_fva_max] ...
        = koRanker(cobra_model, biomass_reaction_id, target_reaction_id, important_reactions);
    
    % Exit iterative loop if no candidate knock-out reactions returned:
    if isempty(ranked_reactions)
       break; 
    end
    
    
    %remove any ranked reactions which might be lethal.
    if ~isempty(lethal_KO_list)

        for kk = 1:length(lethal_KO_list)
        removeRR = find((ranked_reactions(:,1)) == lethal_KO_list(kk));
        ranked_reactions(removeRR,:) = [];
        end
    end
    
    
    % Write results to Excel / plot FVAp distributions:
    writeRankedReactions(cobra_model, ranked_reactions, strcat(results_directory, '/', 'ko_ranked_reactions_', num2str(i), '.xls'), 1);
    fva_plot_directory_name = strcat('fva_plots_', num2str(i));
    fva_plot_directory_path = strcat(results_directory, '/', fva_plot_directory_name);
    
    if ~exist(fullfile(cd, fva_plot_directory_path), 'file')
        mkdir(results_directory, fva_plot_directory_name);
    end
    

    % Optimise biomass while making target at different percentages (from 0 to % 99%) of maximum...
    % (This is purely for plotting).
    [optimal_biomass_defined_target_fva_min, optimal_biomass_defined_target_fva_max] = multiFluxVariability(cobra_model, ranked_reactions(:,1), target_reaction_id, biomass_reaction_id, range);
    optimal_biomass_defined_target_fva_min = fliplr(optimal_biomass_defined_target_fva_min);
    optimal_biomass_defined_target_fva_max = fliplr(optimal_biomass_defined_target_fva_max);
    
    fva_plot(ranked_reactions(:,1), range, optimal_target_defined_biomass_fva_max, optimal_target_defined_biomass_fva_min, optimal_biomass_defined_target_fva_max, optimal_biomass_defined_target_fva_min, cobra_model, fva_plot_directory_path);
    
    
    %Find top ranked reaction
    
    %identify if there is an MCT identified reaction with a positive
    %ranking value
    MCTReaction = find(ranked_reactions(:,2) == 1 & ranked_reactions(:,3) >0);
    
    if ~isempty(MCTReaction) && MCT ==1
    reactionNo = ranked_reactions(MCTReaction(1),1); 
    else 
    reactionNo = ranked_reactions(1,1);
    end
    %KO reactions or gene?
    
    %check reaction has a gene associated
    hasGene = (find(full(cobra_model.rxnGeneMat(reactionNo,:))==1));
    
    if ~isempty(hasGene) && geneKO == 1
        
            geneReactions = [];
        
            for k = 1:length(hasGene)
                reactionList = find(full(cobra_model.rxnGeneMat(:,hasGene(k))==1));
                geneReactions = [geneReactions;reactionList];
            end
            associatedReactions = unique(geneReactions);
        
        %KO all reactions associated with the gene
        koTest_model = cobra_model;
        for k = 1:length(associatedReactions)
            koTest_model = changeRxnBounds(koTest_model, cobra_model.rxns(associatedReactions(k)), 0, 'b');
        end
        
        %test for lethality
        
        FBAsolution = optimizeCbModel(koTest_model,'max');
        if FBAsolution.f == 0
           
            fprintf('Selected KO (%s) was rendered lethal, recomputing KO',char(cobra_model.rxns(reactionNo)));
            lethal_KO_list = [lethal_KO_list;associatedReactions];
            
        else
        
            fprintf( 'Knockout %d: %s with associated genes %s \n', i, char(cobra_model.rxns(reactionNo)),char(cobra_model.genes(hasGene)));
            
            %allow KO to be implemented
            cobra_model = koTest_model;
            noKO = noKO + 1;
        end
        
        
    else
        
        %if gene KO selected, but reaction has no gene association
        if geneKO == 1
            fprintf('Warning: There is no gene association with selected KO, only reaction will be knocked out \n');
        end
        
        %KO the associated reaction in the SBML
        koTest_model = cobra_model;
        koTest_model = changeRxnBounds(koTest_model, cobra_model.rxns(reactionNo), 0, 'b');
        
         %test for lethality
        FBAsolution = optimizeCbModel(koTest_model,'max');
        if FBAsolution.f == 0
            fprintf('Selected KO (%s) was rendered lethal, recomputing KO \n',char(cobra_model.rxns(reactionNo)));
            lethal_KO_list = [lethal_KO_list;reactionNo];
            
        else
            fprintf('Knockout %d: %s \n', i, char(cobra_model.rxns(reactionNo)));
            
            %allow KO to be implemented
            cobra_model = koTest_model;
            noKO = noKO + 1;
        end
        
        
    end

end

% Write updated SBML model including knockouts:
if max_knockouts > 0
    model_output_filepath = strcat(results_directory, '/', strrep(model_filename, '.xml', '_ko.xml'));
    writeCbToSBML(cobra_model, model_output_filepath);
end


%% Overexpression analysis
fprintf('computing overxpression')
[ranked_reactions, range, optimal_target_defined_biomass_fva_min, optimal_target_defined_biomass_fva_max, optimal_biomass_defined_target_fva_min, optimal_biomass_defined_target_fva_max] ...
    = overexpressRanker(cobra_model, biomass_reaction_id, target_reaction_id);

% Write results to Excel / plot FVAp distributions:
writeRankedReactions(cobra_model, ranked_reactions, strcat(results_directory, '/', 'overexpression_ranked_reactions.xls'), 0);

fva_plot_directory_name = 'fva_plots_overexpression';
fva_plot_directory_path = strcat(results_directory, '/', fva_plot_directory_name);

if ~exist(fullfile(cd, fva_plot_directory_path), 'file')
    mkdir(results_directory, fva_plot_directory_name);
end

fva_plot(ranked_reactions(:,1), range, optimal_target_defined_biomass_fva_max, optimal_target_defined_biomass_fva_min, optimal_biomass_defined_target_fva_max, optimal_biomass_defined_target_fva_min, cobra_model, fva_plot_directory_path);

% Output to console:
fprintf( 'Overexpression: %d candidate reactions\n', length(ranked_reactions));
end