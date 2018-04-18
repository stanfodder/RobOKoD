function fva_plot(importantReactions, range, biomass_max, biomass_min, target_max, target_min, model, directory)

for k = 1:length(importantReactions)
    
    reaction_index = (importantReactions(k));

    h = figure('visible', 'off');
    
    % target variability when making biomass:
    shadedplot(range, biomass_max(k,:), biomass_min(k,:), 'r', 'r');
    
    hold('on');
    
    % biomass variability when making target:
    shadedplot(range, target_max(k,:), target_min(k,:), 'b');
    
    alpha(0.1)
    title(model.rxnNames(reaction_index));
    ylabel('flux');
    xlim([1 100]);
    set(gca, 'Xtick', [1,25,50,75,100], 'XTickLabel', {'100% target', '', '', '', '100% bm'});
    
    
    s = char(strcat(directory, '/', model.rxns(reaction_index), '.png'));
    saveas(h, s);
    close(h);
end