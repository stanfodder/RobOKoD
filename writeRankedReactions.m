function writeRankedReactions(cobra_model, rankedReactions, filename, ko)
%PRINTRANKEDREACTIONS Summary of this function goes here
%   Detailed explanation goes here
    fileID = fopen(filename, 'w');
    
    if ko
        % Knockout:
        fprintf(fileID, '%s\t%s\t%s\t%s\t%s\n', 'Reaction id', 'Reaction name', 'Gene association', 'High flux loss?', 'KO ranking');
    else
        % Overexpression:
        fprintf(fileID, '%s\t%s\t%s\t%s\n', 'Reaction id', 'Reaction name', 'Gene association', 'Overexpression ranking');
    end
    
    for row = 1:length(rankedReactions)
        reaction_index = rankedReactions(row, 1);
        fprintf(fileID, '%s\t%s\t%s\t', char(cobra_model.rxns(reaction_index)), char(cobra_model.rxnNames(reaction_index)), char(cobra_model.grRules(reaction_index)));;
        
        if ko
            fprintf(fileID, '%d\t%f\n', rankedReactions(row, 2), rankedReactions(row, 3));
        else
            fprintf(fileID, '%f\n', rankedReactions(row, 2));
        end
    end
    
    fclose(fileID);
end