clc;
clear;

% General parameters"
model_filename = 'iNS125.xml';
biomass_reaction_id = 'Biomass_Ecoli_core_w_GAM';
target_reaction_id = 'er_027';
max_knockouts = 5;
results_directory = 'iNS125_butanol_results';
geneKO = 0;
MCT = 0;

% Constraints for OptKnock and RobustKnock:
constraint_reaction_list = {'EX_glc(e)', 'ATPM'}; 
constraint_reaction_values = [-10, 8.39];


% Run RobOKoD:
fprintf('%s\n', 'RobOKoD');
fprintf('%s\n', '-------');

robokod(model_filename, biomass_reaction_id, target_reaction_id, max_knockouts, results_directory,geneKO,MCT);


% Run OptKnock:
fprintf('\n');
fprintf('%s\n', 'OptKnock');
fprintf('%s\n', '--------');

optKnockSolution = runOptKnock(model_filename, biomass_reaction_id, target_reaction_id, ...
    constraint_reaction_list, constraint_reaction_values, max_knockouts, 0);

fprintf('%s\n', 'OptKnock predicts the following knockouts:');

for i = 1:length(optKnockSolution.rxnList) 
    fprintf('%s\n', optKnockSolution.rxnList{i}); 
end


% Run OptKnock (excluding exchange reactions):
fprintf('\n');
fprintf('%s\n', 'OptKnock (exchange reactions excluded)');
fprintf('%s\n', '--------------------------------------');

optKnockSolution = runOptKnock(model_filename, biomass_reaction_id, target_reaction_id, ...
    constraint_reaction_list, constraint_reaction_values, max_knockouts, 1);

fprintf('%s\n', 'OptKnock predicts the following knockouts:');

for i = 1:length(optKnockSolution.rxnList) 
    fprintf('%s\n', optKnockSolution.rxnList{i}); 
end


% Run RobustKnock:
fprintf('\n');
fprintf('%s\n', 'RobustKnock');
fprintf('%s\n', '-----------');

robustKnockSolutionRxnList = runRobustKnock(model_filename, biomass_reaction_id, target_reaction_id, ...
    constraint_reaction_list, constraint_reaction_values, max_knockouts);

fprintf('%s\n', 'RobustKnock predicts the following knockouts:');

for i = 1:length(robustKnockSolutionRxnList) 
    fprintf('%s\n', robustKnockSolutionRxnList{i}); 
end