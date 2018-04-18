function consModel = createModel(cobra_model, chemicalInd, biomassInd, organismObjectiveInd, constraintsList)

    % Added by natalie stanford so that a COBRA model can be loaded.
    % Has a similar format to the original.
    consModel = cobra_model; 
    [metNum, rxnNum] = size(consModel.S);
    consModel.row_lb = zeros(metNum, 1);
    consModel.row_ub = zeros(metNum, 1);

    % Set target chemical:
    consModel.C_chemical = zeros(rxnNum, 1);
    consModel.C_chemical(chemicalInd) = 1;
    consModel.chemicalInd = chemicalInd;

    % Set minimum biomass / growth:
    consModel.biomassInd = biomassInd;
    consModel.lb(biomassInd) = 0.1;

    % Set cellular objective (not necessarily growth, but may be):
    consModel.organismObjective = zeros(rxnNum, 1);
    consModel.organismObjective(organismObjectiveInd) = 1;
    consModel.organismObjectiveInd = organismObjectiveInd;

    % Remove small reactions:
    consModel.S(consModel.S > -(10^-3) & consModel.S < 0) = 0;
    consModel.S(consModel.S < (10^-3) & consModel.S > 0) = 0;

    % Is this necessary? The bounds should already be in the model.
    % This has the effect of fixing glucose uptake to 10, and ATP
    % maintenance to 8.39, as opposed to specifying these values as their
    % lower bounds...
    for k = 1:length(constraintsList.reactions)
        rctnNo = find(ismember(consModel.rxns, constraintsList.reactions(k)));
        if isempty(rctnNo)
            error('No matching reactions')
        elseif length(rctnNo) > 1
            error('Too many matching reactions')
        else
            consModel.lb(rctnNo) = constraintsList.lb(k);
            consModel.ub(rctnNo) = constraintsList.ub(k);
        end
    end