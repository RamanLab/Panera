function panModel = createPanGenusModels(spNames,models,panPath,genusName,WithoutFutileCycle)
% The function constructs the Pan-Genus Metabolic Model (PGMM) from the
% existing GSMM for a genus
%
%   INPUTS:
%       spNames    : Information of the species needs to be included
%       models      : A struct containing all the species models
%       panPath     : Path to save the PanModel
%       genusName   : Name of the genus
%       WithoutFutileCycle : Binary input to decide whether thefutile cycle
%                            needs to be removed
%
%   OPTIONAL INPUTS:
%       WithoutFutileCycle : Binary value indicating whether to remove the futile reactions (default: )
%
%   OUTPUT:
%       panModel    : The reconstructed panModel
%
%   Author: Indumathi Palanikumar, 2023

% Define default input parameters if not specified
if ~exist('WithoutFutileCycle','var')
    WithoutFutileCycle = 1;
end

% Builds a pan-model from all models corresponding to strains that belong
% to the respective taxon.
tol = 1e-5;
% List Western diet constraints
dietConstraints = readtable('WesternDietAGORA2.txt');
dietConstraints = table2cell(dietConstraints);
dietConstraints(:, 2) = cellstr(num2str(cell2mat(dietConstraints(:, 2))));

% Databases
metaboliteDatabase = readtable('MetaboliteDatabase.txt', 'Delimiter', '\t', 'ReadVariableNames', false,'format','%s%s%s%s%s%s%s%s%s%s%s%s%s');
metaboliteDatabase = table2cell(metaboliteDatabase);
bioMetaboliteDatabase = table2cell(readtable('biomassMetaboliteDatabase.txt', 'Delimiter', '\t', 'ReadVariableNames', false,'format','%s%s%s%s%s%s%s%s%s%s%s%s%s'));
database.metabolites = [metaboliteDatabase; bioMetaboliteDatabase];
reactionDatabase = readtable('ReactionDatabase.txt', 'Delimiter', '\t', 'ReadVariableNames', false);
bioReactionDatabase = table2cell(readtable('biomassReactionDatabase.txt','Delimiter','\t','ReadVariableNames', false));
reactionDatabase = table2cell(reactionDatabase);
database.reactions = [reactionDatabase; bioReactionDatabase];


if size(models, 1) == 1
    % Loading model
    model = models{1,1};
    % rename biomass reaction to be comparable with other pan-models
    bio = find(contains(model.rxns, [spNames{1}, '_bio']));
    model.rxns{bio, 1} = 'biomassPan';
    orgRxns = model.rxns;
elseif size(models, 1) > 1
    for k = 1:size(models, 1)
        model = models{k,1};

        % Change the biomass equation and mets into orgSp
        bioMetInd = find(contains(model.mets,'biomass'));
        bioMets(k,1) = model.mets(bioMetInd);

        % Analyzing the reactions in the organisms
        bioRxnInd(k) = find(model.c);
        bioRxns(k,1) = model.rxns(find(model.c));
        orgRxns{k,1} = model.rxns;
        bio = find(strncmp(model.rxns, 'bio', 3));
        if k == 1
            panModel.rxns = model.rxns;
            panModel.grRules = model.grRules;
            panModel.rxnNames = model.rxnNames;
            panModel.subSystems = model.subSystems;
            panModel.lb = model.lb;
            panModel.ub = model.ub;
            forms = printRxnFormula(model, model.rxns, false, false, false, [], false);
            panModel.formulas = forms;
            % biomass products and substrates with coefficients
            bioPro = model.mets(find(model.S(:, bio) > 0), 1);
            bioProSC = full(model.S(find(model.S(:, bio) > 0), bio));
            bioSub = model.mets(find(model.S(:, bio) < 0), 1);
            bioSubSC = full(model.S(find(model.S(:, bio) < 0), bio));
        else
            panModel.rxns = [panModel.rxns; model.rxns];
            panModel.grRules = [panModel.grRules; model.grRules];
            panModel.rxnNames = [panModel.rxnNames; model.rxnNames];
            panModel.subSystems = [panModel.subSystems; model.subSystems];
            panModel.lb = [panModel.lb; model.lb];
            panModel.ub = [panModel.ub; model.ub];
            forms = printRxnFormula(model, model.rxns, false, false, false, [], false);
            panModel.formulas = [panModel.formulas; forms];
            % biomass products and substrates with coefficients
            bioPro = [bioPro; model.mets(find(model.S(:, bio) > 0), 1)];
            bioProSC = [bioProSC; full(model.S(find(model.S(:, bio) > 0), bio))];
            bioSub = [bioSub; model.mets(find(model.S(:, bio) < 0), 1)];
            bioSubSC = [bioSubSC; full(model.S(find(model.S(:, bio) < 0), bio))];
        end
    end
    % take out biomass reactions
    bio = [find(contains(panModel.rxns, 'biomass'));find(contains(panModel.rxns, 'bio1'))];
    bioEx = find(contains(panModel.rxns, 'EX_biomass'));
    bioRxnFormulas = panModel.formulas(setdiff(bio,bioEx));
    panModel.rxns(bio) = [];    panModel.grRules(bio) = [];
    panModel.rxnNames(bio) = [];    panModel.subSystems(bio) = [];
    panModel.lb(bio) = [];    panModel.ub(bio) = [];
    panModel.formulas(bio) = [];
    % set up data matrix for rBioNet
    [uniqueRxns, oldInd] = unique(panModel.rxns);
    rbio.data = cell(size(uniqueRxns, 1), 14);
    rbio.data(:, 1) = num2cell(ones(size(rbio.data, 1), 1));
    rbio.data(:, 2) = uniqueRxns;
    rbio.data(:, 3) = panModel.rxnNames(oldInd);
    rbio.data(:, 4) = panModel.formulas(oldInd);
    rbio.data(:, 6) = panModel.grRules(oldInd);
    rbio.data(:, 7) = num2cell(panModel.lb(oldInd));
    rbio.data(:, 8) = num2cell(panModel.ub(oldInd));
    rbio.data(:, 10) = panModel.subSystems(oldInd);
    rbio.description = cell(7, 1);
    % build model with rBioNet
    model = data2model(rbio.data, rbio.description, database);
    % build biomass reaction through a linear combination of species
    % biomass
    bioForm = '';
    for x = 1: length(bioRxns)
        model = addReaction(model,bioRxns{x,1},'reactionFormula',bioRxnFormulas{x});
    end

    % add pan Biomasss as a ratio of other biomass reactions
    % coefs = [0.25,0.5,0.25];
    % add biomass reaction to pan model
    for x = 1: length(bioRxns)
        bioForm = [bioForm, num2str(1/k), ' ', bioMets{x} , ' + '];
    end
    bioForm = bioForm(1:end - 3);
    bioForm = [bioForm, ' -> '];
    bioForm = [bioForm, num2str(1), ' ', 'biomass[c]', ' + '];
    bioForm = bioForm(1:end - 3);
    model = addReaction(model, 'biomassPan','reactionFormula',bioForm);
    model = addReaction(model, 'EX_biomass(e)','reactionFormula',' biomass[c] -> ');
    model.comments{end + 1, 1} = '';
    model.citations{end + 1, 1} = '';
    model.rxnConfidenceScores{end + 1, 1} = '';
    model.rxnECNumbers{end + 1, 1} = '';
    model.rxnKEGGID{end + 1, 1} = '';
end
% update some fields to new standards
model.osenseStr = 'max';
if isfield(model, 'rxnConfidenceScores')
    model = rmfield(model, 'rxnConfidenceScores');
end
model.rxnConfidenceScores = zeros(length(model.rxns), 1);
for k = 1:length(model.rxns)
    model.subSystems{k, 1} = cellstr(model.subSystems{k, 1});
    model.rxnKEGGID{k, 1} = '';
    model.rxnECNumbers{k, 1} = '';
end

for k = 1:length(model.mets)
    if strcmp(model.metPubChemID{k, 1}, '[]') || isempty(model.metPubChemID{k, 1})
        model.metPubChemID{k, 1} = string;
    end
    if strcmp(model.metChEBIID{k, 1}, '[]') || isempty(model.metChEBIID{k, 1})
        model.metChEBIID{k, 1} = string;
    end
    if strcmp(model.metKEGGID{k, 1}, '[]') || isempty(model.metKEGGID{k, 1})
        model.metKEGGID{k, 1} = string;
    end
    if strcmp(model.metInChIString{k, 1}, '[]') || isempty(model.metInChIString{k, 1})
        model.metInChIString{k, 1} = string;
    end
    if strcmp(model.metHMDBID{k, 1}, '[]') || isempty(model.metHMDBID{k, 1})
        model.metHMDBID{k, 1} = string;
    end
end
model.metPubChemID = cellstr(model.metPubChemID);
model.metChEBIID = cellstr(model.metChEBIID);
model.metKEGGID = cellstr(model.metKEGGID);
model.metInChIString = cellstr(model.metInChIString);
model.metHMDBID = cellstr(model.metHMDBID);
% fill in descriptions
model = rmfield(model, 'description');
model.description.organism = ['pan',genusName];
model.description.name = ['pan',genusName];
model.description.author = 'https://vmh.life';
model.description.date = date;


% Rebuild model consistently
%model = rebuildModel(model,database);
model=changeObjective(model,'biomassPan');

% remove duplicate reactions
% Will remove reversible reactions of which an irreversible version is also
% there but keep the irreversible version.
[modelRD, removedRxnInd, keptRxnInd] = checkDuplicateRxn(model);
removedRxns = model.rxns(removedRxnInd);
keptRxns = model.rxns(keptRxnInd);
% test if the model can still grow
modelRD = useDiet(modelRD,dietConstraints);
FBA=optimizeCbModel(modelRD,'max');
if FBA.f > tol
    model=model;
else
    toRM={};
    modelTest=model;
    for k=1:length(removedRxnInd)
        modelTest=removeRxns(modelTest,modelTest.rxns(removedRxnInd(k)));
        FBA=optimizeCbModel(modelTest,'max');
        if FBA.f > tol
            toRM{k} =  modelTest.rxns{removedRxnInd(k)};
            model
        else
            toRM{k} =  modelTest.rxns{keptRxnInd(k)};
        end
        modelTest=removeRxns(modelTest,toRM{k});
    end
    model=removeRxns(model,toRM);
end

% If there is only one model, it will be kept as such except for the
% biomass equation
if size(models, 1) == 1
    panModel = model
    PanRxnPresence = ones(length(panModel.rxns),1);
else
    modelPath = [panPath filesep 'pan' genusName '.mat'];
    if WithoutFutileCycle == 1
        panModel = modelWithoutFutileCycle(model,modelPath,dietConstraints,database);
    else
        panModel = model;
    end
    %PanRxnPresence construction
    [RxnPresenceMatrix,PanRxnPresence] = createPanRxnPresence(orgRxns,length(spNames),removedRxns,keptRxns,panModel,1);
end

panModel.rxnPresenceMat = PanRxnPresence;
% Adding species list to the panModels
panModel.spList = spNames;

savePath = [panPath filesep 'pan' genusName '.mat'];
save(savePath, 'panModel');

end
