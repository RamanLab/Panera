function model = modelWithoutFutileCycle(model,modelPath,dietConstraints,database)
% Function to remove the curated futile reactions from the reconstructed
% model
%
%   INPUTS:
%       model       : Reconstructed model
%       modelPath   : Path to store the refined model
%       dietConstraints : Table with the diet constraints 
%       database    : Metabolite and reaction database to add or remove the
%       reactions
%
%   OUTPUT:
%       model   : Modified model
%
%   Author: Indumathi Palanikumar, 2023

% Create table with information on reactions to replace to remove futile
% cycles.
tol = 1e-5;

% Loading futile reactions to remove and replace
reactionsToReplace = table2cell(readtable('futileReactionCuration.csv','ReadVariableNames',false));

model = useDiet(model, dietConstraints);
model = changeObjective(model, 'DM_atp_c_');
FBA = optimizeCbModel(model, 'max');
% Ensure that pan-models can still produce biomass
model = changeObjective(model, 'biomassPan');
if FBA.f > 50
    for j = 2:size(reactionsToReplace, 1)
        rxns = strsplit(reactionsToReplace{j, 1}, ' AND ');
        go = true;
        for k = 1:size(rxns, 2)
            RxForm = database.reactions{find(ismember(database.reactions(:, 1), rxns{k})), 3};
            if contains(RxForm,'[e]') && any(contains(model.mets,'[p]'))
                newName=[rxns{k} 'pp'];
                % make sure we get the correct reaction
                newForm=strrep(RxForm,'[e]','[p]');
                rxnInd=find(ismember(database.reactions(:, 1), {newName}));
                if ~isempty(rxnInd)
                    dbForm=database.reactions{rxnInd, 3};
                    if checkFormulae(newForm, dbForm) && any(contains(model.mets,'[p]'))
                        rxns{k}=newName;
                    end
                end
            end
            if isempty(find(ismember(model.rxns, rxns{k})))
                go = false;
            end
        end
        if go
            % account for periplasmatic versions
            replacePP=0;
            RxForm = database.reactions{find(ismember(database.reactions(:, 1), reactionsToReplace{j, 2})), 3};
            if contains(RxForm,'[e]') && any(contains(model.mets,'[p]'))
                newName=[reactionsToReplace{j, 2} 'pp'];
                % make sure we get the correct reaction
                newForm=strrep(RxForm,'[e]','[p]');
                dbForm=database.reactions{find(ismember(database.reactions(:, 1), {newName})), 3};
                replacePP=1;
            end
            % Only make the change if biomass can still be produced
            if replacePP
                modelTest = removeRxns(model, newName);
            else
                modelTest = removeRxns(model, reactionsToReplace{j, 2});
            end
            if ~isempty(reactionsToReplace{j, 3})
                RxForm = database.reactions{find(ismember(database.reactions(:, 1), reactionsToReplace{j, 3})), 3};
                if replacePP
                    % create a new formula
                    RxForm = database.reactions{find(ismember(database.reactions(:, 1), reactionsToReplace{j, 3})), 3};
                    if contains(RxForm,'[e]') && any(contains(model.mets,'[p]'))
                        newName=[reactionsToReplace{j, 3} 'ipp'];
                        % make sure we get the correct reaction
                        newForm=strrep(RxForm,'[e]','[p]');
                        rxnInd=find(ismember(database.reactions(:, 1), {newName}));
                        if ~isempty(rxnInd)
                            dbForm=database.reactions{rxnInd, 3};
                            if checkFormulae(newForm, dbForm) && any(contains(model.mets,'[p]'))
                                RxForm=dbForm;
                            end
                        else
                            % if not present already, add to database
                            RxForm=newForm;
                            database.reactions(size(database.reactions,1)+1,:)={reactionsToReplace{j, 3},newName,RxForm,'0','','','','','','','',''};
                        end
                    end
                    modelTest = addReaction(modelTest, newName, RxForm);
                else
                    modelTest = addReaction(modelTest, reactionsToReplace{j, 3}, RxForm);
                end
            end
            FBA = optimizeCbModel(modelTest, 'max');
            if FBA.f > tol
                model = modelTest;
                % account for periplasmatic versions
                % fix some special cases
                if ~isempty(intersect(model.rxns,'CITt2ipp'))
                    model.rxns=strrep(model.rxns,'CITt2ipp','CITt2pp');
                end
                if ~isempty(intersect(model.rxns,'CITCAtiipp'))
                    model.rxns=strrep(model.rxns,'CITCAtiipp','CITCAtipp');
                end
            end
        end
    end
end
% set back to unlimited medium
model = changeRxnBounds(model, model.rxns(strmatch('EX_', model.rxns)), -1000, 'l');
% Rebuild model consistently
%         model = rebuildModel(model,database);
save(modelPath, 'model');
end
