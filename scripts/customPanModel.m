function modifiedPanModel = customPanModel(panModel,speciesProb)
% Function to customize the reconstructed PGMM with the provided species
% probability vector
%
%   INPUT:
%       panModel    :   Reconstructed PGMM
%       speciesProb  :  Species probability vector
%
%   OUTPUT:
%       modifiedPanModel    : Customized PGMM
%
%   Author: Indumathi Palanikumar, 2023

% Preprocesing species probability vector
RxnProbabilities = panModel.rxnPresenceMat*speciesProb;
rxnIndex = find(RxnProbabilities > 0);
rxnPresence = zeros(size(panModel.rxns,1),1);
rxnPresence(rxnIndex) = 1; rxnPresence(~rxnIndex) = 0;
PanRxnAbsence = find(rxnPresence== 0);
panModelAbun = panModel;

% Updating biomass reaction according to the species probability
bioInd = findRxnIDs(panModel,'biomassPan');
bioSub = panModel.mets(find(panModel.S(:, bioInd) < 0), 1);
bioPro = panModel.mets(find(panModel.S(:, bioInd) > 0), 1);
bioForm=[];
for x = 1: length(bioSub)
    bioForm = [bioForm, num2str(speciesProb(x)), ' ', bioSub{x} , ' + '];
end
bioForm = bioForm(1:end - 3);
bioForm = [bioForm, ' -> '];
bioForm = [bioForm, num2str(1), ' ',bioPro{1}, ' + '];
bioForm = bioForm(1:end - 3);
panModelAbun = addReaction(panModelAbun, 'biomassPan','reactionFormula',bioForm);

% Selecting the reactions in Pan model
panModelAbun = changeRxnBounds(panModelAbun,panModelAbun.rxns(PanRxnAbsence),0,'b');
panModelAbun = changeObjective(panModelAbun,'biomassPan');

% Checking whether the model works
panModelSol = optimizeCbModel(panModelAbun);
if panModelSol.f > 1e-5
    modifiedPanModel = panModelAbun;
else
    disp('The model cannot grow with this abundance ratio')
    modifiedPanModel = {};
end
end
