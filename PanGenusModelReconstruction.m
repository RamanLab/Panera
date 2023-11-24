function panModel = PanGenusModelReconstruction(modelPath,panPath,genusInfo,infoFilePath,dietApplied,dietFilePath)
% The function prepares the input to construct the Pan-Genus Metabolic Model (PGMM) from the
% existing GSMM for a genus
%
%   INPUTS:
%       modelPath  : Path to GSMMs to include in PGMM
%       panPath     : Path tp save the PanModel
%       genusInfo   : Information about the genus that needs PGMM
%
%   OPTIONAL INPUTS:
%       infoFilePath : Path to the information file containing taxonomic information
%       dietApplied : Binary value indicating whether to incorporate diet
%       information to the model (default: )
%       dietFilePath : If the diet needs to be applied, path to the diet
%       file
%
%   OUTPUT:
%       panModel    : A struct of the produced panModels
%
%   Author: Indumathi Palanikumar, 2023


% Check for the genus name and file path for the species information
if ~exist(infoFilePath)
    AGORAinfoFile = '/home/indumathi/Desktop/PhD/Research work/Resources/AGORA_infoFile.xlsx';
    InfoFile = readtable(AGORAinfoFile);
else
    InfoFile = readtable(infoFilePath);
end
if ~isempty(genusInfo)
    genusNames = genusInfo;
else
    genusNames = unique(InfoFile.Genus);
end

if ~(isempty(dir(panPath)))
    dInfo = dir(panPath);
    modelList={dInfo.name};
    modelList=modelList';
    modelList=strrep(modelList,{'.mat'},'');
    modelList=strrep(modelList,{'pan'},'');
    genusNames=setdiff(genusNames,modelList);
end
genusNames = intersect(genusNames,InfoFile.Genus);

for i = 1: length(genusNames)
    genusName = genusNames{i};

    % Information about the species in the genus
    SpeciesInfo = InfoFile(strcmp(InfoFile.Genus,genusName),:);

    % Loading the models from the genus
    orgTable = table2cell(SpeciesInfo(:,1));
    nameTags = regexprep(orgTable,'_',' ');
    orgModels = {};

    % Models need to be changed
    for i =1: length(orgTable)
        modelName = [modelPath,orgTable{i},'.mat'];
        orgModels{i,1} = readCbModel(modelName);
        % Change Reaction names to species specific
        bioRxnInd(i,1) = find(orgModels{i,1}.c);
        bioRxn(i,1) = orgModels{i,1}.rxns(bioRxnInd(i));
        orgModels{i,1}.rxns(bioRxnInd(i,1)) = {[orgTable{i},'_',bioRxn{i,1}]};
        % Change metabolite names to species specific
        bioMetInd(i,1) = find(contains(orgModels{i,1}.mets,'biomass'));
        bioMet(i,1) = orgModels{i,1}.mets(bioMetInd(i,1));
        orgModels{i,1}.mets(bioMetInd(i,1)) = {[orgTable{i},'_',bioMet{i,1}]};
        %orgRxns{i,1} = orgModels{i}.rxns;
    end

    if dietApplied == 1
        if ~isempty(dietFilePath)
            Diet = adaptVMHDietToAGORA(dietFilePath,'AGORA');
        else
            Diet = adaptVMHDietToAGORA('/home/indumathi/Desktop/PhD/Research work/Resources/Diet files/EUdiet','AGORA');
        end
    end

    %% Growth analysis of individual species on specific diet condition
    % default - EU diet
    for i =1: length(orgTable)
        modelSp = orgModels{i,1};

        if dietApplied == 1
            % Analyze the growth of the organisms
            modelSp = useDiet(modelSp,Diet);
        end
        sol = optimizeCbModel(modelSp);
        indFlux(i,1) = sol.f;
    end

    %% Creation of Pan genus model
    [panModel] = createPanGenusModels(orgTable,orgModels,panPath,genusName,1);
end
end