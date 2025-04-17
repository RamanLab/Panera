% The script provides tutorial on creating PGMM to use in community
% modelling using MMT in CobraToolbox. The steps including abundance file
% processing to comparison of microbial community output

%% Step -1: Processing abundance files 
% Initializing Cobratoolbox after installation and addition to path
initCobraToolbox(false)

%Setting a solver
solverOK=changeCobraSolver('ibm_cplex','LP');

% Path to AGORA2 Mat files
% modPath = '/path/to/AGORA2/mat files/';

% Path to the downloaded PaneraFolder
% folderPath = '/path/to/folder/';

% Adding path
addpath(genpath(folderPath))

% Translate abundances to AGORA names
[translatedAbundances,normalizedAbundances,unmappedRows] = modifiedtranslateMetagenome2AGORA([folderPath '/tutorial/data/sampleAbundanceFile.csv'],'Species');

%% Step -2: Create Genus specific models

% Identify the required genus models
genusNames = normalizedAbundances(strmatch('pan',normalizedAbundances(:,1)),1);
genusNames = strrep(genusNames,'pan','');

agoraPath = modPath;
panPath=[pwd filesep '/PanGenusModels_AGORA2'];
infoFileName = [folderPath 'AGORA2_infoFile.xlsx'];

if (~exist(dir(panPath)))
    mkdir(panPath)
end

% Reconstructing panModel - panModel created will be assigned equal
% probabilities to the species present in a genus
PanGenusModelReconstruction(agoraPath,panPath,genusNames,infoFileName,0,[]);

%% Step-3: Using PGMM in microbial community modelling

dietFilePath = [folderPath 'tutorial/data/EUdiet'];
analysisPath = [folderPath 'Analysis/EU/'];
resPath = analysisPath;
metaboliteInfo=table2cell(readtable([folderPath 'dependencies/MetaboliteInformation.csv'], 'ReadVariableNames', false));
computeProfiles = true;
saveConstrModels = true;
normAbunFilePath = [folderPath 'tutorial/data/normalizedAbundances_Genus.txt'];
numWorkers = 12;
[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(agoraPath, normAbunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers);

%% Optional Step : Modelling communities using only GSMM for comparison

% Translate abundances to AGORA names
[translatedAbundances,normalizedAbundances,unmappedRows] = modifiedtranslateMetagenome2AGORA([folderPath '/tutorial/data/sampleAbundanceFile.csv'],'Species');

% Code to simulate the GSMM-based communities. Uncomment the following lines and run to simulte
% analysisPath = [folderPath 'Analysis/GSMM/'];
% resPath = analysisPath;
% normAbunFilePath = [folderPath 'tutorial/data/normalizedAbundances.txt'];
% [init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(agoraPath, normAbunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers);

% Loading net secretion fluxes
pgmm_flux = table2cell(readtable([mainFolderPath 'Analysis/EU/inputDiet_net_secretion_fluxes.csv']));
gsmm_flux = table2cell(readtable([mainFolderPath 'Analysis/GSMM/inputDiet_net_secretion_fluxes.csv']));

flux = zeros(796,20);
flux(:,1:10) = cell2mat(pgmm_flux(:,2:11));
[~,ind] = ismember(gsmm_flux(:,1),pgmm_flux(:,1));
flux(ind, 11:20) = cell2mat(gsmm_flux(:,2:11));

% PCoA plotting
JD = pdist(flux','euclidean'); % Similarities are not accepted; Needs dissimilarity metrix
[Y, eigvals] = cmdscale(JD);
P = [eigvals eigvals / max(abs(eigvals))];
expr = [eigvals/sum(eigvals)];
sampNames = [arrayfun(@(x) sprintf('P%d',x), 1:10, 'UniformOutput', false), arrayfun(@(x) sprintf('G%d',x), 1:10, 'UniformOutput', false)]';

figure
plot(Y(:, 1), Y(:, 2), 'bx')
xlabel(strcat('PCoA1: ',num2str(round(expr(1)*100,2)),'% of explained variance'));
ylabel(strcat('PCoA2: ',num2str(round(expr(2)*100,2)),'% of explained variance'));
text(Y(:,1),Y(:,2),sampNames,'HorizontalAlignment','left');%to insert numbers
title('PCoA of net secretion profiles');

% Quantitative and qualitative variation
flux_bi = flux > 0;
JD = pdist(flux_bi','euclidean');

