% Tutorial for reconstructing Pan-Genus Metabolic Model (PGMM) with the provided species information
% The tutorial generates PGMM for the genus Escherichia 

% Define needed variables and paths
genusName = 'Escherichia';
infoFilePath = 'AGORA_infoFile.xlsx';

% Extracting species information for the genera from the information file
taxaInfo = readtable(infoFilePath);
SpeciesInfo = table2cell(taxaInfo(strcmp(taxaInfo.Genus,genusName),:));
nSp = length(SpeciesInfo(:,1));

% Path to the GSMM 
modelPath = 'AGORA reconstructions/mat files/'; % Path to the models
panPath = [pwd filesep 'PGMM_tutorial'];      % Path to store panModel
mkdir(panPath);

% Reconstructing PGMM
PanGenusModelReconstruction(modelPath,panPath,genusName,infoFilePath,[]);

% Loading the pan model
panModel = readCbModel([panPath filesep 'panEscherichia.mat']);

% Optimizing the panModel growth
solPanModel = optimizeCbModel(panModel);

% Customizing PGMM
% Creating a random vector of species probability
spProb = normalize([rand(nSp,1)],'norm',1);
customModel = customPanModel(panModel,spProb);

