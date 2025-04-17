function [translatedAbundances,normalizedAbundances,unmappedRows]=modifiedtranslateMetagenome2AGORA(MetagenomeAbundancePath,sequencingDepth)
% Translates organism identifiers in a published metagenomic or 16S rRNA
% data file with organism abundances (retrieved e.g., from  MetaPhlAn) to
% AGORA pan-model IDs. This will not catch every case since the format of
% input files with abundance data greatly varies. Feel free to modify this
% function and submit a pull request to enable more input files to be
% translated to AGORA. Moreover, slight spelling variations in taxa across
% input files may  result to taxa not being mapped. Check the unmappedRows
% output to identify these cases and modify the function accordingly.
% Pan-models that can be used to create microbiome models in mgPipe can be
% created with the function createPanModels.
%
% USAGE:
%
%   [translatedAbundances,normalizedAbundances,unmappedRows]=translateMetagenome2AGORA(MetagenomeAbundancePath,sequencingDepth)
%
% INPUT:
%   MetagenomeAbundancePath   String containing the path to csv file with
%                             organism abundance data retrieved from
%                             16S rRNA or metagenomic samples (example:
%                             'SRP065497_taxonomy_abundances_v3.0.tsv').
%
% OPTIONAL INPUT:
%   sequencingDepth           Sequencing depth on the taxonomical level
%                             in the input data (e.g., genus, species).
%                             Allowed inputs are 'Species','Genus',
%                             'Family','Order', 'Class', 'Phylum'. 
%                             Default: 'Species'.
%
% OUTPUTS:
%   translatedAbundances      Abundances with organism names from the
%                             input file translated to AGORA pan-model IDs
%   normalizedAbundances      Translated abundances normalized so they sum
%                             up to 1 for each sample
%   unmappedRows              Taxa on the selected taxonomical level that
%                             could not be mapped to AGORA pan-models
%
% .. Author: Almut Heinken, 02/2019
%                           10/2019: edited mapping based on an output file
%                           of the NG-Tax pipeline (Ramiro-Garcia et al.,
%                           F1000Res 2016). Made input more intuitive.

if nargin <2
    sequencingDepth='Species';
end
% read the csv file with the abundance data
metagenome_abundance = readtable(MetagenomeAbundancePath, 'ReadVariableNames', false ,'FileType','text');
metagenome_abundance = table2cell(metagenome_abundance);
%metagenome_abundance{1,1}='';

% list the taxonomical levels
taxLevels={
    'Phylum'
    'Class'
    'Order'
    'Family'
    'Genus'
    'Species'
    };

%[~, infoFile, ~] = xlsread('AGORA_infoFile.xlsx');  % create the pan-models
[~, infoFile, ~] = xlsread('AGORA2_infoFile.xlsx');  % create the pan-models

% format the input file so it can be matched
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),';','|');

% match different phrasing
% This may not catch all cases! Feel free to add additional cases of
% mismatching nomenclature and submit a pull request.
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Candidatus_','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_Family_XI_Incertae_Sedis',' Incertae Sedis XI');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_Family_XIII_Incertae_Sedis',' Incertae Sedis XIII');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'typhimurium','enterica');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_1','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_2','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_3','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_4','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_5','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_6','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_7','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_8','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_9','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),' family','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),' order','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Ruminococcus_gauvreauii_group','Ruminococcus');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Ruminococcus_gnavus_group','Blautia');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Ruminococcus_torques_group','Blautia');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Eubacterium_coprostanoligenes_group','Eubacterium');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Eubacterium_eligens_group','Eubacterium');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Eubacterium_hallii_group','Anaerobutyricum');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Eubacterium_ruminantium_group','Eubacterium');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Eubacterium_ventriosum_group','Eubacterium');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Eubacterium_xylanophilum_group','Eubacterium');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'Clostridium_innocuum_group','Erysipelatoclostridium');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_sensu_stricto','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_sensu_stricto','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'[','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),']','');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),' ','_');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'.','_');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'-','_');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'__','_');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'___','_');

% Remove unclassified organisms
metagenome_abundance(~cellfun(@isempty, strfind(metagenome_abundance(:,1),'_uncl')),:) = [];

% replace unknown species with pan-genus model
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),'_cf','');

% summarize duplicate entries
[uniqueA,i,j] = unique(metagenome_abundance(:,1));
n  = accumarray(j(:),1);
Dupes=uniqueA(find(n>1));
delArray=[];
cnt=1;
for i=1:length(Dupes)
    indexToDupes = find(strcmp(metagenome_abundance(:,1),Dupes{i}));
    for j=1:length(indexToDupes)
        for k=2:size(metagenome_abundance,2)
            metagenome_abundance{indexToDupes(1),k}=num2str(str2double(metagenome_abundance{indexToDupes(1),k})+str2double(metagenome_abundance{indexToDupes(j),k}));
        end
        delArray(cnt,1)=indexToDupes(j);
        cnt=cnt+1;
    end
end
metagenome_abundance(delArray,:)=[];

% match and find overlapping IDs
%metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),mpsDepth,'');
metagenome_abundance(:,1)=strrep(metagenome_abundance(:,1),' ','_');

% Identify only convert the genus specific data
if sequencingDepth == 'Species'
    taxaNames = metagenome_abundance(:,1);
    testTaxa = cellfun(@(x) strsplit(x, '_'), taxaNames, 'UniformOutput', false);
    index4Pan = find(cellfun(@(x) length(x) == 1, testTaxa));
    index4sp = find(cellfun(@(x) length(x) == 2, testTaxa));
    panTaxa = unique(infoFile(2:end, find(strcmp(infoFile(1, :), 'Genus'))));
    panTaxa(:,1)=strcat('pan',panTaxa(:,1));
end
metagenome_abundance(index4Pan,1)=strcat('pan',metagenome_abundance(index4Pan,1));

% remove and report taxa in abundance data table that could not be found in
% AGORA and/or are not of the desired sequencing depth

findTaxCol = find(strcmp(infoFile(1, :), sequencingDepth));
allTaxa = unique(infoFile(2:end, findTaxCol));
allTaxa(:,1)=strrep(allTaxa(:,1),' ','_');
%allTaxa(:,1)=strcat('pan',allTaxa(:,1));
[unmappedRows,IA] = setdiff(metagenome_abundance(:,1),[allTaxa; panTaxa]);
if ~isempty(IA)
    IA(1,:)=[];
    metagenome_abundance(IA,:)=[];
end

if ~isempty(index4sp)
    % Adding AGORA names to the abundances file
    spNames = metagenome_abundance(index4sp,1);
    AGORAID =  cellfun(@(x) strmatch(x, strrep(infoFile(:,5),' ','_')), spNames, 'UniformOutput', false);
    AGORANames = cellfun(@(x) infoFile(x(1),1), AGORAID);
    metagenome_abundance(index4sp,1) = AGORANames;
end

translatedAbundances=metagenome_abundance;

% normalize the abundances so that sum for each individual is 1
normalizedAbundances=translatedAbundances;
for i=2:size(translatedAbundances,2)
    if sum(cell2mat(translatedAbundances(2:end,i)))>0
        for j=2:size(translatedAbundances,1)
            normalizedAbundances{j,i}=num2str(cell2mat(translatedAbundances(j,i))/sum(cell2mat(translatedAbundances(2:end,i))));
        end
    else
        for j=2:size(translatedAbundances,1)
            normalizedAbundances{j,i}=translatedAbundances{j,i};
        end
    end
end

% export the translated abundances
writetable(cell2table(translatedAbundances),'translatedAbundances','FileType','text','WriteVariableNames',false,'Delimiter','tab');
writetable(cell2table(normalizedAbundances),'normalizedAbundances','FileType','text','WriteVariableNames',false,'Delimiter','tab');

end