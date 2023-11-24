function [RxnPresenceMatrix,PanRxnPresence] = createPanRxnPresence(Rxns,n,removedRxns,keptRxns,panModel,withoutFutile)
% Function to generate matrix of reaction presence/absence in the PGMM and
% to inform the models which reactions to include during customization
%
%   INPUTS:
%       Rxns    : All the reactions within a genus
%       n       :  Number of species incorporated in the model
%       removedRxns : Reactions removed from the model (either duplicated or futile reactions)
%       keptRxns    : Final set of reactions included in PGMM
%       panModel    : Reconstructed PGMM
%       withoutFutile  :  Binary value indicating whether to remove the futile reactions (default: 1)
%       
%   OUTPUT:
%       RxnPresenceMatrix   : Binary matrix representing the presence of
%       reactions from specific species
%       PanRxnPresence      : Binary matrix representing the presence of
%       species-specific reactions in PGMM
%
%   Author: Indumathi Palanikumar, 2023

if ~exist('WithoutFutile','var')
    WithoutFutile = 1;
end

% Change the reactions which has been deleted while creating a pan Model
for a = 1:n
    for i =1:length(removedRxns)
        indR = find(ismember(Rxns{a,1},removedRxns{i}));
        if isempty(indR)==0
            Rxns{a,1}{indR,1} = keptRxns{i};
        end
    end
end

% To make the reactions similar to Pan-genus model
panRxns = panModel.rxns;
PanModelRxns = [];RxnPresenceMatrix = [];
for j = 1:n
    rxns4pan = Rxns{j};
    % Pan reactions across the strains
    PanModelRxns = [PanModelRxns; rxns4pan];
    j = j+1;
end
PanRxnsUnique = unique(PanModelRxns);
for k=2:n+1
    RxnPresenceMatrix(:,k) = ismember(PanRxnsUnique,Rxns{k-1});
    k = k+1;
end
RxnPresenceMatrix = num2cell(RxnPresenceMatrix);
RxnPresenceMatrix(:,1) = PanRxnsUnique;

PanRxnPresence = zeros(1,size(RxnPresenceMatrix,2)-1);
for x = 1: length(panRxns)
        ind = find(ismember(RxnPresenceMatrix(:,1),panRxns(x)));
        if isempty(ind) == 0
            PanRxnPresence(x,:) = cell2mat(RxnPresenceMatrix(ind,2:end));
        else
            PanRxnPresence(x,:) = zeros(1,size(PanRxnPresence,2));
        end
end

if withoutFutile == 1
    replacedRxns ={'LYSt3','FDH','GLYO1i','EAR40x','PROt4','NO2t2','NHFRBO','L_LACD','PIt7ir','ABUTt2','ABTA','Kt3','CYTDt2','ASPte','ASPt2','FUMt','SUCCt','MALFADOi','MALFADOi','ALDD8x','PHPB2i','PPCK','PPCK','ACOAD1i','26DAPLLATi','DAPDAi','NAt3','MALt2','URIt2','HYXNti','URAt2','URAt2','CSNt2','XANt2','XPPT','XANt2','ARGt2','r1144','GLYt2','L_LACt2','G3PD8i','ACOAD1fi','D_GLY3PRi','NPR','BTNCLi','r0318i','FDX_NAD_NADP_OXi','CITt2','TRPS2','GUAt2','INSt2i','PRKINi','G16BPSi','G3PD8i','PROt2','MTHFCi','MNLt6i','FXXRDOi','OOR2','POR4i','ICDHxi','ICDHxi','GCALDD','AKGt2','TRPS2','FACOAL180i','CITCAti','CITCAti','FRDO','GNOXmqi','GNOXuqi','SHSL4','TRDR','OCOAT1','ALAt2','DCMPDAi','NAt3','PIt6bi','LLEUD','L_ILE3MRi','TRSAR','THRAi','GLYATi','SUCD1i','POR4i','FDOXRi','ICDHxi','PYNP1','ASPKi','GLUt2','DURADi','G16BPSi','G1PPTi','PPIte','PPIte','MACPMT','PPC','PPC','PPC','PPC','CBMK','TRDR','AMPSO3OXi','GALt2_2i','HISt2','L_LACD','ALAt4','UCO2Li','RIBFLVt2','FDOXRi','G3PD5','PPCOAOci','G16BPSi','FACOAL140i','R5PATi','R5PATi','MCCC','FRUt1r','ALAt2','ALAt2','r2471','G16BPSi','ILEt2','VALt2','DHLPHEORi','G3PD5','NTMAORi','PIt7ir','THMt3i','PROt4','GLUORi','G1PPi','FDOXRi','FRDO','ASP4DCi','NZP_NRei','NFORGLUAHi','FDOXRi','ACGAMt2','OIVD1','CITCAti','4ABZt2','FUMt2','TARCGLYLi','SULRi','CITt15i','GLFRDOi','THMDt2','HXANt2','ETOHt2','GSNt2','FUCt2_1i','GALt4i','PHEt2','ACOAD2fi','TSULt2i','MALFADOi','AKGt2','5ASAt2','MALt2','MALt4','MALt4','MTHFR2rev','NPR','PROD3i','DGORi','GNOXmqi','GNOXuqi','LPCDHi','CITCAti'};
    repInd = unique(findRxnIDs(panModel,replacedRxns))';
    if repInd(1) == 0
        PanRxnPresence(repInd(2:end),:) = 1;
    else
        PanRxnPresence(repInd,:) = 1;
    end
end

indBio = find(ismember(panRxns,'biomassPan'));
PanRxnPresence(indBio,:) = 1;
end