# Panera: Pan-genera metabolic reconstruction
MATLAB algorithm to reconstruct pan-genus metabolic model (PGMM) from the existing genome-scale metabolic models (GSMMs).

## Motivation
A major limitation of requiring species information to leverage GSMMs for _in-silico_ community modelling exists while using 16S rRNA sequencing data due to their inherent ambiguity of lower taxonomic resolution. Most of the available tools frequently restrict taxonomic input to either genus or species level. In this study, we propose a unique framework to address two critical challenges in microbial community characterization: (i) uncertainty associated with taxonomic assignment in amplicon sequencing and (ii) scarcity of representative genus-level models. We introduce _**‘Panera’**_, an innovative framework designed to model microbial communities under this uncertainty and yet perform metabolic inferences using pan-genus metabolic models (PGMMs).

## Instruction to the users:
Our method allows the user to tailor the model using probabilities to their requirements. This proposed flexibility caters to two user types:  

- **Users with prior species information:** The potential species configuration for a specific context can be calculated from the previous studies and can be incorporated into the genus models using `customPanModel.m` for _in-silico_ microbial community construction in the [manuscript](https://doi.org/10.1016/j.isci.2024.110358)
  
- **Users exploring the metabolic potential of a genus:** The method enables the optimisation of species configuration within the genus model for the desired metabolic output. In addition, this approach facilitates the exploration of functional landscape of the whole genus. The adaptability of PGMM allows the users to investigate the extreme species combinations and their functionalities. Monitoring the resulting variation in metabolite production under this synthetic species combination could provide deeper understanding of the genus metabolic capabilities and potentially reveal the hidden metabolic niches in the genus.  

In case the user does not have any established information and is required to represent the genus, employing equal probabilities for all the species within a genus can be a valid option.

## Scripts
 - **Scripts required for building PGMMs**
     - `PanGenusModelReconstruction.m`  - Process all the inputs required to build PGMMs
     - `createPanGenusModel.m` - Builds PGMM for a given genera
     - `customPanModel.m` - Personalizes the PGMMs using species probability vector
     - `tutorial_PGMMreconstruction.m` - Quick tutorial on building PGMM for _Escherichia_ and customize the model with random vector
       
 - **Data and script to reproduce the tutorial**
     - `tutorial_script.m` - Script including the abundance data processing, PGMM construction and their application in _in-silico_ microbial community generation
     - `data` folder contains all the necessary abundance files and diet files for the simulation.
       
    Running `tutorial_script.m` in MATLAB will help in understanding the PGMM reconstruction and their applicability in hybrid community modelling. 
  
 - The folder `dependency` contains all the resource files required for the model building using Panera

## Prerequisites
  - MATLAB v2020b or later
  - COBRA Toolbox
  - IBM CPLEX v12.8 or later
  - [AGORA](https://www.vmh.life/files/reconstructions/AGORA/1.03/AGORA-1.03-With-Mucins.zip) / [AGORA2](https://www.vmh.life/files/reconstructions/AGORA2) model repositories


<details>
<summary> <b> Formulation </b> </summary>
Reconstruction of PGMM from species-specific GSMMs of a selected genus can be performed using the ‘Panera’ algorithm. The reconstruction pipeline employs three steps to produce a flexible PGMM: (i) Building a unified model from the reactions in all the species of a genus, (ii) Formulating biomass to represent all the species in a genus model, and (iii) Adding fields to accommodate the variation in species composition. The steps included in the PGMM reconstruction are detailed in this section.  

### Building a unified model from all the species GSMM of a genus 
- A database of all metabolites and reactions in Virtual Metabolic Human (VMH) models is retrieved from the Demeter pipeline 63. A separate database for the biomass reactions and metabolites of the species models is generated for the reconstruction (Table S1: Information of the species biomass reactions used in the model reconstruction).
- Reactions from the selected species GSMM models of a specific genus are extracted, and unique reactions (set of all the reactions) are identified to build a model.
- Unique reactions, except species biomass reactions, are integrated into a model using rBioNet. The fields such as rxnNames (reaction names), grRules (gene reaction association), compNames (Compartment where the reaction takes place - cytosol or Extracellular) and subsystems are added using a reaction and metabolite database. 

### Formulating biomass to represent the species in a genus model 

- The biomass reaction for the pan-genera model is formulated as the linear combination of biomass reactions of individual species in the genus: 

$$ v_{panBiomass} =  \sum_{i = 1}^{n} v_{bio}^{i}*s_{i}   $$

where $v_{panBiomass}$ is the biomass flux of the pan-genera model (Objective function), n is the number of species in the genus, 
$v_{bio}^{i}$ is the biomass flux of the $i^{th}$ species and $s_{i}$ is the coefficient for $i^{th}$ species, which implies the relative abundance or proportion of the microbial species in a community. The $s_{i}$ values can be adjusted to study the influence of a particular species in a genus. The reactions and metabolites associated with the ‘panBiomass’ and species biomass reactions are incorporated using biomass reaction and metabolite database. The default values of coefficients of species biomass ($s_{i}$) will be set to $\frac{1}{n}$. The default setting establishes an equal contribution from each species, and the coefficients can be adjusted to explore the distinct impact of a species.  
- Duplicate reactions or metabolites and reactions/metabolites involved in futile cycles are removed from the PGMM if the removal does not impact the growth of the model.
- The refined pan-genus model is examined for growth by optimising the model with biomass as an objective while constraining to a provided media condition. 

### Adding fields to accommodate the species composition variation 

- After PGMM refinement, a “reaction-species matrix”, a binary matrix representing whether the reaction is present (1) or absent (0) for an individual species, is combined as a field (‘rxnPresenceMat’) with the model.
- An ‘spList‘ field is incorporated into the model. Both ‘rxn-species matrix’ and ‘spList’ along with normalised ‘species probability vector’ will help filter the reactions to include in PGMM.

PGMM can be customised for a user-defined species composition using two key variables: (i) ‘species probability vector’, a vector of length n, a user-defined vector to reflect the estimated abundances of species in a community; and (ii) ‘rxn-species matrix’, a predefined matrix that encodes the reaction presence within a species. The product of these two variables determines whether the reaction is active in the model. A non-zero product indicates that the corresponding reaction is present in at least one species with a non-zero abundance, allowing it to be active within the model. Furthermore, the species probability vector plays a crucial role in incorporating compositional constraints into the biomass formulation. This formulation, in turn, influences the flux of internal and exchange reactions within the model.  
</details>


