%%% Model building Scipt
%%% 11/19/2020

% Script contibutions by Dr. Douglas Chung, Dr. Tongli Zhang, Dr. Robert
% Sheehan
% © 2020 Rosa & Co. LLC

%%% Tongli Zhang

%%% Note%%%%%
%%% this script was adopted from previous ones by Doug
%%% As of 11/19/2020, this has been only used once on TGR019
%%% hence, this is an 'alpha' version, any comments/suggestions welcomed

% Updated 3/3/2021 by Renee Myers
% Updated version includes notes for generated parameters; saves an updated 
% SB project file with notes and parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% How to Use this Script %%%%%%%%%%%%%%%
% % % % % % 1. load the model project;
% % % % % % 2. Run the script to generate a file 'Before_Manually_modification_Reaction_Component.xlsx', which has a row for 'Modifier'
% % % % % % 3. Manaually choose the 'Activator' and 'Inhibitor', and put them in corresponding rows;
% % % % % % 4. Manually add reaction names in the 'Name' column;
% % % % % % 5. Save the manually modified file as 'After_Manually_modification_Reaction_Component.xlsx';
% % % % % % 6. Run the script again to complete equation writing
% % % % % % 7. Manually check the equations and units (REMOVE all User defined units)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% with Simbiology diagram as input;
%%%%% to output differential equations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Goal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Reduce the time and potential error associated with writing
%%%%% equations manually;
%%%%% make it easier to follow Rosa Standard in naming

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TGR019 project as a test case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Generate 'Before_Manually_modification_Reaction_Component.xlsx' spreadsheet
% Preamble
clc; close all; clear;

% Import SB project as 'model'
projfilename = '20201007 SFS001 SLE PhysioPD 2018b v1_07_testEqns.sbproj';
models = sbioloadproject(projfilename);
model = models.m1;

% Only choose active reactions to analyze, the inactive ones are for notes only %%%%
%%% notation only %%%%
Active_Status = cell2mat(get(model.Reactions, 'Active'));
ActiveReactions=model.Reactions(Active_Status);

% Set basic production and clearance reactions
Pro_Cle_Reaction=[]; % Store all production/clearance reactions

for nn=1:length(ActiveReactions)
    
    % Get reactants and products for indexed reaction
    Reactants_of_reaction=ActiveReactions(nn).Reactants;
    Product_of_reaction=ActiveReactions(nn).Products;
    reactionObj = ActiveReactions(nn);
    
    % Account for null -> products (production)
    if isempty(Reactants_of_reaction)  %%% produced from null 
        if ~strcmp(reactionObj.Name,[reactionObj.product.Name,'_production'])
            reactionObj.Name=[reactionObj.product.Name,'_production'];
        end
        
        % Set reaction object note 
        if ~strcmp(reactionObj.Note,[reactionObj.product.Name,'_production'])
            reactionObj.Note=[reactionObj.product.Name,'_production'];
        end
        
        % Set parameter object for null -> products as zero order production
        parameterObj1 = addparameter (model, sprintf('%s_rate_k', reactionObj.Name), 1);
        set (reactionObj, 'ReactionRate', sprintf('%s', parameterObj1.Name)); %%% zero order production

        Pro_Cle_Reaction=[Pro_Cle_Reaction,nn];
    end
    
    % Account for reactants -> null (clearance)
    if isempty(Product_of_reaction) %%% into null
        
        if ~strcmp(reactionObj.Name,[reactionObj.Reactants.Name,'_clearance'])
            reactionObj.Name=[reactionObj.Reactants.Name,'_clearance'];
        end
        
        % Set reaction object note 
        if ~strcmp(reactionObj.Note,[reactionObj.Reactants.Name,'_clearance'])
            reactionObj.Note=[reactionObj.Reactants.Name,'_clearance'];
        end
        
        % Set parameter object for reactants -> null as first order clearance
        parameterObj1 = addparameter (model, sprintf('%s_rate_k', ActiveReactions(nn).Name), 1);
        set (reactionObj, 'ReactionRate', sprintf('%s*%s', parameterObj1.Name, reactionObj.Reactants.Name)); % first order clearance
        Pro_Cle_Reaction=[Pro_Cle_Reaction,nn];
                
        % Set parameter note for reactants -> null (clearance)
        set(parameterObj1, 'Note', sprintf('Clearance rate constant for %s', reactionObj.Reactants.Name));
    end  
end

% Remove the Basic production and clearance reactions
Reaction_Indexes=1:length(ActiveReactions);
Reaction_Index_of_interest = ~ismember(Reaction_Indexes, Pro_Cle_Reaction); % Find indexes for all valid reactions
Remained_Reactions = ActiveReactions(Reaction_Index_of_interest);

%Initialize filename and table to track reaction components
%Using tables instead of array to facilitate use of writetable command
%instead of writematrix to insure backwards compatibility

% Add reactions to table
for mm=1:length(Remained_Reactions) 
    % Get all reactions/reactants/products
    reaction = Remained_Reactions(mm);
    reac = reaction.reactants;
    prod = reaction.products;
    
    % Store names of all reactant/products in cell array (R2020 only?)
    % reacspec = reshape({reac.name},length(reac),1);
    % prodspec = reshape({prod.name},length(prod),1);

    % Store names of all reactants/products in cell array
    reacspec = cell(length(reac),1);
    prodspec = cell(length(prod),1);
    for idx = 1:length(reac)
        spec = reac(idx);
        reacspec{idx} = spec.name;
    end
    for jdx = 1:length(prod)
        spec = prod(jdx);
        prodspec{jdx} = spec.name;
    end
    
    % Get reactants, modifiers, and products
    reactant{mm} = setdiff(reacspec,prodspec);
    modifiers{mm} = intersect(reacspec,prodspec); %modifiers are both reactants and products
    product{mm} = setdiff(prodspec,reacspec);
    newRow{mm} = {num2str(mm)}; % Index label for table
end

% Make table to track reactants, products, and modifiers
reactionTable = cell2table([newRow', reactant',product',modifiers']); % Transpose into column vectors
reactionTable.Properties.VariableNames = {'Reaction_Index','Reactant','Product','Modifier'};

% Turn multiple reactants/modifiers/products into comma separated list and extract nested cell data
for nn = 2:width(reactionTable) % Loop over reactants, products, modifiers
    for mm = 1:height(reactionTable) % Loop over each row of table
        if iscell(reactionTable{mm,nn}{:}) && length(reactionTable{mm,nn}{:}) > 1
            reactionTable{mm,nn}{:} = strjoin(reactionTable{mm,nn}{:}, ', ');
        elseif isempty(reactionTable{mm,nn}{:})
            reactionTable(mm,nn) = {''};
        elseif iscell(reactionTable{mm,nn}{:})
             reactionTable{mm,nn}{:} = reactionTable{mm,nn}{:}{:};
        end
    end
end

%Add table columns where user will define which modifiers are activators
%and inhibitors
newColumns = cell2table(cell(size(reactionTable,1),3));
newColumns.Properties.VariableNames = {'Activator','Inhibitor','Name'};
reactionTable = [reactionTable newColumns];

% Output reaction components to an Excel File
filename='Reaction_Component.xlsx';

%Write reaction table to file
if verLessThan('matlab','9.8') % Uses xlswrite for versions before R2020a
    xlswrite(filename, [reactionTable.Properties.VariableNames; table2array(reactionTable)],'Sheet',1)
else % Use writetable for versions 2020a and beyond
    writetable(reactionTable, filename,'Sheet',1);
end
movefile 'Reaction_Component.xlsx' 'Before_Manually_modification_Reaction_Component.xlsx'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Note: this file needs manual modification
%%% 1). separate activator and inhibitor columns
%%% 2). add reaction names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load the manually modified file to write equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Part 2: Generate final spreadsheet

% Load manually modified file to write equations
Filename = 'After_Manually_modification_Reaction_Component.xlsx';
File_To_Use=readtable(Filename);

% Initialize table to store file data
T=File_To_Use;
T.Reaction_Rates = cell(size(T,1),1);

% Error check file to make sure it contains all necessary data
%--------------------------------------------------------------------------
check_names = T.Name;
if any(cellfun(@isempty,check_names)) % Check if any name entries are blank
    warndlg('Make sure input spreadsheet contains entries for all reaction names.')
    return
end
if length(check_names) ~= length(unique(check_names)) % Check if any name entries are repeated
    warndlg('Make sure input spreadsheet does not contain repeated reaction names.')
    return
end
if any(contains(check_names, ' ')) % Check if any name entries contain spaces
    warndlg('Make sure input spreadsheet does not contain spaces in reaction names. Spaces will be removed.')
    T.Name = strrep(check_names, ' ', '_'); % Replace any spaces with underscores
end
% Replace any NaN values in table with empty cell arrays
if ~isa(File_To_Use.Activator, 'cell')
    File_To_Use.Activator = num2cell(File_To_Use.Activator);
    File_To_Use.Activator(cellfun(@isnan,File_To_Use.Activator))={''};
end

if ~isa(File_To_Use.Inhibitor, 'cell')
    File_To_Use.Inhibitor = num2cell(File_To_Use.Inhibitor);
    File_To_Use.Inhibitor(cellfun(@isnan,File_To_Use.Inhibitor))={''};
end
%--------------------------------------------------------------------------

% Loop over all remained reactions
for nn=1:length(Remained_Reactions)
    
    reactionObj = Remained_Reactions(nn);
    
    reaction_Name = char(File_To_Use.Name(nn)); %Leads to error if no name or repeat name given.
                                                %Name is also incorporated into parameter names. Means no spaces.
    Activators_All = File_To_Use.Activator(nn);
    Inhibitors_All = File_To_Use.Inhibitor(nn); 
    reactantName = File_To_Use.Reactant(nn);
    
    root=reaction_Name;
    simbiomodel=model;
    reactantName = File_To_Use.Reactant(nn);

    parameterObj1 = addparameter (model, sprintf('%s_rate_k', reaction_Name), 1);
    
    % Handle logic for tracking modifier lists in all three cases of multiple modifiers: 
        % (1) activators and inhibitors, (2) activators only, (3) inhibitors only
        
     if numel(Activators_All{1}) > 0 
        if numel(Inhibitors_All{1}) > 0 % Check if there are activator and inhibitors
            equationType = 'activation_and_inhibition';
        else % Check if there are activators only
           equationType = 'activation_equation';
        end
    elseif numel(Inhibitors_All{1}) > 0 % Check if there are inhibitors only
        equationType = 'inhibition_equation';
    end
        
    %Handle logic for converting modifier lists to eqns in all three cases of mutliple modifiers:
    %activators and inhibitors, activators only, inhibitors only
    %Ensures Activators use Emax/EC50 and are additive, Inhibitors use
    %Imax/IC50 and are multiplicative. Correctly tracks extra parentheses
    %for inhibition reactions.

    switch equationType
        case 'activation_and_inhibition' % Activators and inhibitors
           % Write equation
           if ~isempty(char(reactantName))
                activation_and_inhibition_equation = sprintf('%s*%s*(1', parameterObj1.Name,char(reactantName)); %initialize
           else
                activation_and_inhibition_equation = sprintf('%s*(1', parameterObj1.Name); %initialize
           end  
           for idx = 1:length(Activators_All)
                modifier = Activators_All{idx};
                activation = strcat(modifier,'_',root,'_Emax*',modifier,'/(',modifier,'_',root,'_EC50+',modifier,')'); %activation term for each modifier
                activation_and_inhibition_equation = strcat(activation_and_inhibition_equation,'+',activation); %(1 + SUMOF activators)
                newparam1 = addparameter(simbiomodel,strcat(modifier,'_',root,'_Emax')); %create parameters with default value
                newparam2 = addparameter(simbiomodel,strcat(modifier,'_',root,'_EC50'));
                
                % Generate notes for each parameter
                newparam1.Note = sprintf('Maximum effect of %s on %s',modifier, root);
                newparam2.Note = sprintf('Level of %s at which its effect on %s is half maximal',modifier, root);
           end
           activation_and_inhibition_equation = strcat(activation_and_inhibition_equation,')');
           for idx = 1:length(Inhibitors_All)
                modifier = Inhibitors_All{idx};
                inhibition = strcat('*(1-',modifier,'_',root,'_Imax*',modifier,'/(',modifier,'_',root,'_IC50+',modifier,')'); %activation term for each modifier
                activation_and_inhibition_equation = strcat(activation_and_inhibition_equation,inhibition); %(1 + SUMOF activators)
                newparam1 = addparameter(simbiomodel,strcat(modifier,'_',root,'_Imax')); %create parameters with default value
                newparam2 = addparameter(simbiomodel,strcat(modifier,'_',root,'_IC50'));
                
                % Generate notes for each parameter
                newparam1.Note = sprintf('Maximum fractional inhibition of %s by %s',root, modifier);
                newparam2.Note = sprintf('Level of %s at which its effect on %s is half maximal',modifier, root);
           end
           modifier_equation = strcat(activation_and_inhibition_equation,')');
           T.Reaction_Rates{nn}=modifier_equation;
           
        case 'activation_equation' % Activators only
            % Write equation
            if ~isempty(char(reactantName))
                activation_equation = sprintf('%s*%s*(1', parameterObj1.Name,char(reactantName)); %initialize
            else
                activation_equation = sprintf('%s*(1', parameterObj1.Name); %initialize
            end
            for idx = 1:length(Activators_All)
                modifier = Activators_All{idx};
                activation = strcat(modifier,'_',root,'_Emax*',modifier,'/(',modifier,'_',root,'_EC50+',modifier,')'); %activation term for each modifier
                activation_equation = strcat(activation_equation,'+',activation,')'); %(1 + SUMOF activators)
                newparam1 = addparameter(simbiomodel,strcat(modifier,'_',root,'_Emax')); %create parameters with default value
                newparam2 = addparameter(simbiomodel,strcat(modifier,'_',root,'_EC50'));
                
                % Generate notes for each parameter
                newparam1.Note = sprintf('Maximal effect of %s on %s',modifier, root);
                newparam2.Note = sprintf('Half maximal effective concentration of %s on %s',modifier, root);
                
            end
            modifier_equation = strcat(activation_equation,')');
            T.Reaction_Rates{nn}=modifier_equation;
        
        case 'inhibition_equation' % Inhibitors only
            % Write equation
            if ~isempty(char(reactantName))
                inhibition_equation = sprintf('%s*%s', parameterObj1.Name,char(reactantName)); %initialize
            else
                inhibition_equation = sprintf('%s', parameterObj1.Name); %initialize
            end
            for idx = 1:length(Inhibitors_All)
                modifier = Inhibitors_All{idx};
                inhibition = strcat('*(1-',modifier,'_',root,'_Imax*',modifier,'/(',modifier,'_',root,'_IC50+',modifier,')'); %activation term for each modifier
                inhibition_equation = strcat(inhibition_equation,inhibition); %(1 + SUMOF activators)
                addparameter(simbiomodel,strcat(modifier,'_',root,'_Imax')); %create parameters with default value
                addparameter(simbiomodel,strcat(modifier,'_',root,'_IC50'));
                
                % Generate notes for each parameter
                newparam1.Note = sprintf('Maximal inhibition capacity of %s on %s',modifier, root);
                newparam2.Note = sprintf('Half maximal inhibition concentration of %s on %s',modifier, root);
            end
            modifier_equation = strcat(inhibition_equation,')');
            T.Reaction_Rates{nn}=modifier_equation;  
    end
   
    if ~strcmp(reactionObj.Note,reaction_Name)
        reactionObj.Note=reaction_Name;
    end
    set(reactionObj, 'ReactionRate', modifier_equation);

    %clear reaction-specific variables before next loop iteration
    clear activation_and_inhibition_equation; clear activation_equation; clear inhibition_equation; clear equationType;
end

% Write down the final equations
if verLessThan('matlab','9.8') % Uses xlswrite for versions before R2020a
    xlswrite('Final_Rea_file.xlsx', [T.Properties.VariableNames; table2array(T)],'Sheet',1)
else
    writetable(T, 'Final_Rea_file.xlsx');
end

% Save an updated SB project file with notes and parameters
% New project file is saved as 'originalFileName_newparams.sbproj'
m1 = model;
newProjectFilename = [extractBefore(projfilename, '.sbproj') '_newparams.sbproj'];
sbiosaveproject(newProjectFilename, 'm1')

