%---------------------------------------------------------%

% Raghavan Thiagarajan, DanStem, Copenhagen, 24th November 2020.

% This script is used for reorganising the Tissue Analyser (TA) data. After extracting the "cell data" and "bond_data" from the TA,
% using these two data sheets are inputs, this script: (1) obtains the two adjacent cells corresponding to each bond; (2) gets their
% intensities; (3) normalises the bond intensity by the average of the adjacent cell intensities; (4) Then it places the bond intensities
% against the corresponding cell in the same format as in the "local_id_of_bonds" column in the "cell_data" sheet obtained from the TA.
% This script was mainly written to generate such a .csv file so that theorists working on Guilherme's project can use these
% intensities to map the forces in the tissue.

clear all
clc

% importing the bond data and the cell data from the excel sheets of Tissue Analyzer
bond_data = readtable('./bond_data.csv');
[bondrow,bondcol] = size(bond_data); % get the no of rows and columns in the imported table

cell_data = readtable('./cell_data.csv');
[cellrow,cellcol] = size(cell_data); % get the no of rows and columns in the imported table

% creating new columns for loading two data: (1) bond_id of the bonds that can be found in the 
% bond_data file; and (2) bond_intensities of the corresponding bonds.
cell_data.local_id_of_bonds_new(cellrow,1) = {nan}; % new column for bond_id 
cell_data.local_int_of_bonds_new(cellrow,1) = {nan}; % new column for bond intensities

% initialisation of start values for browsing the indices of bond and cell data
increment = 1;
bond_increment = 1;

% This while loop sets the range of the data in terms of rows in the data
% sheet of "cell_data"
while increment ~= cellrow
    
    % section 1 ---------------------------------------------------------%
    % The logic below makes sure we search and obtain the data only from a single
    % time point in the data sheet of "cell_data".
    cellrow_1 = increment;
    start_frame = cell_data.frame_nb(increment);
    cellrow_next = cell_data.frame_nb(increment + 1);
    count = increment;
    while (start_frame == cellrow_next && count < cellrow)
        count = count + 1;
        cellrow_next = cell_data.frame_nb(count);
    end
    end_frame = cell_data.frame_nb(count - 1);
    increment = count;
    cellrow_last = count - 1;
    
    % section 2 ---------------------------------------------------------%
    % The logic below makes sure we search and obtain the data only belonging to a single time point in the data sheet of "bond_data".
    % This is written in such a way that the same timepoint selected in the "cell_data" sheet is also
    % selected here.
    bond_count = bond_increment;
    while (bond_data.frame_nb(bond_count) == start_frame && bond_count < bondrow)
        bond_row_1 = bond_increment;
        bond_count = bond_count + 1;
    end
    bond_row_last = bond_count - 1;
    bond_increment = bond_count;
    
    % section 3 ---------------------------------------------------------%
    % The below for loop will go through every cell (every row) of the current time point fixed by the while loop above
    % in section 1. Look at the data sheet of "cell_data"
    for single_timepoint = cellrow_1 : cellrow_last
        
        % section 4 ---------------------------------------------------------%
        % below logic makes sure the "border cells" and the "border cells plus one"
        % cells are not included while reorganising the data.
        border_cell = cell_data.is_border_cell(single_timepoint);
        border_cell_plus_1 = cell_data.is_border_cell_plus_one(single_timepoint);
        if strcmp (border_cell{1,1}, 'FALSE') && strcmp(border_cell_plus_1{1,1}, 'FALSE')
            
            % section 5 ---------------------------------------------------------%
            % Below we fetch the bond numbers associated with the cells from the "cell_data" table. Since this is obtained in a
            % combined format and as character type, the bond numbers are first split and then converted to double type to be used as
            % numbers. Look at the  "local_id_of_bonds" column in the "cell_data" sheet.
            bond_id_set_0 = cell_data.local_id_of_bonds(single_timepoint);
            bond_id_set_1 = bond_id_set_0{1,1}(:,:);
            bond_id_set_2 = split(bond_id_set_1, "#");
            len_bond_id_set_2 = length(bond_id_set_2);
            
            % section 6 ---------------------------------------------------------%
            % We then go through every item in the list of obtained bond
            % numbers, and use these bond numbers to fetch their intensities and size from the "bond_data" sheet, and caluculate the mean
            % intensities. Then we fetch the cell numbers of the two neighboring cells (that are adjacent to the bond). We then
            % use these cell numbers to fetch the cytoplasmic intensities from the "cell_data" sheet. Then the bond intensity is
            % normalised by the cell cytoplasm intensity.
            for bondnos = 1 : len_bond_id_set_2
                
                bond_id_set_3 = bond_id_set_2{bondnos,1};
                bond_id = str2double(bond_id_set_3); % converting string to double
                
                bond_id_idx = find(bond_data.local_id_bonds(bond_row_1:bond_row_last,1) == bond_id); % getting the index of the corresponding bond id within the current time point.
                
                bondnos_for_array_generation = bondnos; % this value is to get the array size of the bond numbers belonging to a particular cell.
                
                % this if loop checks whether this bond_id can be found in the bond_data sheet. This check is done for the following
                % reason: sometimes, some erroneous bonds were introduced into the cell_data sheet by the TA but similar bonds
                % could not be identified in the bond_data sheet. The logic below makes sure that only bonds that are found in both
                % cell_data sheet and bond_data sheet are used for the analysis.
                if isempty(bond_id_idx)
                    % this value is to get the array size of the bond numbers belonging to a particular cell.
                    % here the value is changed if this for loop is skipped when the if condition is successfull.
                    bondnos_for_array_generation = bondnos_for_array_generation - 1;   
                    continue
                end
                
                bond_size = bond_data.bond_size_in_px(bond_id_idx); % bond size
                bond_int = bond_data.sum_px_int_vertices_excluded_12bits(bond_id_idx); % total intensity of the bond
                bond_int_mean = bond_int / bond_size; % mean intesnity of the bond
                
                cell_1_id = bond_data.cell_id_around_bond1(bond_id_idx); % getting the first adjacent cell
                % the following line is to be used only if the input of the "cell numbers" are in cell array format
                cell_1_id = cell2mat(cell_1_id); cell_1_id = str2double(cell_1_id); 
                cell_1_idx = find(cell_data.local_id_cells(cellrow_1:cellrow_last,1) == cell_1_id); % getting the index of that cell
                cell_1_area = cell_data.area_cells(cell_1_idx); % cell area
                cell_1_int = cell_data.sum_px_intensity_cells_12_bits(cell_1_idx); % total intensity of the cell
                cell_1_int_mean = cell_1_int / cell_1_area; % mean intensity of first cell
                
                
                cell_2_id = bond_data.cell_id_around_bond2(bond_id_idx); % getting the second adjacent cell
                % the following line is to be used only if the input of the "cell numbers" are in cell array format
                cell_2_id = cell2mat(cell_2_id); cell_2_id = str2double(cell_2_id);
                cell_2_idx = find(cell_data.local_id_cells(cellrow_1:cellrow_last,1) == cell_2_id); % getting the index of that cell
                cell_2_area = cell_data.area_cells(cell_2_idx); % cell area
                cell_2_int = cell_data.sum_px_intensity_cells_12_bits(cell_2_idx); % total intensity of the cell
                cell_2_int_mean = cell_2_int / cell_2_area; % mean intensity of second cell
                
                % normalising the bond intensity with the average of the two neighboring cell intensities.
                bond_int_norm(bondnos_for_array_generation,1) = bond_int_mean / ((cell_1_int_mean + cell_2_int_mean)/2);
                bond_id_array(bondnos_for_array_generation,1) = bond_id;
                
                % the following is done to create bond id in the format "bond_id#bond_id#bond_id" 
                bond_int_norm_conv_1 = num2str(bond_int_norm(:,1)');
                bond_int_norm_conv_2 = split(bond_int_norm_conv_1);
                bond_int_for_writing = join(bond_int_norm_conv_2,"#");
                % the following is done to create bond intensities in the format "bond_intensity#bond_intensity#bond_intensity"
                bond_id_array_conv_1 = num2str(bond_id_array(:,1)');
                bond_id_array_conv_2 = split(bond_id_array_conv_1);
                bond_id_for_writing = join(bond_id_array_conv_2,"#");
                
            end
            % the following arrays are cleared just to avoid errors.
            clear bond_int_norm; clear bond_id_array;
            % here we write the bond id and the bond intensities into the new columns
            cell_data.local_id_of_bonds_new(single_timepoint) = bond_id_for_writing;
            cell_data.local_int_of_bonds_new(single_timepoint) = bond_int_for_writing;         
            
        else
            % whenever the bond_ids of the border cells and border plus one cells are omitted, we put a zero to indicate that these cells
            % should not be considered.
            cell_data.local_id_of_bonds_new(single_timepoint) = {000};
            cell_data.local_int_of_bonds_new(single_timepoint) = {000};
        end
        
    end
    
end

% saving the table as .csv file. This file has two additional columns: (1)
% those bonds that are present in the bond_data file (sometimes bonds were missing; if bonds are missing in the 
% bond_data file, these bonds are not considered); (2) normalised bond intensities of those bonds.
writetable(cell_data,'cell_data_modified.csv');
clear all;
disp('Finished !');




