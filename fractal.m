clear all;clc

global border_length 

 %% Calculate Molecular Fractal Dimension of Structures
i1 = 0;
total = 1+max(t)*1e4/1e6;
atmass = [1 12.011; 2 12.011; 3 12.011; 4 12.011; 5 12.011; 6 12.011; 7 15.999; 8 12.011; 9 12.011; 10 12.011; 11 15.999; 12 1.008; 13 15.999; 14 32.06; 15 12.011; 16 12.011; 17 12.011; 18 12.011; 19 14.007; 20 1.008; 21 1.008];
hydrogen = [12 20 21]; % Hydrogen Atom
H_bonds = [2,6,8,12,15,17,22,24,26,27];
bond_data_all = importdata("bond.txt");
bond_data = bond_data_all(~any(ismember(bond_data_all(:,2),H_bonds),2),:);
bond_atom_1 = bond_data(:,3);
bond_atom_2 = bond_data(:,4);

for i2 = 0:1e6:1e8
    i1 = i1+1;

    dump_name = "dump."+num2str(i2)+".txt";
    dump_read = importdata(dump_name,' ',9);
    dump_data_full = dump_read.data;
    mol_index = unique(dump_data_full(:,1));
    num_mol = length(mol_index);
    dump_data = dump_data_full(~ismember(dump_data_full(:,2),hydrogen),:); % Remove hydrogen from xyz data
    
    %% Replace the Atom Type with Atom Mass
    for iall = 1:1:size(dump_data,1)
        index4 = find(atmass(:,1)==dump_data(iall,2));
        dump_data(iall,2)=atmass(index4,2);
    end

    fid_xyz = fopen(dump_name);
    for line_no = 1:5
        a = fgetl(fid_xyz);
    end
    clear a
    boundary = str2num(fgetl(fid_xyz));
    boxMin = min(boundary);
    boxMax = max(boundary);
    border_length = boxMax-boxMin;

    for i3 = 1:1:num_mol
        clear current_xyz current_mol
        current_mol = [dump_data(dump_data(:,1)==mol_index(i3),9) dump_data(dump_data(:,1)==mol_index(i3),2) dump_data(dump_data(:,1)==mol_index(i3),3:5)];

%% Determine Bonding Correlations in Molecules
        
        clear bond_matrix
        bond_matrix = zeros(size(current_mol,1),5);
        bond_matrix(:,1) = current_mol(:,1);
        for i4 = 1:1:size(current_mol,1)
            index1 = find(bond_atom_1==current_mol(i4,1));
            index2 = find(bond_atom_2==current_mol(i4,1));
            i5 = 2;
            if ~isempty(index1)
                for itemp = 1:1:length(index1)
                    bond_matrix(i4,i5) = bond_atom_2(index1(itemp));
                    i5 = i5+1;
                end
            end
            if ~isempty(index2)
                for itemp = 1:1:length(index2)
                    bond_matrix(i4,i5) = bond_atom_1(index2(itemp));
                    i5 = i5+1;
                end
            end
        end
        clear itemp


%% Sort Molecule and Calculate Fractal Dimension
    atomdata = [current_mol(:,1) current_mol(:,3:5) bond_matrix(:,2:5)];
    clear current_xyz
    current_xyz = ConnectTree(atomdata);% XYZ data of current molecule, after connection
    df(i3,i1) = FractalDimension(current_xyz);



%% Calculate Radius of Gyration
    Rg(i3,i1) = RadiusGyration(current_xyz,current_mol);

    end

%% Output Progress Percentage
    %"Fractal Progress "+num2str((i1/total)*100)+"%"
end

df_ave = mean(df,2);
Rg_ave = mean (Rg,2);
Rgdf = Rg.^df;
Rgdf_ave = mean(Rgdf,2);

ave_structure = [df_ave,Rg_ave,Rgdf_ave];

%% Output
save result_fractal.txt -ascii df
save result_rg.txt -ascii Rg
save result_rgdf.txt -ascii Rgdf
save result_struc_ave.txt -ascii ave_structure

