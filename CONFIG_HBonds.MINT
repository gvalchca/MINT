Mode: TrajHbonds #Single or Traj or Download
out_name: example_Hbonds
working_dir: ./example/
file_name: ./example/example.pdb
files_dcd: example.dcd
chains_names:  #eg. A B C
nucleotides:
first_frame: 0
last_frame: -1
stride:1
cutoff:20
cutoff_stacking:10.5
vdw_cutoff_stacking: -0.5
OP_stacking_distance_cutoff:4.5
h_bond_atom:donor #donor/hydrogen
h_bond_l:3.5
h_bond_angle: 150
vmd:0
table_nucleotides:data/nucleotides_modified.csv
table_charges:data/charges_and_VDW_modified.csv
list_of_modified_nucs:data/unknown_modified.fa
force_field:AMBER
margin:0.8
time_cutoff:0.02
max_memory_GB:1.5
threads:4
only_analysis:0
create_csvs:1
pairing_cutoff: 0.1
