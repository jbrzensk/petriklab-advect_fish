%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
move_Y1_enc 	= true;
move_Y1_ingest 	= true;
move_Y1_mort 	= true;
move_Y1_nu 		= true;
move_Y1_prey 	= true;
move_spin50_enc 	= false;
move_spin50_ingest 	= false;
move_spin50_mort 	= false;
move_spin50_nu 		= false;
move_spin50_prey 	= false;
move_core_enc = false;
move_core_ingest = false;
move_core_mort = false;
move_core_nu = false;
move_core_prey = false;
spin50 = false;
core = true;


tic
	
if spin50
	Spinup_fished_gfdl_nomove_newdt
end
if core
	Historic_fished_gfdl_core_newdt
end	
	

if move_Y1_enc
	Spinup_fished_gfdl_move_newdt_enc_Y1
end
if move_Y1_ingest
	Spinup_fished_gfdl_move_newdt_ingest_Y1
end			
if move_Y1_mort
	Spinup_fished_gfdl_move_newdt_mort_Y1
end
if move_Y1_nu
	Spinup_fished_gfdl_move_newdt_nu_Y1
end
if move_Y1_prey
	Spinup_fished_gfdl_move_newdt_prey_Y1
end
	
if move_spin50_enc
	Spinup_fished_gfdl_move_newdt_enc
end
if move_spin50_ingest
	Spinup_fished_gfdl_move_newdt_ingest
end			
if move_spin50_mort
	Spinup_fished_gfdl_move_newdt_mort
end
if move_spin50_nu
	Spinup_fished_gfdl_move_newdt_nu
end
if move_spin50_prey
	Spinup_fished_gfdl_move_newdt_prey
end
													
if move_core_enc
	Historic_fished_gfdl_move_newdt_enc
end
if move_core_ingest
	Historic_fished_gfdl_move_newdt_ingest
end			
if move_core_mort
	Historic_fished_gfdl_move_newdt_mort
end
if move_core_nu
	Historic_fished_gfdl_move_newdt_nu
end
if move_core_prey
	Historic_fished_gfdl_move_newdt_prey
end										
											
toc
