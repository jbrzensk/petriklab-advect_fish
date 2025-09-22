%%%%!! RUN SPINUP FOR ALL LOCATIONS
function Spinup_fished_gfdl_nomove_newdt()

% Add your specific subfunctions to the path
addpath(genpath('colleen_functions'));

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set time step
param.DThr = 6.0;

%! Set fishing rate
param.frate = 0.3;

param.dfrateF = nan;
param.dfrateP = nan;
param.dfrateD = nan;

%! Make core parameters/constants
param = make_param_newdt(param);

%! Grids
%vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
vpath = '/project/Feisty/GCM_Data/CORE-forced/';

%1-D
load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
GRD1 = GRD;
clear GRD

param.NX = length(GRD1.Z);
param.ID = 1:param.NX;
NX = length(GRD1.Z);
ID = 1:param.NX;

%! How long to run the model
YEARS = 200;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
opath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';
%opath = '/project/Feisty/NC/Matlab_new_size/';
exper = 'Spinup1988_nomove_dt6h_newdt';
[fname,simname,sname] = sub_fname_spin_core(param,opath,exper);

%! Storage variables
Year_Bent_bio = zeros(NX,DAYS);

Year_Sml_f = zeros(NX,DAYS);
Year_Sml_p = zeros(NX,DAYS);
Year_Sml_d = zeros(NX,DAYS);
Year_Med_f = zeros(NX,DAYS);
Year_Med_p = zeros(NX,DAYS);
Year_Med_d = zeros(NX,DAYS);
Year_Lrg_p = zeros(NX,DAYS);
Year_Lrg_d = zeros(NX,DAYS);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);

% % Last month of spinup without movement
% load([sname '_' simname '.mat']); 
% BENT.mass = BENT.bio;
% [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%! Dims of netcdf file
nt = 12 * YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% %%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_sml_f.nc'];
file_sml_p = [fname,'_sml_p.nc'];
file_sml_d = [fname,'_sml_d.nc'];
file_med_f = [fname,'_med_f.nc'];
file_med_p = [fname,'_med_p.nc'];
file_med_d = [fname,'_med_d.nc'];
file_lrg_p = [fname,'_lrg_p.nc'];
file_lrg_d = [fname,'_lrg_d.nc'];
file_bent  = [fname,'_bent.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');
ncidB  = netcdf.create(file_bent,'NC_WRITE');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~15 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,time_dim]);
% vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,time_dim]);
% vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidbioSD    = netcdf.defVar(ncidSD,'biomass','double',[xy_dim,time_dim]);
% vidprodSD   = netcdf.defVar(ncidSD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
% vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,time_dim]);
% vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,time_dim]);
% vidprodMD   = netcdf.defVar(ncidMD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
% vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
% vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

xy_dim     = netcdf.defDim(ncidB,'nid',NX);
time_dim   = netcdf.defDim(ncidB,'ntime',nt);
vidbioB    = netcdf.defVar(ncidB,'biomass','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
netcdf.endDef(ncidB);

%% %%%%%%%%%%%%%%%%%%%% Run the Model

addpath('matlab_functions');

load([vpath,'Data_ocean_cobalt_daily_1988.mat'],'COBALT');

MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    ti = num2str(YR)

    for DAY = param.DTday:param.DTday:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(DY,COBALT,GRD1,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

        %! Store
        Year_Bent_bio(:,DY) = BENT.mass;

        Year_Sml_f(:,DY) = Sml_f.bio;
        Year_Sml_p(:,DY) = Sml_p.bio;
        Year_Sml_d(:,DY) = Sml_d.bio;
        Year_Med_f(:,DY) = Med_f.bio;
        Year_Med_p(:,DY) = Med_p.bio;
        Year_Med_d(:,DY) = Med_d.bio;
        Year_Lrg_p(:,DY) = Lrg_p.bio;
        Year_Lrg_d(:,DY) = Lrg_d.bio;

    end %Days


    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker

        %! Put vars of netcdf file
        netcdf.putVar(ncidB,vidbioB,[0 MNT-1],[NX 1],mean(Year_Bent_bio(:,a(i):b(i)),2));
        netcdf.putVar(ncidB,vidTB,MNT-1,1,MNT);

        netcdf.putVar(ncidSF,vidbioSF,[0 MNT-1],[NX 1],mean(Year_Sml_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidbioSP,[0 MNT-1],[NX 1],mean(Year_Sml_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidbioSD,[0 MNT-1],[NX 1],mean(Year_Sml_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidbioMF,[0 MNT-1],[NX 1],mean(Year_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidbioMP,[0 MNT-1],[NX 1],mean(Year_Med_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidbioMD,[0 MNT-1],[NX 1],mean(Year_Med_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidbioLP,[0 MNT-1],[NX 1],mean(Year_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidbioLD,[0 MNT-1],[NX 1],mean(Year_Lrg_d(:,a(i):b(i)),2));

    end %Monthly mean

end %Years

% save([fname,'_Y1.mat'],'Year_Bent_bio','Year_Sml_f','Year_Sml_p','Year_Sml_d','Year_Med_f',...
%     'Year_Med_p','Year_Med_d','Year_Lrg_p','Year_Lrg_d','GRD1','GRD2','exper');

%! Close save
netcdf.close(ncidSF);
netcdf.close(ncidSP);
netcdf.close(ncidSD);
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);
netcdf.close(ncidB);



end
