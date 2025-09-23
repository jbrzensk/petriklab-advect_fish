% Visualize output of FEISTY Historic-CORE at single locations
% 1988-2007, monthly means saved

clear 
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
load([vpath 'ocean_cobalt_grid.mat'],'geolon_t','geolat_t');
load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
load([vpath 'core_grid_360x200_id_locs_area_dep.mat'],'ids','abbrev','T');

spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';
longnames = T.Location;

nspots = length(spots);

%% CORE_Hindcast1988_no_move_All_fish03_newdt_locs.mat
exper2 = '1988_no_move_All_fish03_newdt';
load([fpath 'CORE_Hindcast',exper2,'_locs.mat'])

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

dSP = Mo_Sml_p;
dSF = Mo_Sml_f;
dSD = Mo_Sml_d;
dMP = Mo_Med_p;
dMF = Mo_Med_f;
dMD = Mo_Med_d;
dLP = Mo_Lrg_p;
dLD = Mo_Lrg_d;
dCO = Mo_Cobalt;

dF = dSF +dMF ;
dP = dSP +dMP +dLP ;
dD = dSD +dMD +dLD ;
    

%% CORE_Hindcast1988_no_move_All_fish03_locs.mat
exper1 = '1988_no_move_All_fish03_newdt';
load([fpath 'CORE_Hindcast',exper1,'_locs.mat'])

SP = Mo_Sml_p;
SF = Mo_Sml_f;
SD = Mo_Sml_d;
MP = Mo_Med_p;
MF = Mo_Med_f;
MD = Mo_Med_d;
LP = Mo_Lrg_p;
LD = Mo_Lrg_d;
CO = Mo_Cobalt;

F = SF + MF ;
P = SP + MP +LP ;
D = SD + MD +LD ;

t=1:size(SP,1);

%load([fpath 'CORE_Hindcast',exper,'_locs_lastyr_sum_mean_biom.mat']);
    
%%
for s=1:nspots
    loc = spots{s};
    lname = [loc '_'];
    loclong = longnames{s};
    
    %% TIME SERIES ----------------------------------------------------
    y=(t/12)+1988;
    figure(1)
    clf
    plot(y,log10(F(:,1,s)),'r','Linewidth',2); hold on;
    plot(y,log10(P(:,1,s)),'b','Linewidth',2); hold on;
    plot(y,log10(D(:,1,s)),'k','Linewidth',2); hold on;
    plot(y,log10(dF(:,1,s)),'--m','Linewidth',2); hold on;
    plot(y,log10(dP(:,1,s)),'--c','Linewidth',2); hold on;
    plot(y,log10(dD(:,1,s)),'--','Color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    legend('F','P','D','dF','dP','dD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-2 2])
    xlabel('Year')
    ylabel('log10 Biomass (g m^-^2)')
    title(['CORE ' loclong])
    stamp('')
    print('-dpng',[ppath 'CORE_Hindcast1988_no_move_All_fish03_newdt_olddt_ts_logmean_biomass_types_' loc '.png'])
    
    
    figure(2)
    clf
    plot(y,(dF(:,1,s) - F(:,1,s)),'r','Linewidth',2); hold on;
    plot(y,(dP(:,1,s) - P(:,1,s)),'b','Linewidth',2); hold on;
    plot(y,(dD(:,1,s) - D(:,1,s)),'k','Linewidth',2); hold on;
    legend('F','P','D')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-2 2])
    xlabel('Year')
    ylabel('Biomass (g m^-^2)')
    title(['CORE newdt - olddt ' loclong])
    stamp('')
    print('-dpng',[ppath 'CORE_Hindcast1988_no_move_All_fish03_newdt_olddt_ts_diff_mean_biomass_types_' loc '.png'])
    
    %  TIME SERIES ----------------------------------------------------
     
    
end



