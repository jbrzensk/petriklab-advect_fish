%% Output visualization
close all; clear all; clc;

% Set file to read in
filename = 'test_output_D_happy_D_0.1_maskLap_600k.mat';

% output image will have same name as input file
[~, name, ~] = fileparts(filename);
picfilename = [name '.png'];


load(filename);

title_str = 'Diffusion 0.1 test, Happy, day 30';

day = 30;


Sf = S_Sml_f(:,day);
Sp = S_Sml_p(:,day);
Sd = S_Sml_d(:,day);
Mf = S_Med_f(:,day);
Mp = S_Med_p(:,day);
Md = S_Med_d(:,day);
Lp = S_Lrg_p(:,day);
Ld = S_Lrg_d(:,day);

preySf = sub_1Dto2D(GRD1,Sf,param);
preySp = sub_1Dto2D(GRD1,Sp,param);
preySd = sub_1Dto2D(GRD1,Sd,param);
preyMf = sub_1Dto2D(GRD1,Mf,param);
preyMp = sub_1Dto2D(GRD1,Mp,param);
preyMd = sub_1Dto2D(GRD1,Md,param);
preyLp = sub_1Dto2D(GRD1,Lp,param);
preyLd = sub_1Dto2D(GRD1,Ld,param);

fignum = 2;

plotNineFish( fignum, title_str, Sf, Sp, Sd, Mf, Mp, Md, Lp, Ld, GRD1, param)
print(gcf, picfilename, '-dpng', '-r400');
