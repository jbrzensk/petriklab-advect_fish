%%%% Fish 1-D vector of ocean only cells to 2-D matrix of full grid
function nu2D = sub_happy_1Dto2D(GRD,nu,param)

nu1D = nu;
nu1D(nu>0) = ones;
nu1D(nu<=0) = zeros;

nu2D = nan*ones(param.ni,param.nj);
nu2D(GRD.ID) = nu1D;

end