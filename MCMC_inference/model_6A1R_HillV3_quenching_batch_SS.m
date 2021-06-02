function ss = model_6A1R_HillV3_quenching_batch_SS(params,data,local)
% sum-of-square function for the 6A1R_HillV3_compete model

nbatch = length(data);
% cumulate the sum-of-squares over the data sets (batches)

ss = 0;
for i=1:nbatch
  
  theta     = params(local==0|local==i);
%   tspan  = data{i}.ydata(:,1);
  ydata  = data{i}.ydata;
  TF = data{i}.xdata;
  params_fixed = data{i}.params_fixed;
  
  ymodel = model_6A1R_HillModel_V3_quenching_fixedBcdRNAPparams(theta, TF, params_fixed);
  
  ss = ss + sum((ydata - ymodel).^2);
  
end