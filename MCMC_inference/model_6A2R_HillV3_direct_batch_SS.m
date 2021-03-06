function ss = model_6A2R_HillV3_direct_batch_SS(params,data,local)
% sum-of-square function for the 6A1R_HillV3_direct model

nbatch = length(data);
% cumulate the sum-of-squares over the data sets (batches)

ss = 0;
for i=1:nbatch
  % theta is params filtered out the local parameters not associated with this
  % batch
  theta  = params(local==0|local==i); % either global or local to this batch.
  if i==1
      theta = params([1,3,4]);
  elseif i==2
      theta = params([1,2,4]);
  elseif i==3
      theta = params([1,2,3]);
  end
%   tspan  = data{i}.ydata(:,1);
  ydata  = data{i}.ydata;
  TF = data{i}.xdata;
  params_fixed = data{i}.params_fixed;
  
  ymodel = model_6A2R_HillModel_V3_direct_fixed_Kr_Bcd_RNAPparams(theta, TF, params_fixed);
  
  ss = ss + sum((ydata - ymodel).^2);
  
end