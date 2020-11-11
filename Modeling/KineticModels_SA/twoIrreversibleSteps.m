function y = twoIrreversibleSteps(k1,k2,t,D1,D2)



y= ((k1+D1)+(-1).*(k2+D2)).^(-1).*((k1+D1)+(-1).*exp(1).^((-1).*(k2+D2).*t).*(k1+D1)+((-1)+ ...
  exp(1).^((-1).*(k1+D1).*t)).*(k2+D2));

end