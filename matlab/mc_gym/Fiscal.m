function fuzzylogicalculate = Fiscal(F_R, F_D ,F_E)
%	fis 
%  input: Reliability, Latency, Energy Consumption
% output: /omega value for trigger

fis_r = readfis('test_fuzzy3.fis');
fuzzylogicalculate = evalfis([F_R F_D F_E], fis_r);


end

