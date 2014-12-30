// Options test script
// you can change the model used by changing the .opt.optionModelApplied variable to either `Black76 or `BlackScholes

/ Load data
paramsTable:.opt.bls[`ParamsTable];
paramsTable,:("FFFFFFFFSS";" ")0: `options.txt;

/ Calculate the options
/ Asian options only working with futures with cost of carry = 0, i.e. yield = rate.
.opt.bls[`CalculateOptionPrice] each paramsTable
